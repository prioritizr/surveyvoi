#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_posterior_probability_matrix.h"
#include "rcpp_predict_missing_rij_data.h"
#include "rcpp_expected_value_of_action.h"

double expected_value_of_decision_given_survey_scheme(
  Eigen::MatrixXd &rij, // observed presence/absence matrix
  Eigen::MatrixXd &pij, // prior matrix
  Eigen::MatrixXd &wij, // site weight matrix
  std::vector<bool> &survey_features, // features that we want to survey, 0/1
  Eigen::VectorXd &survey_sensitivity,
  Eigen::VectorXd &survey_specificity,
  std::vector<bool> &pu_survey_solution, // planning units to survey, 0/1
  Eigen::VectorXd &pu_survey_status,  // has planning unit been surveyed? 0/1
  Eigen::VectorXd &pu_survey_costs,   // cost of surveying planning units
  Eigen::VectorXd &pu_purchase_costs, // cost of purchasing planning units
  Eigen::VectorXd &pu_purchase_locked_in,  // planning units that are locked in
  MatrixXfRM &pu_env_data, // environmental data
  std::vector<std::vector<std::string>> &xgb_parameter_names, // xgboost
                                                              // parameter names
  std::vector<std::vector<std::string>> &xgb_parameter_values, // xgboost
                                                               // parameter
                                                               // values
  std::vector<std::size_t> &n_xgb_nrounds,  // xgboost training rounds
  std::vector<std::vector<std::vector<std::size_t>>> &xgb_train_folds,
  std::vector<std::vector<std::vector<std::size_t>>> &xgb_test_folds,
  Eigen::VectorXd &obj_fun_preweight,  // objective function calculation term
  Eigen::VectorXd &obj_fun_postweight,  // objective function calculation term
  Eigen::VectorXd &obj_fun_target,  // objective function calculation term
  std::size_t n_approx_obj_fun_points, // number of approximate points
  double total_budget, // total budget for surveying + monitor costs
  double optim_gap    // optimality gap for prioritizations
) {
  // initialization
  /// constant variables
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = rij.rows();
  const std::size_t n_f_survey =
    std::accumulate(survey_features.begin(), survey_features.end(), 0);
  const std::size_t n_pu_surveyed_in_scheme =
    std::accumulate(pu_survey_solution.begin(), pu_survey_solution.end(), 0);
  const std::size_t n_pu_surveyed_already =
    static_cast<std::size_t>((pu_survey_status.array()).sum());

  /// integer over-flow checks, highest std::size_t value is 1e+18
  if (n_pu_surveyed_in_scheme > 20)
    Rcpp::stop("number of planning units in selected in survey scheme is too large (i.e. >20)");

  /// calculate number of outcomes for a given feature
  const std::size_t n_f_outcomes = n_states(n_pu_surveyed_in_scheme) + 1;

  /// calculate toal number of outcomes across all features
  mpz_class n_outcomes;
  n_states(n_pu_surveyed_in_scheme * n_f_survey, n_outcomes);
  n_outcomes = n_outcomes + 1; // increment to include final outcome

  // data integrity checks
  /// calculate remaining budget
  double remaining_budget = total_budget;
  for (std::size_t j = 0; j < n_pu; ++j)
    if (pu_survey_solution[j])
      remaining_budget -= pu_survey_costs[j];

  /// return zero if there is no remaining budget left after accounting
  /// for the survey costs
  if (std::abs(remaining_budget) < 1.0e-6)
    return 0.0;
  /// throw error if surveys cost more than the total budget
  double locked_in_cost =
    (pu_purchase_locked_in.array() * pu_purchase_costs.array()).sum();
  if ((remaining_budget - locked_in_cost) < 0) {
    Rcpp::stop("cost of surveying the planning units selected in scheme exceeds the total budget");
    return 0.0;
  }
  // create sparse representations of the boolean vectors
  //// store indices of planning units that are selected for surveying
  std::vector<std::size_t> pu_survey_solution_idx;
  pu_survey_solution_idx.reserve(n_pu_surveyed_in_scheme);
  for (std::size_t i = 0; i < n_pu; ++i)
    if (pu_survey_solution[i])
      pu_survey_solution_idx.push_back(i);

  //// store indices of planning units that have already been surveyed
  std::vector<std::size_t> pu_survey_status_idx;
  pu_survey_status_idx.reserve(n_pu_surveyed_already);
  for (std::size_t i = 0; i < n_pu; ++i)
    if (pu_survey_status[i] < 0.5)
      pu_survey_status_idx.push_back(i);

  //// store indices of planning units that need feature probs. predicted
  std::vector<std::size_t> pu_model_prediction_idx;
  pu_model_prediction_idx.reserve(n_pu);
  for (std::size_t i = 0; i < n_pu; ++i)
    if ((pu_survey_status[i] < 0.5) && (!pu_survey_solution[i]))
      pu_model_prediction_idx.push_back(i);
  pu_model_prediction_idx.shrink_to_fit();
  const std::size_t n_pu_model_prediction = pu_model_prediction_idx.size();

  //// store indices of features that need surveying
  std::vector<std::size_t> survey_features_idx;
  survey_features_idx.reserve(n_f);
  for (std::size_t i = 0; i < n_f; ++i)
    if (survey_features[i])
      survey_features_idx.push_back(i);
  survey_features_idx.shrink_to_fit();

  /// store indices for features in sparse format for reverse lookup
  /// e.g. if survey_features = [0, 1, 0, 1]
  ///         survey_features_idx = [1, 3]
  ///         survey_features_rev_idx = [?, 0, ?, 1], note ? = 0 as a default
  ///                                                 but could be anything
  ///                                                 since it is not used
  std::vector<std::size_t> survey_features_rev_idx(n_f, 0);
  {
    std::size_t counter = 0;
    for (std::size_t i = 0; i < n_f; ++i) {
      if (survey_features[i]) {
        survey_features_rev_idx[i] = counter;
        ++counter;
      }
    }
  }

  /// declare temporary variables used in the main loop
  double curr_value;
  double curr_expected_value_of_action_given_outcome;
  double curr_probability_of_outcome;
  double curr_expected_value_of_decision =
    std::numeric_limits<double>::infinity();
  Eigen::MatrixXd curr_state(n_pu_surveyed_in_scheme * n_f_survey, 1);
  Eigen::MatrixXd curr_oij = rij;
  Eigen::MatrixXd curr_mij = rij;
  Eigen::MatrixXd curr_pij(n_f, n_pu);
  curr_pij.setConstant(-100.0);
  Eigen::MatrixXd curr_total_probability_of_model_positive(n_f_survey, n_pu);
  Eigen::MatrixXd curr_total_probability_of_model_negative(n_f_survey, n_pu);
  std::vector<bool> curr_solution(n_pu);
  std::vector<std::size_t> feature_outcome_idx(n_f_survey);

  // preliminary processing
  /// create subset of prior matrix for just the species that need surveys
  Eigen::MatrixXd pij_survey_species_subset(n_f_survey, n_pu);
  for (std::size_t i = 0; i < n_f_survey; ++i)
    pij_survey_species_subset.row(i) = pij.row(survey_features_idx[i]);

  /// calculate the total probabilities of positive and negative outcomes
  /// from the surveys
  Eigen::MatrixXd total_probability_of_survey_positive(n_f, n_pu);
  Eigen::MatrixXd total_probability_of_survey_negative(n_f, n_pu);
  Eigen::MatrixXd total_probability_of_survey_positive_log;
  Eigen::MatrixXd total_probability_of_survey_negative_log;
  total_probability_of_positive_result(
    pij, survey_sensitivity, survey_specificity,
    total_probability_of_survey_positive);
  total_probability_of_negative_result(
    pij, survey_sensitivity, survey_specificity,
    total_probability_of_survey_negative);
  total_probability_of_survey_positive_log =
    total_probability_of_survey_positive;
  total_probability_of_survey_negative_log =
    total_probability_of_survey_negative;
  log_matrix(total_probability_of_survey_positive_log);
  log_matrix(total_probability_of_survey_negative_log);

  /// initialize prioritization object
  Prioritization prioritize(
    rij.cols(), rij.rows(), pu_purchase_costs, pu_purchase_locked_in,
    obj_fun_preweight, obj_fun_postweight, obj_fun_target,
    n_approx_obj_fun_points, remaining_budget, optim_gap);

  /// overwrite missing data for feature we are not interested in surveying
  /// using the prior data
  for (std::size_t j = 0; j < n_pu; ++j)
    if (!pu_survey_status[j])
      for (std::size_t i = 0; i < n_f; ++i)
        if (!survey_features[i])
          curr_oij(i, j) = pij(i, j);

  /// store indices for cells in the rij matrix that will be used for
  /// simulating different outcomes
  std::vector<std::size_t> rij_outcome_idx;
  {
    std::size_t k = 0;
    rij_outcome_idx.reserve(n_pu_surveyed_in_scheme * n_f_survey);
    for (std::size_t j = 0; j < n_pu; ++j) {
      for (std::size_t i = 0; i < n_f; ++i) {
        if (pu_survey_solution[j] && survey_features[i]) {
          rij_outcome_idx.push_back(k);
        }
        ++k;
      }
    }
  }

  /// fit distribution models for feature under different outcomes
  /// for each feature separately
  Eigen::Array<std::vector<BoosterHandle>, Eigen::Dynamic, Eigen::Dynamic>
    model_beta(n_f_survey, n_f_outcomes);
  Eigen::MatrixXd model_sensitivity(n_f_survey, n_f_outcomes);
  Eigen::MatrixXd model_specificity(n_f_survey, n_f_outcomes);
  Eigen::MatrixXd model_fit_oij = curr_oij;
  Eigen::MatrixXd curr_outcome(1, n_pu_surveyed_in_scheme);
  std::vector<std::size_t> curr_outcome_idx(n_pu_surveyed_in_scheme);
  std::iota(curr_outcome_idx.begin(), curr_outcome_idx.end(), 0);
  mpz_class tmp;
  for (std::size_t o = 0; o < n_f_outcomes; ++o) {
    /// generate the i'th outcome from surveying the planning units
    tmp = o;
    nth_state_sparse(tmp, curr_outcome_idx, curr_outcome);

    /// create a copy of the oij matrix with prior-filled data
    /// with all the planning units for features that need surveying
    /// updated
    for (std::size_t j = 0; j < n_pu_surveyed_in_scheme; ++j)
      for (std::size_t i = 0; i < n_f_survey; ++i)
        model_fit_oij(survey_features_idx[i], pu_survey_solution_idx[j]) =
          curr_outcome(j);

    /// fit model and store results
    fit_xgboost_models_and_assess_performance(
      model_fit_oij, wij, pu_env_data, survey_features_idx,
      o, xgb_parameter_names, xgb_parameter_values,
      n_xgb_nrounds, xgb_train_folds, xgb_test_folds,
      model_beta, model_sensitivity, model_specificity);
  }

  // validate calculations
  assert_valid_probability_data(model_sensitivity,
                                "issue calculating model sensitivites");
  assert_valid_probability_data(model_specificity,
                                "issue calculating model specificities");

  // main processing
  mpz_class o = 0;
  while (cmp(o, n_outcomes) < 0) {
    /// generate the o'th outcome from surveying the planning units across
    /// all species
    nth_state_sparse(o, rij_outcome_idx, curr_oij);

    // find out which model coefficients should be used for making predictions,
    // in other words, find out what outcome has been generated for the feature.
    // this involves finding the cell numbers of rij that correspond to
    // the planning units that have been selected for surveying
    which_feature_state(curr_oij, survey_features_idx, pu_survey_solution_idx,
                        feature_outcome_idx);

    /// generate modelled predictions for species we are interested in surveying
    predict_missing_rij_data(curr_oij, pu_env_data, survey_features_idx,
                             feature_outcome_idx, pu_model_prediction_idx,
                             model_beta);
    assert_valid_probability_data(curr_oij, "issue predicting missing data");

    /// calculate total probability of models' positive results
    total_probability_of_positive_model_result(
      pij_survey_species_subset, model_sensitivity, model_specificity,
      feature_outcome_idx, curr_total_probability_of_model_positive);
    assert_valid_probability_data(curr_total_probability_of_model_positive,
                                  "issue calculating total model positives");

    /// calculate total probability of models' negative results
    total_probability_of_negative_model_result(
      pij_survey_species_subset, model_sensitivity, model_specificity,
      feature_outcome_idx, curr_total_probability_of_model_negative);
    assert_valid_probability_data(curr_total_probability_of_model_negative,
                                  "issue calculating total model negatives");

    /// generate posterior data
    posterior_probability_matrix(
      rij, pij, curr_oij,
      pu_survey_solution,
      survey_features,
      feature_outcome_idx, survey_features_rev_idx,
      survey_sensitivity, survey_specificity,
      total_probability_of_survey_positive,
      total_probability_of_survey_negative,
      model_sensitivity, model_specificity,
      curr_total_probability_of_model_positive,
      curr_total_probability_of_model_negative,
      curr_pij);

    assert_valid_probability_data(curr_pij,
                                  "issue calculating posterior probabilities");

    /// generate prioritisation
    prioritize.add_rij_data(curr_pij);
    prioritize.solve();
    prioritize.get_solution(curr_solution);

    /// calculate expected value of the prioritisation
    curr_expected_value_of_action_given_outcome =
      std::log(expected_value_of_action(
        curr_solution, curr_pij, obj_fun_preweight, obj_fun_postweight,
        obj_fun_target));

    /// calculate likelihood of outcome
    curr_probability_of_outcome = log_probability_of_outcome(
      curr_oij, total_probability_of_survey_positive_log,
      total_probability_of_survey_negative_log, rij_outcome_idx);

    /// calculate expected value of action
    if (std::isinf(curr_expected_value_of_decision)) {
      curr_expected_value_of_decision =
        curr_expected_value_of_action_given_outcome +
        curr_probability_of_outcome;
    } else {
      curr_expected_value_of_decision =
        log_sum(curr_expected_value_of_decision,
                curr_expected_value_of_action_given_outcome +
                curr_probability_of_outcome);
    }

    /// reset oij matrix so that -1s are present for planning units/features
    /// that need surveying
    for (std::size_t j = 0; j < n_pu; ++j)
      if (!pu_survey_status[j])
        for (std::size_t i = 0; i < n_f; ++i)
          if (survey_features[i])
            curr_oij(i, j) = -1.0;
    /// increment o loop variable
    o = o + 1;
  }

  // clean-up
  for (std::size_t i = 0; i < (n_f_outcomes * n_f_survey); ++i)
    for (std::size_t k = 0; k < model_beta(i).size(); ++k)
      XGBoosterFree(model_beta(i)[k]);

  // exports
  return std::exp(curr_expected_value_of_decision);
}

// [[Rcpp::export]]
double rcpp_expected_value_of_decision_given_survey_scheme(
  Eigen::MatrixXd rij,
  Eigen::MatrixXd pij,
  Eigen::MatrixXd wij,
  std::vector<bool> survey_features,
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  std::vector<bool> pu_survey_solution,
  Eigen::VectorXd pu_survey_status,
  Eigen::VectorXd pu_survey_costs,
  Eigen::VectorXd pu_purchase_costs,
  Eigen::VectorXd pu_purchase_locked_in,
  Eigen::MatrixXf pu_env_data,
  Rcpp::List xgb_parameters,
  Rcpp::List xgb_train_folds,
  Rcpp::List xgb_test_folds,
  std::vector<std::size_t> n_xgb_nrounds,
  Eigen::VectorXd obj_fun_preweight,
  Eigen::VectorXd obj_fun_postweight,
  Eigen::VectorXd obj_fun_target,
  std::size_t n_approx_obj_fun_points,
  double total_budget,
  double optim_gap) {

  // constant parameters
  const std::size_t n_f = rij.rows();
  const std::size_t n_f_survey =
    std::accumulate(survey_features.begin(), survey_features.end(), 0);

  // format xgboost parameters
  std::vector<std::vector<std::string>> xgb_parameter_names;
  std::vector<std::vector<std::string>> xgb_parameter_values;
  extract_xgboost_parameters(xgb_parameters,xgb_parameter_names,
                             xgb_parameter_values);

  // format xgboost fold indices
  std::vector<std::vector<std::vector<std::size_t>>>
    xgb_train_folds2;
  std::vector<std::vector<std::vector<std::size_t>>>
    xgb_test_folds2;
  extract_k_fold_indices(xgb_train_folds, xgb_train_folds2);
  extract_k_fold_indices(xgb_test_folds, xgb_test_folds2);

  // convert environmental data to row major format
  MatrixXfRM pu_env_data2 = pu_env_data;

  // calculate value of information
  return expected_value_of_decision_given_survey_scheme(
    rij, pij, wij,
    survey_features,
    survey_sensitivity, survey_specificity,
    pu_survey_solution, pu_survey_status, pu_survey_costs,
    pu_purchase_costs, pu_purchase_locked_in, pu_env_data2,
    xgb_parameter_names, xgb_parameter_values, n_xgb_nrounds,
    xgb_train_folds2, xgb_test_folds2,
    obj_fun_preweight, obj_fun_postweight, obj_fun_target,
    n_approx_obj_fun_points, total_budget, optim_gap);
}
