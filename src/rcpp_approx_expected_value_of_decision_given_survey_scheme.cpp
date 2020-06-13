#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_sample_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_posterior_probability_matrix.h"
#include "rcpp_predict_missing_rij_data.h"
#include "rcpp_expected_value_of_action.h"

Rcpp::NumericVector approx_expected_value_of_decision_given_survey_scheme(
  Eigen::MatrixXd &rij, // observed presence/absence matrix
  Eigen::MatrixXd &pij, // prior matrix
  Eigen::MatrixXd &wij, // site weight matrix
  std::vector<bool> &survey_features, // features that we want to survey, 0/1
  Eigen::VectorXd &survey_sensitivity,
  Eigen::VectorXd &survey_specificity,
  std::vector<bool> &pu_survey_solution, // planning units to survey, 0/1
  std::vector<std::vector<std::size_t>> &pu_model_prediction_idx,
  // planning units needing prediction for rij matrix
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
  double total_budget, // total budget for surveying + monitor costs
  double optim_gap,    // optimality gap for prioritizations
  std::size_t n_approx_replicates,
  std::size_t n_approx_outcomes_per_replicate,
  std::string method_approx_outcomes) {
  // initialization
  /// constant variables
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = rij.rows();
  const std::size_t n_f_survey =
    std::accumulate(survey_features.begin(), survey_features.end(), 0);
  const std::size_t n_pu_surveyed_in_scheme =
    std::accumulate(pu_survey_solution.begin(), pu_survey_solution.end(), 0);
  const std::size_t n_vars = pu_env_data.cols();
   std::vector<std::size_t> n_pu_model_prediction(n_f);
   for (std::size_t i = 0; i < n_f; ++i)
    n_pu_model_prediction[i] = pu_model_prediction_idx[i].size();

  /// clamp number of approximation outcomes to total number of outcomes across
  /// all features
  mpz_class n_outcomes;
  mpz_class n_approx_outcomes_per_replicate2 = n_approx_outcomes_per_replicate;
  if ((n_pu_surveyed_in_scheme * n_f_survey) < 30) {
    n_states(n_pu_surveyed_in_scheme * n_f_survey, n_outcomes);
    n_outcomes = n_outcomes + 1; // increment to include final outcome
    if (cmp(n_approx_outcomes_per_replicate2, n_outcomes) > 0)
      n_approx_outcomes_per_replicate = n_outcomes.get_ui();
  }

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
  double curr_expected_value_of_decision;
  Eigen::VectorXd out(n_approx_replicates);
  Eigen::VectorXd outcome_probabiliies(n_approx_outcomes_per_replicate);
  Eigen::VectorXd outcome_values(n_approx_outcomes_per_replicate);
  std::size_t curr_n_folds;
  double curr_expected_value_of_action_given_outcome;
  double curr_probability_of_outcome;
  Eigen::MatrixXd curr_oij = rij;
  Eigen::MatrixXd curr_pij(n_f, n_pu);
  curr_pij.setConstant(-100.0);
  Eigen::MatrixXd curr_total_probability_of_model_positive(n_f_survey, n_pu);
  Eigen::MatrixXd curr_total_probability_of_model_negative(n_f_survey, n_pu);
  std::vector<bool> curr_solution(n_pu);
  std::vector<mpz_class> feature_outcome_idx(n_f_survey);
  std::vector<mpz_class> outcomes(n_approx_outcomes_per_replicate);
  model_yhat_map model_yhat;
  model_yhat.reserve(n_f_survey * 1000);
  model_performance_map model_performance;
  model_performance.reserve(n_f_survey * 1000);
  std::vector<std::vector<BoosterHandle> *> curr_models(n_f_survey);
  Eigen::VectorXd curr_model_sensitivity(n_f_survey);
  Eigen::VectorXd curr_model_specificity(n_f_survey);

  // preliminary processing
  /// create subset of prior matrix for just the species that need surveys
  Eigen::MatrixXd pij_survey_species_subset(n_f_survey, n_pu);
  for (std::size_t i = 0; i < n_f_survey; ++i)
    pij_survey_species_subset.row(i) = pij.row(survey_features_idx[i]);

  // subset environmental data for planning unit predictions
  std::vector<MatrixXfRM> pu_predict_env_data(n_f_survey);
  std::size_t curr_n;
  for (std::size_t i = 0; i < n_f_survey; ++i) {
    /// prepare matrix
    curr_n = n_pu_model_prediction[survey_features_idx[i]];
    pu_predict_env_data[i].resize(curr_n, n_vars);
    /// store environmental values for species needing predictions
    for (std::size_t j = 0; j < curr_n; ++j)
      pu_predict_env_data[i].row(j).array() =
        pu_env_data.row(pu_model_prediction_idx[i][j]).array();
  }

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
    remaining_budget, optim_gap);

  /// overwrite outcome data with prior data for features we are
  /// not interested in surveying: planning units needing model predictions
  for (std::size_t i = 0; i < n_f; ++i)
    if (!survey_features[i])
      for (std::size_t j = 0; j < n_pu_surveyed_in_scheme; ++j)
        curr_oij(i, pu_survey_solution_idx[j]) =
          pij(i, pu_survey_solution_idx[j]);

  /// overwrite outcome data with prior data for features we are
  /// not interested in surveying: planning units in the survey scheme
  for (std::size_t i = 0; i < n_f; ++i)
    if (!survey_features[i])
      for (std::size_t j = 0; j < n_pu_model_prediction[i]; ++j)
        curr_oij(i, pu_model_prediction_idx[i][j]) =
          pij(i, pu_model_prediction_idx[i][j]);

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
  /// store prior probabilities for surveyed planning units and survey features,
  /// to generate the outcomes
  Eigen::MatrixXd survey_pij(n_f_survey * n_pu_surveyed_in_scheme, 1);
  for (std::size_t i = 0; i < (n_f_survey * n_pu_surveyed_in_scheme); ++i)
    survey_pij(i) = pij(rij_outcome_idx[i]);

  // main processing
  for (std::size_t r = 0; r < n_approx_replicates; ++r) {
    /// initialize current value
    curr_expected_value_of_decision = std::numeric_limits<double>::infinity();

    /// generate outcomes
    sample_n_states(
      n_approx_outcomes_per_replicate, survey_pij, method_approx_outcomes,
      outcomes);

    // calculate values based outcomes
    for (std::size_t o = 0; o < n_approx_outcomes_per_replicate; ++o) {
      /// generate the o'th outcome from surveying the planning units across
      /// all species
      nth_state_sparse(outcomes[o], rij_outcome_idx, curr_oij);

      // find out the outcome for each feature seperately
      which_feature_state(
        curr_oij, survey_features_idx, pu_survey_solution_idx,
        feature_outcome_idx);

      // fit models for the feature's outcomes if needed
      fit_xgboost_models_and_assess_performance(
        curr_oij, wij,
        pu_env_data, pu_predict_env_data,
        survey_features_idx, feature_outcome_idx,
        xgb_parameter_names, xgb_parameter_values, n_xgb_nrounds,
        xgb_train_folds, xgb_test_folds,
        model_yhat, model_performance,
        curr_model_sensitivity, curr_model_specificity);

      /// update model sensitivity and speciifcity values based on
      /// survey sensitivity and speciifcity values, because the
      /// model values are dependent on the survey data
      for (std::size_t i = 0; i < n_f_survey; ++i)
        curr_model_sensitivity[i] *= survey_sensitivity[survey_features_idx[i]];
      for (std::size_t i = 0; i < n_f_survey; ++i)
        curr_model_specificity[i] *= survey_specificity[survey_features_idx[i]];

      /// generate modelled predictions for survey species
      predict_missing_rij_data(
        curr_oij, survey_features_idx, feature_outcome_idx,
        pu_model_prediction_idx, model_yhat);
      assert_valid_probability_data(curr_oij, "issue predicting missing data");

      /// calculate total probability of models' positive results
      total_probability_of_positive_model_result(
        pij_survey_species_subset,
        curr_model_sensitivity, curr_model_specificity,
        curr_total_probability_of_model_positive);
      assert_valid_probability_data(curr_total_probability_of_model_positive,
                                    "issue calculating total model positives");

      /// calculate total probability of models' negative results
      total_probability_of_negative_model_result(
        pij_survey_species_subset,
        curr_model_sensitivity, curr_model_specificity,
        curr_total_probability_of_model_negative);
      assert_valid_probability_data(curr_total_probability_of_model_negative,
                                    "issue calculating total model negatives");

      /// generate posterior data
      posterior_probability_matrix(
        rij, pij, curr_oij,
        pu_survey_solution,
        survey_features, survey_features_rev_idx,
        survey_sensitivity, survey_specificity,
        total_probability_of_survey_positive,
        total_probability_of_survey_negative,
        curr_model_sensitivity, curr_model_specificity,
        curr_total_probability_of_model_positive,
        curr_total_probability_of_model_negative,
        curr_pij);
      assert_valid_probability_data(
        curr_pij, "issue calculating posterior probabilities");

      /// generate prioritisation
      prioritize.add_rij_data(curr_pij);
      prioritize.solve();
      prioritize.get_solution(curr_solution);

      /// calculate expected value of the prioritisation
      curr_expected_value_of_action_given_outcome =
        expected_value_of_action(
          curr_solution, curr_pij, obj_fun_preweight, obj_fun_postweight,
          obj_fun_target);

      /// calculate likelihood of outcome
      curr_probability_of_outcome = log_probability_of_outcome(
        curr_oij, total_probability_of_survey_positive_log,
        total_probability_of_survey_negative_log, rij_outcome_idx);

      /// calculate values of action
      outcome_values[o] = curr_expected_value_of_action_given_outcome;
      outcome_probabiliies[o] = curr_probability_of_outcome;
    }

    // apply correction to approximate the true value
    outcome_probabiliies.array() -= log_sum(outcome_probabiliies);
    outcome_probabiliies.array() += outcome_values.array().log();

    // store result
    out[r] = log_sum(outcome_probabiliies);
  }

  // exports
  out.array() = out.array().exp();
  return Rcpp::wrap(out);
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_approx_expected_value_of_decision_given_survey_scheme(
  Eigen::MatrixXd rij,
  Eigen::MatrixXd pij,
  Eigen::MatrixXd wij,
  std::vector<bool> survey_features,
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  std::vector<bool> pu_survey_solution,
  Rcpp::List pu_model_prediction,
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
  double total_budget,
  double optim_gap,
  std::size_t n_approx_replicates,
  std::size_t n_approx_outcomes_per_replicate,
  std::string method_approx_outcomes) {

  // constant parameters
  const std::size_t n_f = rij.rows();
  const std::size_t n_pu = rij.cols();
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

  // format pu model prediction indices
  std::vector<std::vector<std::size_t>> pu_model_prediction_idx;
  extract_list_of_list_of_indices(pu_model_prediction, pu_model_prediction_idx);

  // convert environmental data to row major format
  MatrixXfRM pu_env_data2 = pu_env_data;

  // increment model weights for planning units selected in survey scheme
  // and species that will be considerd in future surveys
  for (std::size_t i = 0; i < n_f; ++i)
    for (std::size_t j = 0; j < n_pu; ++j)
      wij(i, j) +=
        static_cast<double>(survey_features[i] && pu_survey_solution[j]);

  // calculate value of information
  return approx_expected_value_of_decision_given_survey_scheme(
    rij, pij, wij,
    survey_features,
    survey_sensitivity, survey_specificity,
    pu_survey_solution, pu_model_prediction_idx,
    pu_survey_costs, pu_purchase_costs, pu_purchase_locked_in, pu_env_data2,
    xgb_parameter_names, xgb_parameter_values, n_xgb_nrounds,
    xgb_train_folds2, xgb_test_folds2,
    obj_fun_preweight, obj_fun_postweight, obj_fun_target,
    total_budget, optim_gap,
    n_approx_replicates, n_approx_outcomes_per_replicate,
    method_approx_outcomes);
}
