#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_sample_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_posterior_probability_matrix.h"
#include "rcpp_predict_missing_rij_data.h"
#include "rcpp_expected_value_of_action.h"

// [[Rcpp::export]]
Rcpp::NumericVector
  rcpp_approx_expected_value_of_decision_given_survey_scheme_n_states(
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
  double optim_gap,
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
  const std::size_t n_pu_surveyed_already =
    static_cast<std::size_t>((pu_survey_status.array()).sum());

  /// format xgboost parameters
  std::vector<std::vector<std::string>> xgb_parameter_names;
  std::vector<std::vector<std::string>> xgb_parameter_values;
  extract_xgboost_parameters(xgb_parameters, xgb_parameter_names,
                             xgb_parameter_values);

  /// format xgboost fold indices
  std::vector<std::vector<std::vector<std::size_t>>>
    xgb_train_folds2;
  std::vector<std::vector<std::vector<std::size_t>>>
    xgb_test_folds2;
  extract_k_fold_indices(xgb_train_folds, xgb_train_folds2);
  extract_k_fold_indices(xgb_test_folds, xgb_test_folds2);

  /// convert environmental data to row major format
  MatrixXfRM pu_env_data2 = pu_env_data;

  /// integer over-flow checks, highest std::size_t value is 1e+18
  if (n_pu_surveyed_in_scheme > 20)
    Rcpp::stop("number of planning units in selected in survey scheme is too large (i.e. >20)");

  /// calculate number of outcomes for a given feature
  const std::size_t n_f_outcomes = n_states(n_pu_surveyed_in_scheme) + 1;

  /// clamp number of approximation outcomes to total number of outcomes across
  /// all features
  mpz_class n_outcomes;
  n_states(n_pu_surveyed_in_scheme * n_f_survey, n_outcomes);
  n_outcomes = n_outcomes + 1; // increment to include final outcome
  n_approx_outcomes_per_replicate = std::min(
    n_approx_outcomes_per_replicate, n_outcomes.get_ui());

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
  Eigen::VectorXd out(n_approx_replicates);
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
  model_beta_map model_beta;
  model_performance_map model_performance;
  std::vector<std::vector<BoosterHandle> *> curr_models(n_f_survey);
  Eigen::VectorXd curr_model_sensitivity(n_f_survey);
  Eigen::VectorXd curr_model_specificity(n_f_survey);

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
  /// store prior probabilities for surveyed planning units and survey features,
  /// to generate the outcomes
  Eigen::MatrixXd survey_pij(n_f_survey * n_pu_surveyed_in_scheme, 1);
  for (std::size_t i = 0; i < (n_f_survey * n_pu_surveyed_in_scheme); ++i)
    survey_pij(i) = pij(rij_outcome_idx[i]);

  // main processing
  for (std::size_t r = 0; r < n_approx_replicates; ++r) {
    /// initialize current value
    out[r] = std::numeric_limits<double>::infinity();
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
        curr_oij, wij, pu_env_data2,
        survey_features_idx, feature_outcome_idx,
        xgb_parameter_names, xgb_parameter_values, n_xgb_nrounds,
        xgb_train_folds2, xgb_test_folds2,
        model_beta, model_performance,
        curr_models, curr_model_sensitivity, curr_model_specificity);

      /// generate modelled predictions for survey species
      predict_missing_rij_data(
        curr_oij, pu_env_data2, survey_features_idx, pu_model_prediction_idx,
        curr_models);
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
        std::log(expected_value_of_action(
          curr_solution, curr_pij, obj_fun_preweight, obj_fun_postweight,
          obj_fun_target));

      /// calculate likelihood of outcome
      curr_probability_of_outcome = log_probability_of_outcome(
        curr_oij, total_probability_of_survey_positive_log,
        total_probability_of_survey_negative_log, rij_outcome_idx);

      /// calculate expected value of action
      if (std::isinf(out[r])) {
        out[r] = curr_expected_value_of_action_given_outcome +
                 curr_probability_of_outcome;
      } else {
        out[r] = log_sum(out[r],
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
    }
  }

  // clean-up
  for (auto itr = model_beta.begin(); itr != model_beta.end(); ++itr) {
    curr_n_folds = itr->second->size();
    for (std::size_t k = 0; k < curr_n_folds; ++k)
      XGBoosterFree(itr->second->at(k));
    delete itr->second;
  }

  // exports
  out.array() = out.array().exp();
  return Rcpp::wrap(out);
}
