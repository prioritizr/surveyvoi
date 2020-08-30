#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_sample_states.h"
#include "rcpp_probability.h"
#include "rcpp_heuristic_prioritization.h"
#include "rcpp_posterior_probability_matrix.h"
#include "rcpp_predict_missing_rij_data.h"
#include "rcpp_expected_value_of_action.h"
#include "expected_value_given_survey_outcome.h"

Rcpp::NumericVector approx_expected_value_of_decision_given_survey_scheme(
  Eigen::MatrixXd &dij, // proportion of detections matrix
  Eigen::MatrixXd &nij, // number of surveys matrix
  Eigen::MatrixXd &pij, // prior matrix
  std::vector<bool> &survey_features, // features that we want to survey, 0/1
  Eigen::VectorXd &survey_sensitivity,
  Eigen::VectorXd &survey_specificity,
  std::vector<bool> &pu_survey_solution, // planning units to survey, 0/1
  std::vector<std::vector<std::size_t>> &pu_model_prediction_idx,
  // planning units needing prediction for rij matrix
  Eigen::VectorXd &pu_survey_costs,   // cost of surveying planning units
  Eigen::VectorXd &pu_purchase_costs, // cost of purchasing planning units
  Eigen::VectorXd &pu_purchase_locked_in,  // planning units that are locked in
  Eigen::VectorXd &pu_purchase_locked_out,  // planning units that locked out
  MatrixXfRM &pu_env_data, // environmental data
  std::vector<std::string> &xgb_parameter_names, // xgboost parameter names
  MatrixXs &xgb_parameter_values, // xgboost parameter combination values
  std::vector<std::size_t> &n_xgb_rounds,  // xgboost training rounds
  std::vector<std::size_t> &n_xgb_early_stopping_rounds,  // xgboost early stop
  std::vector<std::vector<std::vector<std::size_t>>> &xgb_train_folds,
  std::vector<std::vector<std::vector<std::size_t>>> &xgb_test_folds,
  Eigen::VectorXi &obj_fun_target,  // objective function calculation term
  double total_budget, // total budget for surveying + monitor costs
  std::size_t n_approx_replicates,
  std::size_t n_approx_outcomes_per_replicate,
  std::string method_approx_outcomes,
  double seed) {
  // initialization
  /// constant variables
  const std::size_t n_pu = dij.cols();
  const std::size_t n_f = dij.rows();
  const std::size_t n_vars = pu_env_data.cols();
  const std::size_t n_f_survey =
    std::accumulate(survey_features.begin(), survey_features.end(), 0);
  const std::size_t n_pu_surveyed_in_scheme =
    std::accumulate(pu_survey_solution.begin(), pu_survey_solution.end(), 0);
   std::vector<std::size_t> n_pu_model_prediction(n_f);
   for (std::size_t i = 0; i < n_f; ++i)
    n_pu_model_prediction[i] = pu_model_prediction_idx[i].size();

  /// clamp number of approximation outcomes to total number of outcomes across
  /// all features
  std::size_t n_approx_survey_outcomes_per_replicate =
    n_approx_outcomes_per_replicate;
  mpz_class n_outcomes;
  mpz_class n_approx_outcomes_per_replicate2 =
    n_approx_survey_outcomes_per_replicate;
  if ((n_pu_surveyed_in_scheme * n_f_survey) < 20) {
    n_states(n_pu_surveyed_in_scheme * n_f_survey, n_outcomes);
    n_outcomes = n_outcomes + 1; // increment to include final outcome
    if (cmp(n_approx_outcomes_per_replicate2, n_outcomes) > 0)
      n_approx_survey_outcomes_per_replicate = n_outcomes.get_ui();
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
  std::size_t curr_n_approx_model_outcomes_per_replicate;
  Eigen::VectorXd out(n_approx_replicates);
  Eigen::VectorXd outcome_probabilities;
  Eigen::VectorXd outcome_values;
  std::size_t curr_n_folds;
  double curr_expected_value_of_action_given_outcome;
  double curr_probability_of_outcome;
  Eigen::MatrixXd curr_dij = dij;
  Eigen::MatrixXd curr_nij = nij;
  Eigen::MatrixXd curr_oij = pij;
  Eigen::MatrixXd curr_pij(n_f, n_pu);
  curr_pij.setConstant(-100.0);
  Eigen::MatrixXd curr_total_probability_of_model_positive(n_f_survey, n_pu);
  Eigen::MatrixXd curr_total_probability_of_model_negative(n_f_survey, n_pu);
  Eigen::MatrixXd curr_total_probability_of_model_positive_log;
  Eigen::MatrixXd curr_total_probability_of_model_negative_log;
  Eigen::MatrixXd curr_modelled_probabilities;
  std::vector<bool> curr_solution(n_pu);
  std::vector<mpz_class> feature_outcome_idx(n_f_survey);
  std::vector<mpz_class> outcomes;
  std::vector<mpz_class> model_outcomes;
  model_yhat_map model_yhat;
  model_yhat.reserve(n_f_survey * 1000);
  model_performance_map model_performance;
  model_performance.reserve(
    n_f_survey * n_approx_replicates * n_approx_survey_outcomes_per_replicate);
  std::vector<std::vector<BoosterHandle> *> curr_models(n_f_survey);
  Eigen::VectorXd curr_model_sensitivity(n_f_survey);
  Eigen::VectorXd curr_model_specificity(n_f_survey);
  std::vector<std::size_t> curr_model_rij_idx;
  std::size_t curr_n_survey_outcomes;

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

  // prepare xgboost data structures for model training
  std::vector<std::vector<Eigen::VectorXf>> train_y;
  std::vector<std::vector<Eigen::VectorXf>> train_w;
  std::vector<std::vector<MatrixXfRM>> train_x;
  extract_k_fold_y_data_from_indices(
    xgb_train_folds, survey_features_idx, train_y);
  extract_k_fold_x_data_from_indices(
    pu_env_data, xgb_train_folds, survey_features_idx, train_x);

  // prepare xgboost data structures for model evaluation
  std::vector<std::vector<Eigen::VectorXf>> test_y;
  std::vector<std::vector<Eigen::VectorXf>> test_w;
  std::vector<std::vector<MatrixXfRM>> test_x;
  extract_k_fold_y_data_from_indices(
    xgb_test_folds, survey_features_idx, test_y);
  extract_k_fold_x_data_from_indices(
    pu_env_data, xgb_test_folds, survey_features_idx, test_x);

  // update survey data with new numbers of surveys
  for (auto itr = rij_outcome_idx.begin(); itr != rij_outcome_idx.end();
       ++itr) {
    ++curr_nij(*itr);
  }

  // main processing
  for (std::size_t r = 0; r < n_approx_replicates; ++r) {
    /// generate outcomes,
    // manually adjust seed to ensure different samples per replicate
    sample_n_states(
      n_approx_survey_outcomes_per_replicate, survey_pij,
      method_approx_outcomes, seed + static_cast<double>(r), outcomes);
    set_seed(seed);
    curr_n_survey_outcomes = outcomes.size();
    outcome_probabilities.resize(curr_n_survey_outcomes);
    outcome_values.resize(curr_n_survey_outcomes);

    // calculate values based outcomes
    for (std::size_t o = 0; o < curr_n_survey_outcomes; ++o) {
      /// generate the o'th outcome from surveying the planning units across
      /// all species
      nth_state_sparse(outcomes[o], rij_outcome_idx, curr_oij);

      /// calculate likelihood of outcome
      curr_probability_of_outcome = log_probability_of_outcome(
        curr_oij, total_probability_of_survey_positive_log,
        total_probability_of_survey_negative_log, rij_outcome_idx);

      // find out the outcome for each feature seperately
      which_feature_state(
        curr_oij, survey_features_idx, pu_survey_solution_idx,
        feature_outcome_idx);

      /// update survey data with new simulated outcomes
      for (auto itr = rij_outcome_idx.begin(); itr != rij_outcome_idx.end();
           ++itr) {
        curr_dij(*itr) =
          ((dij(*itr) * nij(*itr)) + curr_oij(*itr)) / curr_nij(*itr);
      }

      /// initialize posterior probability matrix with new survey results
      /// and the prior data for cells that don't need
      initialize_posterior_probability_matrix(
        nij, pij, curr_oij,
        pu_survey_solution,
        survey_features, survey_features_rev_idx,
        survey_sensitivity, survey_specificity,
        total_probability_of_survey_positive,
        total_probability_of_survey_negative,
        curr_pij);

      // prepare weight data for modelling
      extract_k_fold_train_w_data_from_indices(
        curr_pij, xgb_train_folds, survey_features_idx, train_w);
      extract_k_fold_test_w_data_from_indices(
        curr_dij, curr_nij, xgb_test_folds, survey_features_idx, test_w);

      /// fit models for the feature's outcomes if needed
      fit_xgboost_models_and_assess_performance(
        survey_features_idx,
        survey_sensitivity, survey_specificity,
        feature_outcome_idx,
        xgb_parameter_names, xgb_parameter_values,
        n_xgb_rounds, n_xgb_early_stopping_rounds,
        train_x, train_y, train_w, test_x, test_y, test_w, pu_predict_env_data,
        model_yhat, model_performance,
        curr_model_sensitivity, curr_model_specificity);

      /// generate model predictions for species we are surveying
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

      /// find idices in the rij matrix that use model estimates
      find_rij_idx_based_on_models(
        nij,
        pu_survey_solution,
        survey_features, survey_features_rev_idx,
        survey_sensitivity, survey_specificity,
        curr_model_sensitivity, curr_model_specificity,
        curr_model_rij_idx);

      /// create posterior matrix with most likely model outcomes
      update_model_posterior_probabilities(
        curr_model_rij_idx,
        pij, curr_oij,
        survey_features_rev_idx,
        curr_model_sensitivity, curr_model_specificity,
        curr_total_probability_of_model_positive,
        curr_total_probability_of_model_negative,
        curr_pij);
        assert_valid_probability_data(
        curr_pij, "issue calculating posterior probabilities");

      // generate prioritisation
      greedy_heuristic_prioritization(
        curr_pij, pu_purchase_costs, pu_purchase_locked_in,
        pu_purchase_locked_out, obj_fun_target, remaining_budget,
        curr_solution);

      /// caclualate expected value of prioritisation given survey outcome
      if (curr_model_rij_idx.size() == 0) {
        /// if no modelled probabilities are being used in the prioritisation,
        /// then calculate the expected value using posterior with the
        /// the given survey outcome
        curr_expected_value_of_action_given_outcome =
          std::log(expected_value_of_action(curr_solution, curr_pij,
            obj_fun_target));
      } else if (curr_model_rij_idx.size() <= 5) {
        curr_expected_value_of_action_given_outcome =
          exact_expected_value_given_survey_outcome_with_model_estimates(
            pij,
            curr_pij,
            curr_solution,
            survey_features_rev_idx,
            curr_model_rij_idx,
            curr_model_sensitivity,
            curr_model_specificity,
            curr_total_probability_of_model_positive,
            curr_total_probability_of_model_negative,
            obj_fun_target);
      } else {
        curr_expected_value_of_action_given_outcome =
          approx_expected_value_given_survey_outcome_with_model_estimates(
            pij,
            curr_pij,
            curr_solution,
            survey_features_rev_idx,
            curr_model_rij_idx,
            curr_model_sensitivity,
            curr_model_specificity,
            curr_total_probability_of_model_positive,
            curr_total_probability_of_model_negative,
            obj_fun_target,
            n_approx_outcomes_per_replicate,
            method_approx_outcomes,
            seed + r + o);
        set_seed(seed);
      }
      /// calculate expected value of decision for the o'th survey outcome
      /// given uncertainty in model predictions
      outcome_values[o] = curr_expected_value_of_action_given_outcome;
      outcome_probabilities[o] = curr_probability_of_outcome;
    }

    // standardize survye outcome data
    outcome_probabilities.array() -= log_sum(outcome_probabilities);
    outcome_probabilities.array() += outcome_values.array();

    // store result
    out[r] = log_sum(outcome_probabilities);
  }

  // exports
  out.array() = out.array().exp();
  return Rcpp::wrap(out);
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_approx_expected_value_of_decision_given_survey_scheme(
  Eigen::MatrixXd dij,
  Eigen::MatrixXd nij,
  Eigen::MatrixXd pij,
  std::vector<bool> survey_features,
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  std::vector<bool> pu_survey_solution,
  Rcpp::List pu_model_prediction,
  Eigen::VectorXd pu_survey_costs,
  Eigen::VectorXd pu_purchase_costs,
  Eigen::VectorXd pu_purchase_locked_in,
  Eigen::VectorXd pu_purchase_locked_out,
  Eigen::MatrixXf pu_env_data_raw,
  std::vector<std::string> xgb_parameter_names,
  Rcpp::CharacterMatrix xgb_parameter_values,
  std::vector<std::size_t> n_xgb_rounds,
  std::vector<std::size_t> n_xgb_early_stopping_rounds,
  Rcpp::List xgb_train_folds,
  Rcpp::List xgb_test_folds,
  Eigen::VectorXi obj_fun_target,
  double total_budget,
  std::size_t n_approx_replicates,
  std::size_t n_approx_outcomes_per_replicate,
  std::string method_approx_outcomes,
  double seed) {

  // constant parameters
  MatrixXfRM pu_env_data = pu_env_data_raw;
  const std::size_t n_f = survey_features.size();
  std::vector<std::size_t> survey_features_idx;
  survey_features_idx.reserve(n_f);
  for (std::size_t i = 0; i < n_f; ++i)
    if (survey_features[i])
      survey_features_idx.push_back(i);
  survey_features_idx.shrink_to_fit();
  const std::size_t n_f_survey = survey_features_idx.size();
  const std::size_t n_f_outcomes = 1;
  std::vector<mpz_class> feature_outcome_idx(n_f_survey, 0);

  // extract xgboost parameter values
  MatrixXs xgb_parameter_values2(
    xgb_parameter_values.rows(), xgb_parameter_values.cols());
  for (std::size_t i = 0; i != xgb_parameter_values2.size(); ++i)
    xgb_parameter_values2(i) = Rcpp::as<std::string>(xgb_parameter_values[i]);

  // format xgboost fold indices
  std::vector<std::vector<std::vector<std::size_t>>> xgb_train_folds2;
  std::vector<std::vector<std::vector<std::size_t>>> xgb_test_folds2;
  extract_k_fold_indices(xgb_train_folds, xgb_train_folds2);
  extract_k_fold_indices(xgb_test_folds, xgb_test_folds2);

  // format pu model prediction indices
  std::vector<std::vector<std::size_t>> pu_model_prediction_idx;
  extract_list_of_list_of_indices(pu_model_prediction, pu_model_prediction_idx);

  // calculate value of information
  return approx_expected_value_of_decision_given_survey_scheme(
    dij, nij, pij,
    survey_features,
    survey_sensitivity, survey_specificity,
    pu_survey_solution, pu_model_prediction_idx,
    pu_survey_costs, pu_purchase_costs,
    pu_purchase_locked_in, pu_purchase_locked_out,
    pu_env_data,
    xgb_parameter_names, xgb_parameter_values2,
    n_xgb_rounds, n_xgb_early_stopping_rounds,
    xgb_train_folds2, xgb_test_folds2,
    obj_fun_target, total_budget,
    n_approx_replicates, n_approx_outcomes_per_replicate,
    method_approx_outcomes, seed);
}
