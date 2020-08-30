#include "expected_value_given_survey_outcome.h"

double exact_expected_value_given_survey_outcome_with_model_estimates(
  Eigen::MatrixXd &pij,
  Eigen::MatrixXd &curr_pij,
  std::vector<bool> &curr_solution,
  std::vector<std::size_t> &survey_features_rev_idx,
  std::vector<std::size_t> &curr_model_rij_idx,
  Eigen::VectorXd &curr_model_sensitivity,
  Eigen::VectorXd &curr_model_specificity,
  Eigen::MatrixXd &curr_total_probability_of_model_positive,
  Eigen::MatrixXd &curr_total_probability_of_model_negative,
  Eigen::VectorXi &obj_fun_target) {
  // Initialization
  // declare varaibles
  mpz_class n_model_outcomes;
  const std::size_t n_model_rij_idx = curr_model_rij_idx.size();
  Eigen::MatrixXd curr_total_probability_of_model_positive_log;
  Eigen::MatrixXd curr_total_probability_of_model_negative_log;
  double curr_expected_value = std::numeric_limits<double>::infinity();
  double curr_prob, curr_value;

  // Preliminary processing
  /// calculate total number of possible model outcomes
  n_states(n_model_rij_idx, n_model_outcomes);
  n_model_outcomes = n_model_outcomes + 1; // increment to include final

  /// calculate log of model total probabilities
  curr_total_probability_of_model_positive_log =
    curr_total_probability_of_model_positive;
  log_matrix(curr_total_probability_of_model_positive_log);
  curr_total_probability_of_model_negative_log =
    curr_total_probability_of_model_negative;
  log_matrix(curr_total_probability_of_model_negative_log);

  // Main processing
  /// iterate over each model outcome
  mpz_class oo = 0;
  while (cmp(oo, n_model_outcomes) < 0) {
    /// generate oo'th model outcome
    nth_state_sparse(oo, curr_model_rij_idx, curr_pij);

    /// calculate likelihood of model outcome
    curr_prob =
      log_probability_of_model_outcome(
        curr_pij, survey_features_rev_idx,
        curr_total_probability_of_model_positive_log,
        curr_total_probability_of_model_negative_log,
        curr_model_rij_idx);

    /// generate posterior probability given model outcome
    update_model_posterior_probabilities(
      curr_model_rij_idx,
      pij, curr_pij,
      survey_features_rev_idx,
      curr_model_sensitivity, curr_model_specificity,
      curr_total_probability_of_model_positive,
      curr_total_probability_of_model_negative,
      curr_pij);
      assert_valid_probability_data(
        curr_pij, "issue calculating model posterior probabilities");

    /// calculate expected value of action given model outcome
    curr_value = std::log(expected_value_of_action(
        curr_solution, curr_pij, obj_fun_target));

    /// add conditional value to the runing total across all model outcomes
    if (std::isinf(curr_expected_value)) {
      curr_expected_value = curr_value + curr_prob;
    } else {
      curr_expected_value =
        log_sum(curr_expected_value, curr_value + curr_prob);
    }

    // increment model outcome counter
    oo = oo + 1;
  }

  // Exports
  /// return result
  return curr_expected_value;
}

double approx_expected_value_given_survey_outcome_with_model_estimates(
  Eigen::MatrixXd &pij,
  Eigen::MatrixXd &curr_pij,
  std::vector<bool> &curr_solution,
  std::vector<std::size_t> &survey_features_rev_idx,
  std::vector<std::size_t> &curr_model_rij_idx,
  Eigen::VectorXd &curr_model_sensitivity,
  Eigen::VectorXd &curr_model_specificity,
  Eigen::MatrixXd &curr_total_probability_of_model_positive,
  Eigen::MatrixXd &curr_total_probability_of_model_negative,
  Eigen::VectorXi &obj_fun_target,
  std::size_t n_approx_outcomes_per_replicate,
  std::string &method_approx_outcomes,
  std::size_t seed) {
  // Initialization
  // declare variables
  const std::size_t n_model_rij_idx = curr_model_rij_idx.size();
  mpz_class n1, n2;
  std::size_t n_model_outcomes =
    n_approx_outcomes_per_replicate;
  Eigen::MatrixXd curr_modelled_probabilities(n_model_rij_idx, 1);

  /// prepare probabilities for generating the model outcomes
  for (std::size_t ii = 0; ii < n_model_rij_idx; ++ii)
    curr_modelled_probabilities(ii) = curr_pij(curr_model_rij_idx[ii]);

  /// calculate log of model total probabilities
  Eigen::MatrixXd curr_total_probability_of_model_positive_log =
    curr_total_probability_of_model_positive;
  log_matrix(curr_total_probability_of_model_positive_log);
  Eigen::MatrixXd curr_total_probability_of_model_negative_log =
    curr_total_probability_of_model_negative;
  log_matrix(curr_total_probability_of_model_negative_log);

  // Preliminary processing
  /// determine number of model outcomes
  n1 = n_approx_outcomes_per_replicate;
  if (n_model_rij_idx < 20) {
    n_states(n_model_rij_idx, n2);
    n2 = n2 + 1; // increment to include final outcome
    if (cmp(n2, n1) < 0)
      n_model_outcomes = n2.get_ui();
  }

  /// generate a set of model outcomes
  /// manually adjust seed to ensure different samples per replicate
  std::vector<mpz_class> model_outcomes;
  sample_n_states(
    n_model_outcomes, curr_modelled_probabilities, method_approx_outcomes,
    seed, model_outcomes);
  n_model_outcomes = model_outcomes.size();

  /// create data structures for calculating results
  Eigen::VectorXd curr_prob(n_model_outcomes);
  Eigen::VectorXd curr_value(n_model_outcomes);

  /// iterate over each model outcome and calculate the
  /// probability of the outcome occurring and the expected value
  /// of the prioritisation given the model outcome
  for (std::size_t oo = 0; oo < n_model_outcomes; ++oo) {

    /// generate oo'th model outcome
    nth_state_sparse(model_outcomes[oo], curr_model_rij_idx, curr_pij);

    /// calculate likelihood of model outcome
    curr_prob[oo] =
      log_probability_of_model_outcome(
        curr_pij, survey_features_rev_idx,
        curr_total_probability_of_model_positive_log,
        curr_total_probability_of_model_negative_log,
        curr_model_rij_idx);

    /// generate posterior probability given model outcome
    update_model_posterior_probabilities(
      curr_model_rij_idx,
      pij, curr_pij,
      survey_features_rev_idx,
      curr_model_sensitivity, curr_model_specificity,
      curr_total_probability_of_model_positive,
      curr_total_probability_of_model_negative,
      curr_pij);
      assert_valid_probability_data(
        curr_pij, "issue calculating model posterior probabilities");

    /// calculate expected value of action given model outcome
    curr_value[oo] =
      expected_value_of_action(
        curr_solution, curr_pij, obj_fun_target);
  }

  // Exports
  /// calculate approximate value for the o'th survey outcome
  /// given uncertainty in the model outcomes
  curr_prob.array() -= log_sum(curr_prob);
  curr_prob.array() += curr_value.array().log();
  return log_sum(curr_prob);
}
