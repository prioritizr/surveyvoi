#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_heuristic_prioritization.h"
#include "rcpp_update_posterior_probability_matrix.h"
#include "rcpp_expected_value_of_action.h"

// [[Rcpp::export]]
double rcpp_expected_value_of_decision_given_survey_scheme(
  Eigen::MatrixXd pij, // prior matrix
  std::vector<bool> survey_features, // features that we want to survey, 0/1
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  std::vector<bool> pu_survey_solution, // planning units to survey, 0/1
  Eigen::VectorXd pu_survey_costs,   // cost of surveying planning units
  Eigen::VectorXd pu_purchase_costs, // cost of purchasing planning units
  Eigen::VectorXd pu_purchase_locked_in,  // planning units that are locked in
  Eigen::VectorXd pu_purchase_locked_out,  // planning units that locked out
  Eigen::VectorXi obj_fun_target,  // objective function calculation term
  double total_budget // total budget for surveying + monitor costs
) {
  // initialization
  /// constant variables
  const std::size_t n_pu = pij.cols();
  const std::size_t n_f = pij.rows();
  const std::size_t n_f_survey =
    std::accumulate(survey_features.begin(), survey_features.end(), 0);
  const std::size_t n_pu_surveyed_in_scheme =
    std::accumulate(pu_survey_solution.begin(), pu_survey_solution.end(), 0);

  /// integer over-flow checks, highest std::size_t value is 1e+18
  if (n_pu_surveyed_in_scheme > 20)
    Rcpp::stop("number of planning units in selected in survey scheme is too large (i.e. >20)");

  /// calculate number of outcomes for a given feature
  const std::size_t n_f_outcomes = n_states(n_pu_surveyed_in_scheme) + 1;

  /// calculate total number of outcomes across all features
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

  //// store indices of features that need surveying
  std::vector<std::size_t> survey_features_idx;
  survey_features_idx.reserve(n_f);
  for (std::size_t i = 0; i < n_f; ++i)
    if (survey_features[i])
      survey_features_idx.push_back(i);
  survey_features_idx.shrink_to_fit();

  /// indices for features in sparse format for reverse lookup
  std::vector<std::size_t> survey_features_rev_idx(n_f, 1e5);
  create_reverse_lookup_id(survey_features, survey_features_rev_idx);

  // preliminary processing
  /// create subset of prior matrix for just the species that need surveys
  Eigen::MatrixXd pij_survey_species_subset(n_f_survey, n_pu);
  for (std::size_t i = 0; i < n_f_survey; ++i)
    pij_survey_species_subset.row(i) = pij.row(survey_features_idx[i]);

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

  /// declare temporary variables used in the main loop
  double curr_expected_value_of_action_given_outcome;
  double curr_probability_of_outcome;
  double curr_expected_value_of_decision =
    std::numeric_limits<double>::infinity();
  Eigen::MatrixXd curr_oij = pij;
  Eigen::MatrixXd curr_wij = pij;
  std::vector<bool> curr_solution(n_pu);

  // main processing
  mpz_class o = 0;
  while (cmp(o, n_outcomes) < 0) {
    /// generate the o'th outcome from surveying the planning units
    nth_state_sparse(o, rij_outcome_idx, curr_oij);

    /// calculate likelihood of outcome
    curr_probability_of_outcome = log_probability_of_outcome(
      curr_oij, total_probability_of_survey_positive_log,
      total_probability_of_survey_negative_log, rij_outcome_idx);

    /// update posterior matrix based on survey outcome
    update_posterior_probability_matrix(
      pij, curr_oij,
      rij_outcome_idx,
      survey_sensitivity,
      survey_specificity,
      total_probability_of_survey_positive,
      total_probability_of_survey_negative,
      curr_wij);
      assert_valid_probability_data(
      curr_wij, "issue calculating posterior probabilities");

    /// generate prioritisation given survey outcome
    greedy_heuristic_prioritization(
      curr_wij, pu_purchase_costs, pu_purchase_locked_in,
      pu_purchase_locked_out, obj_fun_target, remaining_budget, curr_solution);

    /// calculate expected value of the action given survey outcome
    curr_expected_value_of_action_given_outcome =
      std::log(expected_value_of_action(
        curr_solution, curr_wij, obj_fun_target));

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

    /// increment o loop variable
    o = o + 1;
  }

  // exports
  return std::exp(curr_expected_value_of_decision);
}
