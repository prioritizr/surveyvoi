r_approx_expected_value_of_decision_given_survey_scheme <- function(
    pij,
    survey_features, survey_sensitivity, survey_specificity,
    pu_survey_solution,
    pu_survey_costs, pu_purchase_costs,
    pu_purchase_locked_in, pu_purchase_locked_out,
    obj_fun_target, total_budget,
    n_approx_replicates, n_approx_outcomes_per_replicate, seed) {
  # generate outcomes
  outcomes <- lapply(seq_len(n_approx_replicates), function(i) {
    set.seed(seed + i - 1)
    rcpp_sample_n_weighted_states_without_replacement(
      n_approx_outcomes_per_replicate,
      pij[survey_features, pu_survey_solution, drop = FALSE],
      seed + i - 1)
  })
  set.seed(seed)

  # run calculations
  value <- sapply(seq_len(n_approx_replicates), function(i) {
    r_approx_expected_value_of_decision_given_survey_scheme_fixed_states(
      pij,
      survey_features, survey_sensitivity, survey_specificity,
      pu_survey_solution,
      pu_survey_costs, pu_purchase_costs,
      pu_purchase_locked_in, pu_purchase_locked_out,
      obj_fun_target, total_budget, outcomes[[i]])
  })
  value
}

r_approx_expected_value_of_decision_given_survey_scheme_fixed_states <-
  function(
    pij,
    survey_features, survey_sensitivity, survey_specificity,
    pu_survey_solution,
    pu_survey_costs, pu_purchase_costs,
    pu_purchase_locked_in, pu_purchase_locked_out,
    obj_fun_target, total_budget, outcomes) {
  # init
  ## constants
  n_pu <- ncol(pij)
  n_f <- nrow(pij)
  n_f_survey <- sum(survey_features)
  # preliminary processing
  ## calculate remaining budget
  remaining_budget <- total_budget - sum(pu_survey_costs * pu_survey_solution)
  ## planning unit indices
  pu_survey_solution_idx <- which(pu_survey_solution > 0.5)
  ## feature indices
  survey_features_idx <- which(survey_features > 0.5)
  survey_features_rev_idx <- rep(0, n_f)
  survey_features_rev_idx[survey_features_idx] <-
    seq_along(survey_features_idx) - 1
  ## subset of prior matrix for features that need surveying
  pij_survey_species_subset <- pij[survey_features_idx, , drop = FALSE]
  ## total probability of positive survey result
  total_probability_of_survey_positive <-
    total_probability_of_positive_result(
      pij, survey_sensitivity, survey_specificity)
  total_probability_of_survey_positive_log <-
    total_probability_of_survey_positive
  total_probability_of_survey_positive_log[] <-
    log(total_probability_of_survey_positive)
  ## total probability of negative survey result
  total_probability_of_survey_negative <-
    total_probability_of_negative_result(
      pij, survey_sensitivity, survey_specificity)
  total_probability_of_survey_negative_log <-
    total_probability_of_survey_negative
  total_probability_of_survey_negative_log[] <-
    log(total_probability_of_survey_negative)
  ## find indices for cells corresponding to planning units and features that
  ## that are specified to be surveyed
  rij_outcome_idx <- c()
  counter <- 0
  for (j in seq_len(n_pu)) {
    for (i in seq_len(n_f)) {
      if ((survey_features[i] > 0.5) && (pu_survey_solution[j] > 0.5))
        rij_outcome_idx <- c(rij_outcome_idx, counter)
      counter <- counter + 1
    }
  }

  # main processing
  outcome_probs <- numeric(length(outcomes))
  outcome_values <- numeric(length(outcomes))
  for (ii in seq_along(outcomes)) {
    ## generate state
    i <- outcomes[ii]
    curr_oij <- rcpp_nth_state_sparse(i, rij_outcome_idx + 1, pij)

    ## calculate likelihood of outcome
    curr_prob <- probability_of_outcome(
      curr_oij, total_probability_of_survey_positive_log,
      total_probability_of_survey_negative_log, rij_outcome_idx + 1);

    ## update prior data with new survey outcome
    curr_wij <- r_update_posterior_probability_matrix(
      pij, curr_oij,
      survey_features, survey_sensitivity, survey_specificity,
      pu_survey_solution)

    ## generate prioritisation using most likely model predictions
    curr_solution <- r_greedy_heuristic_prioritization(
      curr_wij, pu_purchase_costs,
      as.numeric(pu_purchase_locked_in), as.numeric(pu_purchase_locked_out),
      obj_fun_target, remaining_budget)$x

    ## calculate expected value of action
    curr_value <- log(r_expected_value_of_action(
      curr_solution, curr_wij, obj_fun_target))

    ## store values
    outcome_probs[ii] <- curr_prob
    outcome_values[ii] <- curr_value
  }
  # exports
  outcome_probs <- outcome_probs - rcpp_log_sum(outcome_probs)
  outcome_probs <- outcome_probs + outcome_values
  exp(rcpp_log_sum(outcome_probs))
}
