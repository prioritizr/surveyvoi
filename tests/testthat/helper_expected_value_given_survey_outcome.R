exact_expected_value_given_survey_outcome <- function(
  nij,
  pij,
  curr_postij,
  curr_solution,
  pu_survey_solution,
  survey_features,
  curr_rij_model_idx,
  survey_sensitivity, survey_specificity,
  curr_models_sens, curr_models_spec,
  obj_fun_target) {
  ### if modelled probabilities were used to generate the rij matrix,
  ### then we need to account for them
  #### generate model outcome
  n_model_outcomes <- n_states(length(curr_rij_model_idx), 1) - 1
  #### calculate log probabilities of model returning positive or negative
  curr_total_probability_of_model_positive_log <-
    log(total_probability_of_positive_result(
      pij, curr_models_sens, curr_models_spec))
  curr_total_probability_of_model_negative_log <-
    log(total_probability_of_negative_result(
      pij, curr_models_sens, curr_models_spec))
  #### initialize looping variables
  curr_expected_value_of_action_given_outcome <- Inf
  #### iterate over each model outcome
  for (ii in seq(0, n_model_outcomes)) {
    ##### generate model state
    curr_postij2 <- rcpp_nth_state_sparse(
      ii, curr_rij_model_idx + 1, curr_postij)
    ##### calculate probability of model state
    curr_probability_of_model_outcome <-
      probability_of_outcome(
        curr_postij2,
        curr_total_probability_of_model_positive_log,
        curr_total_probability_of_model_negative_log,
        curr_rij_model_idx + 1)
    ##### generate posterior probability given outcome
    curr_postij2 <- rcpp_update_model_posterior_probabilities(
      nij, pij, curr_postij2,
      pu_survey_solution, survey_features,
      survey_sensitivity, survey_specificity,
      curr_models_sens, curr_models_spec, curr_postij2)
    ##### calculate expected value of action given model outcome
    curr_expected_value_of_action_given_model_outcome <-
      log(r_expected_value_of_action(
        curr_solution, curr_postij2, obj_fun_target))
    #### add conditional value to running total across all models states
    if (!is.finite(curr_expected_value_of_action_given_outcome)) {
      curr_expected_value_of_action_given_outcome <-
        curr_expected_value_of_action_given_model_outcome +
        curr_probability_of_model_outcome
    } else {
      curr_expected_value_of_action_given_outcome <-
        log_sum(curr_expected_value_of_action_given_outcome,
                curr_expected_value_of_action_given_model_outcome +
                curr_probability_of_model_outcome)
    }
  }
  curr_expected_value_of_action_given_outcome
}


approx_expected_value_given_survey_outcome <- function(
  nij,
  pij,
  curr_postij,
  curr_solution,
  pu_survey_solution,
  survey_features,
  curr_rij_model_idx,
  survey_sensitivity, survey_specificity,
  curr_models_sens, curr_models_spec,
  obj_fun_target,
  n_approx_outcomes_per_replicate,
  seed) {
  #### determine number of model outcomes
  curr_n_approx_model_outcomes_per_replicate <-
    n_approx_outcomes_per_replicate;
  if (length(curr_rij_model_idx) < 30) {
    curr_n <- rcpp_n_states(length(curr_rij_model_idx)) + 1
    if (curr_n < curr_n_approx_model_outcomes_per_replicate)
      curr_n_approx_model_outcomes_per_replicate <- curr_n
  }
  #### calculate log probabilities of model returning positive or negative
  curr_total_probability_of_model_positive_log <-
    log(total_probability_of_positive_result(
      pij, curr_models_sens, curr_models_spec))
  curr_total_probability_of_model_negative_log <-
    log(total_probability_of_negative_result(
      pij, curr_models_sens, curr_models_spec))
  #### prepare rij data for generating model outcomes
  curr_model_rij <- matrix(curr_postij[curr_rij_model_idx + 1], ncol = 1)
  #### generate model outcomes
  model_outcomes <- rcpp_sample_n_weighted_states_without_replacement(
    curr_n_approx_model_outcomes_per_replicate, curr_model_rij, seed)
  #### initialize looping variables
  curr_expected_value_of_action_given_model_outcome <-
    numeric(curr_n_approx_model_outcomes_per_replicate)
  curr_probability_of_model_outcome <-
    numeric(curr_n_approx_model_outcomes_per_replicate)
  #### iterate over each model outcome and calculate results
  for (oo in seq_along(model_outcomes)) {
    ##### generate model state
    curr_postij2 <- rcpp_nth_state_sparse(
      model_outcomes[oo], curr_rij_model_idx + 1, curr_postij)
    ##### calculate probability of model state
    curr_probability_of_model_outcome[oo] <-
      probability_of_outcome(
        curr_postij2,
        curr_total_probability_of_model_positive_log,
        curr_total_probability_of_model_negative_log,
        curr_rij_model_idx + 1)
    ##### generate posterior probability given outcome
    curr_postij2 <- rcpp_update_model_posterior_probabilities(
      nij, pij, curr_postij2,
      pu_survey_solution, survey_features,
      survey_sensitivity, survey_specificity,
      curr_models_sens, curr_models_spec, curr_postij2)
    ##### calculate expected value of action given model outcome
    curr_expected_value_of_action_given_model_outcome[oo] <-
      log(r_expected_value_of_action(
        curr_solution, curr_postij2, obj_fun_target))
  }
  #### store value of decision for the o'th survey outcome
  #### given uncertainty in the model outcomes
  curr_probability_of_model_outcome <-
    curr_probability_of_model_outcome -
    rcpp_log_sum(curr_probability_of_model_outcome)
  curr_probability_of_model_outcome <-
    curr_probability_of_model_outcome +
    curr_expected_value_of_action_given_model_outcome
  rcpp_log_sum(curr_probability_of_model_outcome)
}
