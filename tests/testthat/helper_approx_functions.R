r_approx_expected_value_of_decision_given_survey_scheme <- function(
    dij, nij, pij, survey_features, survey_sensitivity, survey_specificity,
    pu_survey_solution, pu_model_prediction,
    pu_survey_costs, pu_purchase_costs,
    pu_purchase_locked_in, pu_purchase_locked_out,
    pu_env_data,
    xgb_parameter_names, xgb_parameter_values,
    n_xgb_rounds, n_xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds,
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
      dij, nij, pij, survey_features, survey_sensitivity, survey_specificity,
      pu_survey_solution, pu_model_prediction,
      pu_survey_costs, pu_purchase_costs,
      pu_purchase_locked_in, pu_purchase_locked_out,
      pu_env_data,
      xgb_parameter_names, xgb_parameter_values,
      n_xgb_rounds, n_xgb_early_stopping_rounds,
      xgb_train_folds, xgb_test_folds,
      obj_fun_target, total_budget, n_approx_outcomes_per_replicate,
      i, seed, outcomes[[i]])
  })
  value
}

r_approx_expected_value_of_decision_given_survey_scheme_fixed_states <-
  function(
    dij, nij, pij, survey_features, survey_sensitivity, survey_specificity,
    pu_survey_solution, pu_model_prediction,
    pu_survey_costs, pu_purchase_costs,
    pu_purchase_locked_in, pu_purchase_locked_out,
    pu_env_data,
    xgb_parameter_names, xgb_parameter_values,
    n_xgb_rounds, n_xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds,
    obj_fun_target, total_budget, n_approx_outcomes_per_replicate,
    r, seed, outcomes) {
  # init
  ## constants
  n_pu <- ncol(dij)
  n_f <- nrow(dij)
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
  ## prepare data for analysis
  oij <- pij
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
  ## prepare updated number of surveys
  m <- as.matrix(expand.grid(
    row = which(survey_features),
    col = which(pu_survey_solution)))
  curr_nij <- nij
  curr_nij[m] <- nij[m] + 1
  rm(m)
  # main processing
  outcome_probs <- numeric(length(outcomes))
  outcome_values <- numeric(length(outcomes))
  for (ii in seq_along(outcomes)) {
    ## generate state
    i <- outcomes[ii]
    curr_oij <- rcpp_nth_state_sparse(i, rij_outcome_idx + 1, oij)

    ## calculate likelihood of outcome
    curr_prob <- probability_of_outcome(
      curr_oij, total_probability_of_survey_positive_log,
      total_probability_of_survey_negative_log, rij_outcome_idx + 1);

    ## update prior data with new survey outcome
    curr_pij <- rcpp_initialize_posterior_probability_matrix(
      nij, pij, curr_oij,
      pu_survey_solution, survey_features,
      survey_sensitivity, survey_specificity)

    ## update survey data with outcome
    curr_dij <- dij
    for (i in (rij_outcome_idx + 1)) {
      curr_dij[i] <-
        ((dij[i] * nij[i]) + curr_oij[i]) / curr_nij[i]
    }

    ## fit distribution models to make predictions
    curr_models <- rcpp_fit_xgboost_models_and_assess_performance(
      curr_dij, curr_nij, curr_pij, pu_env_data,
      as.logical(survey_features), survey_sensitivity, survey_specificity,
      xgb_parameter_names, xgb_parameter_values,
      n_xgb_rounds, n_xgb_early_stopping_rounds,
      xgb_train_folds, xgb_test_folds)

    ## generate model predictions for unsurveyed planning units
    curr_oij2 <- rcpp_predict_missing_rij_data(
      curr_dij, curr_nij, curr_pij, pu_env_data,
      survey_features, survey_sensitivity, survey_specificity,
      xgb_parameter_names, xgb_parameter_values,
      n_xgb_rounds, n_xgb_early_stopping_rounds,
      xgb_train_folds, xgb_test_folds, pu_model_prediction)

    ## extract model sensitivity and specificity data
    curr_models_sens <- rep(NA, n_f)
    curr_models_spec <- rep(NA, n_f)
    curr_models_sens[survey_features_idx] <- curr_models$sens
    curr_models_spec[survey_features_idx] <- curr_models$spec

    ## identify cells in rij matrix with model predictions
    curr_rij_model_idx <-
      r_find_rij_idx_based_on_models(
        nij, pu_survey_solution, survey_features, survey_features_rev_idx,
        survey_sensitivity, survey_specificity,
        curr_models$sens, curr_models$spec) - 1

    ## generate posterior probabilities for most likely model predictions
    curr_postij <- rcpp_update_model_posterior_probabilities(
      nij, pij, curr_oij2,
      pu_survey_solution, survey_features,
      survey_sensitivity, survey_specificity,
      curr_models_sens, curr_models_spec, curr_pij)

    ## generate prioritisation using most likely model predictions
    curr_solution <- r_greedy_heuristic_prioritization(
      curr_postij, pu_purchase_costs,
      as.numeric(pu_purchase_locked_in), as.numeric(pu_purchase_locked_out),
      obj_fun_target, remaining_budget)$x

    ## calculate expected value of decision
    if (length(curr_rij_model_idx) == 0) {
      ### if modelled probabilities were not used to generate the rij matrix,
      ### then we don't have to worry about accounting for them
      curr_value <- log(r_expected_value_of_action(
        curr_solution, curr_postij, obj_fun_target))
    } else {
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
      set.seed(seed + r +  - 1)
      model_outcomes <- rcpp_sample_n_weighted_states_without_replacement(
        curr_n_approx_model_outcomes_per_replicate, curr_model_rij,
        seed + r - 1 + ii - 1)
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
      curr_value <- rcpp_log_sum(curr_probability_of_model_outcome)
    }

    ## store values
    outcome_probs[ii] <- curr_prob
    outcome_values[ii] <- curr_value
  }
  # exports
  outcome_probs <- outcome_probs - rcpp_log_sum(outcome_probs)
  outcome_probs <- outcome_probs + outcome_values
  exp(rcpp_log_sum(outcome_probs))
}
