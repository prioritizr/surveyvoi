r_expected_value_of_action <- function(
  solution, prior_data, target) {
    r_conservation_value(prior_data[, solution > 0.5, drop = FALSE], target)
}

r_expected_value_of_decision_given_current_info <- function(
  prior_data, pu_costs, pu_purchase_locked_in, pu_purchase_locked_out, target, budget) {
  # find optimal solution
  solution <- r_greedy_heuristic_prioritization(
    prior_data, pu_costs, as.numeric(pu_purchase_locked_in),
    as.numeric(pu_purchase_locked_out), target, budget)$x
  # calculate expected value
  r_expected_value_of_action(
    solution, prior_data, target)
}

r_expected_value_of_decision_given_survey_scheme <- function(
  dij, nij, pij,
  survey_features,
  survey_sensitivity, survey_specificity,
  pu_survey_solution, pu_model_prediction, pu_survey_costs,
  pu_purchase_costs, pu_purchase_locked_in, pu_purchase_locked_out, pu_env_data,
  xgb_parameter_names, xgb_parameter_values,
  n_xgb_rounds, n_xgb_early_stopping_rounds,
  xgb_train_folds, xgb_test_folds,
  obj_fun_target, total_budget) {
  # init
  ## constants
  n_pu <- ncol(dij)
  n_f <- nrow(dij)
  n_f_survey <- sum(survey_features)
  # preliminary processing
  ## calculate remaining budget
  remaining_budget <- total_budget - sum(pu_survey_costs * pu_survey_solution)
  ## calculate total outcomes
  total_outcomes <- n_states(sum(pu_survey_solution), n_f_survey) - 1
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
  out <- Inf
  for (i in seq(0, total_outcomes)) {
    ## generate state
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
      curr_value <- exact_expected_value_given_survey_outcome(
        nij,
        pij,
        curr_postij,
        curr_solution,
        pu_survey_solution,
        survey_features,
        curr_rij_model_idx,
        survey_sensitivity, survey_specificity,
        curr_models_sens, curr_models_spec,
        obj_fun_target)
    }
    ## calculate expected value of the action
    if (!is.finite(out)) {
      out <- curr_value + curr_prob
    } else {
      out <- log_sum(out, curr_value + curr_prob)
    }
  }
  # exports
  exp(out)
}
