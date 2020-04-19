r_expected_value_of_action <- function(
  solution, prior_data, preweight, postweight, target) {
  # initialization
  print("here")
  print(solution)
  total <- ncol(prior_data)
  sub_prior_data <- prior_data[, solution, drop = FALSE]
  sub_prior_data_log <- log(sub_prior_data)
  sub_prior_data_log1m <- log(1 - sub_prior_data)
  # iterate over each state
  states <- seq(1, rcpp_n_states(length(sub_prior_data)))
  out <- sapply(states, function(i) {
    s <- rcpp_nth_state(i, sub_prior_data)
    v <- r_conservation_benefit_state(s, preweight, postweight, target, total)
    p <- sum(s[] * sub_prior_data_log[]) +
         sum((1 - s[]) * sub_prior_data_log1m[])
    v * exp(p)
  })
  sum(out)
}

r_expected_value_of_decision_given_current_info <- function(
  prior_data, pu_costs, pu_locked_in, preweight, postweight, target,
  n_approx_obj_fun_points, budget, gap) {
  # find optimal solution
  solution <- rcpp_prioritization(
    prior_data, pu_costs, pu_locked_in, preweight, postweight, target,
    n_approx_obj_fun_points, budget, gap, "")$x
  # calculate expected value
  r_expected_value_of_action(
    solution, prior_data, preweight, postweight, target)
}

r_expected_value_of_decision_given_perfect_info <- function(
  prior_data, pu_costs, pu_locked_in, preweight, postweight, target,
  n_approx_obj_fun_points, budget, gap) {
  # calculate log of prior data
  prior_data_log <- log(prior_data)
  prior_data_log1m <- log(1 - prior_data)
  total <- ncol(prior_data)
  # calculate expected value for each state
  states <- seq(1, rcpp_n_states(length(prior_data)))
  out <- sapply(states, function(i) {
    s <- rcpp_nth_state(i, prior_data)
    solution <- r_prioritization(
      s, pu_costs, as.numeric(pu_locked_in), preweight, postweight, target,
      n_approx_obj_fun_points, budget, gap, "")$x
    v <- r_conservation_benefit_state(
      s[, solution, drop = FALSE], preweight, postweight, target, total)
    if (v < 1e-10) return(NA_real_)
    p <- sum(s[] * prior_data_log[]) +
         sum((1 - s[]) * prior_data_log1m[])
    log(v) + p
  })
  out <- out[is.finite(out)]
  exp(rcpp_log_sum(out))
}

r_expected_value_of_decision_given_survey_scheme <- function(
  rij, pij, wij, survey_features, survey_sensitivity, survey_specificity,
  pu_survey_solution, pu_survey_status, pu_survey_costs,
  pu_purchase_costs, pu_purchase_locked_in, pu_env_data,
  xgb_parameters, xgb_nrounds, xgb_train_folds, xgb_test_folds,
  obj_fun_preweight, obj_fun_postweight, obj_fun_target,
  n_approx_obj_fun_points, total_budget, optim_gap) {
  # init
  ## constants
  n_pu <- ncol(rij)
  n_f <- nrow(rij)
  n_f_survey <- sum(survey_features)
  # preliminary processing
  ## calculate remaining budget
  remaining_budget <- total_budget - sum(pu_survey_costs * pu_survey_solution)
  ## calculate total outcomes
  total_outcomes <- n_states(sum(pu_survey_solution), n_f_survey) - 1
  ## planning unit indices
  pu_survey_solution_idx <- which(pu_survey_solution > 0.5)
  pu_survey_status_idx <- which(pu_survey_status > 0.5)
  pu_model_prediction_idx <- which(!pu_survey_status & !pu_survey_solution)
  pu_model_fitting_idx <- which(pu_survey_status | pu_survey_solution)
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
  ## overwrite missing data with prior data for features we are not interested
  ## in surveying
  oij <- rij
  for (i in which(!survey_features))
    oij[i, !pu_survey_status] <- pij[i, !pu_survey_status]
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
  ## create dummy matrix
  dummy_matrix <- matrix(-100, ncol = n_pu, nrow = n_f)
  # main processing
  out <- Inf
  for (i in seq(0, total_outcomes)) {
    ## generate state
    curr_oij <- rcpp_nth_state_sparse(i, rij_outcome_idx + 1, oij)

    ## fit distribution models to make predictions
    curr_models <- rcpp_fit_xgboost_models_and_assess_performance(
      curr_oij, wij, pu_env_data, as.logical(survey_features), xgb_parameters,
      xgb_nrounds, xgb_train_folds, xgb_test_folds)

    ## generate model predictions for unsurveyed planning units
    curr_oij2 <- r_predict_missing_rij_data(
      curr_oij, wij, pu_env_data, survey_features_idx,
      pu_model_prediction_idx, xgb_parameters, xgb_nrounds, xgb_train_folds,
      xgb_test_folds)

    ## generate posterior matrix
    curr_models_sens <- rep(NA, n_f)
    curr_models_spec <- rep(NA, n_f)
    curr_models_sens[survey_features_idx] <- curr_models$sens
    curr_models_spec[survey_features_idx] <- curr_models$spec
    curr_postij <- r_posterior_probability_matrix(
      rij, pij, curr_oij2,
      pu_survey_solution, survey_features,
      survey_sensitivity, survey_specificity,
      curr_models_sens, curr_models_spec)

    ## generate prioritisation
    curr_solution <- r_prioritization(
      curr_postij, pu_purchase_costs, as.numeric(pu_purchase_locked_in),
      obj_fun_preweight, obj_fun_postweight, obj_fun_target,
      n_approx_obj_fun_points, remaining_budget, optim_gap, "")$x

    ## calculate expected value of the prioritisation
    curr_value <- log(r_expected_value_of_action(
      curr_solution, curr_postij, obj_fun_preweight, obj_fun_postweight,
      obj_fun_target))

    ## calculate likelihood of outcome
    curr_prob <- probability_of_outcome(
      curr_oij2, total_probability_of_survey_positive_log,
      total_probability_of_survey_negative_log, rij_outcome_idx + 1);

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
