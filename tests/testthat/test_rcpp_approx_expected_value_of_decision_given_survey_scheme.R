context("rcpp_approx_expected_value_of_decision_given_survey_scheme")

test_that("single species", {
  # data
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_f <- 1
  n_sites <- 20
  n_folds <- 2
  n_vars <- 3
  n_reps <- 3
  n_outcomes_per_rep <- 2
  ## simulate data
  site_data <- simulate_site_data(n_sites, n_f, 0.5, n_vars)
  feature_data <- simulate_feature_data(n_f)
  feature_data$target <- c(10)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(site_data$n1 < 0.5)[1:5]] <- TRUE
  ## create matrices for data
  site_data2 <- sf::st_drop_geometry(site_data)
  dij <- t(as.matrix(site_data2[, paste0("f", seq_len(n_f))], ncol = n_f))
  nij <- t(as.matrix(site_data2[, paste0("n", seq_len(n_f))], ncol = n_f))
  pij <- prior_probability_matrix(
    site_data, feature_data,
    paste0("f", seq_len(n_f)),
    paste0("n", seq_len(n_f)), paste0("p", seq_len(n_f)),
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  ## specify environmental data
  pu_env_data <- as.matrix(site_data2[, paste0("e", seq_len(n_vars))])
  ## specify planning units for predictions
  pu_model_prediction_idx <- lapply(seq_len(n_f), function(i) {
    which((nij[i, ]) < 0.5 & (!site_data$survey))
  })
  ## model fitting parameters
  xgb_folds <- lapply(seq_len(n_f), function(f) {
    fn <- paste0("f", f)
    nn <- paste0("n", f)
    has_data_idx <- which((site_data[[nn]] > 0) | site_data$survey)
    create_site_folds(
      site_data[[fn]][has_data_idx], site_data[[nn]][has_data_idx],
      index = has_data_idx, n_folds)
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  xgb_nrounds <- rep(10, n_f)
  xgb_early_stopping_rounds <- rep(5, n_f)
  tuning_parameters <- expand.grid(
    eta = c(0.1), lambda = c(0.001), objective = "binary:logistic",
    seed = "123")
  tuning_parameters <- as.matrix(tuning_parameters)
  # calculations
  r1 <- rcpp_approx_expected_value_of_decision_given_survey_scheme(
    dij = dij, nij = nij, pij = pij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction_idx,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = pu_env_data,
    xgb_parameter_names = colnames(tuning_parameters),
    xgb_parameter_values = tuning_parameters,
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    method_approx_outcomes = "weighted_without_replacement",
    seed = 500)
  r2 <- r_approx_expected_value_of_decision_given_survey_scheme(
    dij = dij, nij = nij, pij = pij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction_idx,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = pu_env_data,
    xgb_parameter_names = colnames(tuning_parameters),
    xgb_parameter_values = tuning_parameters,
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    seed = 500)
  # tests
  expect_equal(r1, r2)
  expect_false(any(duplicated(r1)))
  expect_false(any(duplicated(r2)))
})

test_that("multiple species", {
  # data
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_f <- 4
  n_sites <- 20
  n_folds <- 2
  n_vars <- 3
  n_reps <- 1
  n_outcomes_per_rep <- 2
  ## simulate data
  site_data <- simulate_site_data(n_sites, n_f, 0.5, n_vars)
  feature_data <- simulate_feature_data(n_f)
  feature_data$target <- c(5)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(site_data$n1 < 0.5)[1:5]] <- TRUE
  ## create matrices for data
  site_data2 <- sf::st_drop_geometry(site_data)
  dij <- t(as.matrix(site_data2[, paste0("f", seq_len(n_f))], ncol = n_f))
  nij <- t(as.matrix(site_data2[, paste0("n", seq_len(n_f))], ncol = n_f))
  pij <- prior_probability_matrix(
    site_data, feature_data,
    paste0("f", seq_len(n_f)),
    paste0("n", seq_len(n_f)), paste0("p", seq_len(n_f)),
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  ## specify environmental data
  pu_env_data <- as.matrix(site_data2[, paste0("e", seq_len(n_vars))])
  ## specify planning units for predictions
  pu_model_prediction_idx <- lapply(seq_len(n_f), function(i) {
    which((nij[i, ]) < 0.5 & (!site_data$survey))
  })
  ## model fitting parameters
  xgb_folds <- lapply(seq_len(n_f), function(f) {
    fn <- paste0("f", f)
    nn <- paste0("n", f)
    has_data_idx <- which((site_data[[nn]] > 0) | site_data$survey)
    create_site_folds(
      site_data[[fn]][has_data_idx], site_data[[nn]][has_data_idx],
      index = has_data_idx, n_folds)
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  xgb_nrounds <- rep(10, n_f)
  xgb_early_stopping_rounds <- rep(5, n_f)
  tuning_parameters <- expand.grid(
    eta = c(0.1), lambda = c(0.001), objective = "binary:logistic",
    seed = "123")
  tuning_parameters <- as.matrix(tuning_parameters)
  # calculations
  r1 <- rcpp_approx_expected_value_of_decision_given_survey_scheme(
    dij = dij, nij = nij, pij = pij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction_idx,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = pu_env_data,
    xgb_parameter_names = colnames(tuning_parameters),
    xgb_parameter_values = tuning_parameters,
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    method_approx_outcomes = "weighted_without_replacement",
    seed = 500)
  r2 <- r_approx_expected_value_of_decision_given_survey_scheme(
    dij = dij, nij = nij, pij = pij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction_idx,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = pu_env_data,
    xgb_parameter_names = colnames(tuning_parameters),
    xgb_parameter_values = tuning_parameters,
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    seed = 500)
  # tests
  expect_equal(r1, r2)
  expect_false(any(duplicated(r1)))
  expect_false(any(duplicated(r2)))
})
