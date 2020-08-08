context("approx_evdsi")

test_that("equal weights", {
  # data
  RandomFields::RFoptions(seed = 501)
  set.seed(500)
  n_f <- 2
  n_sites <- 30
  n_reps <- 3
  n_outcomes_per_rep <- 4
  site_data <- simulate_site_data(n_sites, n_f, 0.5)
  feature_data <- simulate_feature_data(n_f, 0.5)
  feature_data$target <- c(15, 15)
  total_budget <- sum(site_data$management_cost * 0.9)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:5]] <- TRUE
  site_data$w1 <- 1
  site_data$w2 <- 1
  xgb_n_folds <- rep(5, n_f)
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  site_env_columns <- c("e1", "e2", "e3")
  site_wgt_columns <- c("w1", "w2")
  rij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_occ_columns])))
  pij <- prior_probability_matrix(
    site_data, feature_data, site_occ_columns, site_prb_columns,
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  wij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_wgt_columns])))
  ejx <- as.matrix(sf::st_drop_geometry(site_data[, site_env_columns]))
  pu_model_prediction <- lapply(seq_len(nrow(feature_data)), function(i) {
    which(!site_data$survey & is.na(rij[i, ]))
  })
  # prepare xgboost inputs
  xgb_tuning_parameters <-
    list(objective = "binary:logistic", lambda = c(0.01, 0.1, 0.5))
  xgb_full_parameters <- do.call(expand.grid, xgb_tuning_parameters)
  attr(xgb_full_parameters, "out.attrs") <- NULL
  xgb_full_parameters$seed <- as.character(1) # set seed
  ## folds for training and testing models
  xgb_folds <- lapply(seq_len(nrow(feature_data)), function(i) {
    pu_train_idx <- which(site_data$survey | !is.na(rij[i, ]))
    withr::with_seed(seed,
      create_folds(unname(rij[i, pu_train_idx]), xgb_n_folds[i],
                   index = pu_train_idx,
                   na.fail = FALSE,
                   seed = 1))
  })
  ## set NA rij values to -1
  rij[is.na(rij)] <- -1
  # calculations
  r1 <- approx_evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_weight_columns = site_wgt_columns,
    site_survey_scheme_column = "survey",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_tuning_parameters = xgb_tuning_parameters,
    xgb_early_stopping_rounds = rep(5, n_f),
    xgb_n_rounds = rep(10, n_f),
    xgb_n_folds = xgb_n_folds,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    method_approx_outcomes = "weighted_without_replacement",
    seed = 1)
  r2 <- r_approx_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction,
    pu_survey_costs = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameter_names = names(xgb_full_parameters),
    xgb_parameter_values = as.matrix(xgb_full_parameters),
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    seed = 1)
  expect_equal(r1, r2)
  expect_false(any(duplicated(r1)))
})

test_that("variable weights", {
  # data
  RandomFields::RFoptions(seed = 501)
  set.seed(500)
  n_f <- 2
  n_sites <- 30
  n_reps <- 3
  n_outcomes_per_rep <- 4
  site_data <- simulate_site_data(n_sites, n_f, 0.5)
  feature_data <- simulate_feature_data(n_f, 0.5)
  feature_data$target <- c(15, 15)
  total_budget <- sum(site_data$management_cost * 0.9)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:5]] <- TRUE
  site_data$w1 <- runif(nrow(site_data)) + 1
  site_data$w2 <- runif(nrow(site_data)) + 1
  xgb_n_folds <- rep(5, n_f)
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  site_env_columns <- c("e1", "e2", "e3")
  site_wgt_columns <- c("w1", "w2")
  rij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_occ_columns])))
  pij <- prior_probability_matrix(
    site_data, feature_data, site_occ_columns, site_prb_columns,
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  wij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_wgt_columns])))
  ejx <- as.matrix(sf::st_drop_geometry(site_data[, site_env_columns]))
  pu_model_prediction <- lapply(seq_len(nrow(feature_data)), function(i) {
    which(!site_data$survey & is.na(rij[i, ]))
  })
  # prepare xgboost inputs
  xgb_tuning_parameters <-
    list(objective = "binary:logistic", lambda = c(0.01, 0.1, 0.5))
  xgb_full_parameters <- do.call(expand.grid, xgb_tuning_parameters)
  attr(xgb_full_parameters, "out.attrs") <- NULL
  xgb_full_parameters$seed <- as.character(1) # set seed
  ## folds for training and testing models
  xgb_folds <- lapply(seq_len(nrow(feature_data)), function(i) {
    pu_train_idx <- which(site_data$survey | !is.na(rij[i, ]))
    withr::with_seed(seed,
      create_folds(unname(rij[i, pu_train_idx]), xgb_n_folds[i],
                   index = pu_train_idx,
                   na.fail = FALSE,
                   seed = 1))
  })
  ## set NA rij values to -1
  rij[is.na(rij)] <- -1
  # calculations
  r1 <- approx_evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_weight_columns = site_wgt_columns,
    site_survey_scheme_column = "survey",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_tuning_parameters = xgb_tuning_parameters,
    xgb_early_stopping_rounds = rep(5, n_f),
    xgb_n_rounds = rep(10, n_f),
    xgb_n_folds = xgb_n_folds,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    method_approx_outcomes = "weighted_without_replacement",
    seed = 1)
  r2 <- r_approx_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction,
    pu_survey_costs = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameter_names = names(xgb_full_parameters),
    xgb_parameter_values = as.matrix(xgb_full_parameters),
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    seed = 1)
  expect_equal(r1, r2)
  expect_false(any(duplicated(r1)))
})

test_that("sparse", {
  # data
  RandomFields::RFoptions(seed = 501)
  set.seed(500)
  n_f <- 2
  n_sites <- 30
  n_reps <- 3
  n_outcomes_per_rep <- 4
  site_data <- simulate_site_data(n_sites, n_f, 0.5)
  feature_data <- simulate_feature_data(n_f, 0.5)
  feature_data$target <- c(15, 15)
  total_budget <- sum(site_data$management_cost * 0.9)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:5]] <- TRUE
  site_data$w1 <- runif(nrow(site_data)) + 1
  site_data$w2 <- runif(nrow(site_data)) + 1
  xgb_n_folds <- rep(5, n_f)
  # randomly add sparsity to survey data
  for (i in paste0("f", seq_len(n_f))) {
    na_idx <- which(is.na(site_data[[i]]))
    idx <- sample(na_idx, ceiling(length(na_idx) * 0.5))
    site_data[[i]][idx] <- round(runif(length(idx)))
  }
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  site_env_columns <- c("e1", "e2", "e3")
  site_wgt_columns <- c("w1", "w2")
  rij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_occ_columns])))
  pij <- prior_probability_matrix(
    site_data, feature_data, site_occ_columns, site_prb_columns,
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  wij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_wgt_columns])))
  ejx <- as.matrix(sf::st_drop_geometry(site_data[, site_env_columns]))
  pu_model_prediction <- lapply(seq_len(nrow(feature_data)), function(i) {
    which(!site_data$survey & is.na(rij[i, ]))
  })
  # prepare xgboost inputs
  xgb_tuning_parameters <-
    list(objective = "binary:logistic", lambda = c(0.01, 0.1, 0.5))
  xgb_full_parameters <- do.call(expand.grid, xgb_tuning_parameters)
  attr(xgb_full_parameters, "out.attrs") <- NULL
  xgb_full_parameters$seed <- as.character(1) # set seed
  ## folds for training and testing models
  xgb_folds <- lapply(seq_len(nrow(feature_data)), function(i) {
    pu_train_idx <- which(site_data$survey | !is.na(rij[i, ]))
    withr::with_seed(seed,
      create_folds(unname(rij[i, pu_train_idx]), xgb_n_folds[i],
                   index = pu_train_idx,
                   na.fail = FALSE,
                   seed = 1))
  })
  ## set NA rij values to -1
  rij[is.na(rij)] <- -1
  # calculations
  r1 <- approx_evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_weight_columns = site_wgt_columns,
    site_survey_scheme_column = "survey",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_tuning_parameters = xgb_tuning_parameters,
    xgb_early_stopping_rounds = rep(5, n_f),
    xgb_n_rounds = rep(10, n_f),
    xgb_n_folds = xgb_n_folds,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    method_approx_outcomes = "weighted_without_replacement",
    seed = 1)
  r2 <- r_approx_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction,
    pu_survey_costs = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameter_names = names(xgb_full_parameters),
    xgb_parameter_values = as.matrix(xgb_full_parameters),
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    seed = 1)
  expect_equal(r1, r2)
  expect_false(any(duplicated(r1)))
})
