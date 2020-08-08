context("rcpp_expected_value_of_decision_given_survey_scheme")

test_that("equal weights", {
  # set seeds
  set.seed(500)
  RandomFields::RFoptions(seed = 500)
  # data
  n_f <- 2
  n_sites <- 30
  n_folds <- 5
  site_data <- simulate_site_data(n_sites, n_f, 0.5)
  feature_data <- simulate_feature_data(n_f, 0.5)
  feature_data$target <- c(15, 15)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:2]] <- TRUE
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  env_columns <- c("e1", "e2")
  rij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_occ_columns])))
  pij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_prb_columns])))
  wij <- matrix(1, ncol = ncol(pij), nrow = nrow(pij))
  ejx <- as.matrix(sf::st_drop_geometry(site_data[, env_columns]))
  pu_model_prediction <- lapply(seq_len(nrow(feature_data)), function(i) {
    which(!site_data$survey & is.na(site_data[[paste0("f", i)]]))
  })
  ## set NA rij values to -1
  rij[is.na(rij)] <- -1
  ## model fitting parameters
  xgb_folds <- lapply(paste0("f", seq_len(n_f)), function(f) {
    non_na_idx <- which(!(is.na(site_data[[f]] & !site_data$survey)))
    create_folds(site_data[[f]][non_na_idx], index = non_na_idx, n_folds,
                 na.fail = FALSE)
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  ## set xgboost modelling parameters
  xgb_nrounds <- rep(10, n_f)
  xgb_early_stopping_rounds <- rep(5, n_f)
  tuning_parameters <-
    expand.grid(eta = c(0.1, 0.5, 1.0),
                lambda = c(0.001, 0.01, 0.05),
                objective = "binary:logistic",
                seed = "123")
  tuning_parameters <- as.matrix(tuning_parameters)
  # calculations
  r1 <- rcpp_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameter_names = colnames(tuning_parameters),
    xgb_parameter_values = tuning_parameters,
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    optim_gap = 0)
  r2 <- r_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameter_names = colnames(tuning_parameters),
    xgb_parameter_values = tuning_parameters,
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget)
  # tests
  expect_equal(r1, r2)
})

test_that("variables weights", {
  # set seeds
  set.seed(500)
  RandomFields::RFoptions(seed = 500)
  # data
  n_f <- 2
  n_sites <- 30
  n_folds <- 5
  site_data <- simulate_site_data(n_sites, n_f, 0.5)
  feature_data <- simulate_feature_data(n_f, 0.5)
  feature_data$target <- c(15, 15)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:2]] <- TRUE
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  env_columns <- c("e1", "e2")
  rij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_occ_columns])))
  pij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_prb_columns])))
  wij <- matrix(runif(length(rij)) + 0.1, ncol = ncol(pij), nrow = nrow(pij))
  ejx <- as.matrix(sf::st_drop_geometry(site_data[, env_columns]))
  pu_model_prediction <- lapply(seq_len(nrow(feature_data)), function(i) {
    which(!site_data$survey & is.na(site_data[[paste0("f", i)]]))
  })
  ## model fitting parameters
  xgb_folds <- lapply(paste0("f", seq_len(n_f)), function(f) {
    non_na_idx <- which(!(is.na(site_data[[f]] & !site_data$survey)))
    create_folds(site_data[[f]][non_na_idx], index = non_na_idx, n_folds,
                 na.fail = FALSE)
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  ## set xgboost modelling parameters
  xgb_nrounds <- rep(10, n_f)
  xgb_early_stopping_rounds <- rep(5, n_f)
  tuning_parameters <-
    expand.grid(eta = c(0.1, 0.5, 1.0),
                lambda = c(0.001, 0.01, 0.05),
                objective = "binary:logistic",
                seed = "123")
  tuning_parameters <- as.matrix(tuning_parameters)
  ## set NA rij values to -1
  rij[is.na(rij)] <- -1
  # calculations
  r1 <- rcpp_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameter_names = colnames(tuning_parameters),
    xgb_parameter_values = tuning_parameters,
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    optim_gap = 0)
  r2 <- r_expected_value_of_decision_given_survey_scheme(
    rij = rij, pij = pij, wij = wij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_model_prediction = pu_model_prediction,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    pu_env_data = ejx,
    xgb_parameter_names = colnames(tuning_parameters),
    xgb_parameter_values = tuning_parameters,
    n_xgb_rounds = rep(10, n_f),
    n_xgb_early_stopping_rounds = rep(5, n_f),
    xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
    xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
    obj_fun_target = feature_data$target,
    total_budget = total_budget)
  # tests
  expect_equal(r1, r2)
})
