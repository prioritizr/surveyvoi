context("sensible values")

test_that("lower voi when most of budget spent on surveys", {
  # initialize rng
  set.seed(501)
  RandomFields::RFoptions(seed = 501)
  # data
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 8, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- 5
  # manually define costs
  site_data$management_cost <- 0.1
  site_data$survey_cost <- c(4.8, rep(0.01, 7))
  # create surveys
  site_data$survey1 <- FALSE
  site_data$survey1[1] <- TRUE
  site_data$survey2 <- FALSE
  site_data$survey2[c(3, 6)] <- TRUE
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  site_env_columns <- c("e1", "e2", "e3")
  # prepare xgboost inputs
  xgb_n_folds <- rep(5, n_f)
  xgb_parameters <-
    list(list(objective = "binary:logistic", nrounds = 8, eta = 0.1))[rep(1, 2)]
  # calculations
  r1 <- evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_survey_scheme_column = "survey1",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_alpha_column = "alpha",
    feature_gamma_column = "gamma",
    total_budget = total_budget,
    xgb_parameters = lapply(xgb_parameters, append, list(nrounds = 8)),
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1)
  r2 <- evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_survey_scheme_column = "survey2",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_alpha_column = "alpha",
    feature_gamma_column = "gamma",
    total_budget = total_budget,
    xgb_parameters = lapply(xgb_parameters, append, list(nrounds = 8)),
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1)
  # tests
  expect_true(is.finite(r1))
  expect_true(is.finite(r2))
  expect_gt(r1, 0)
  expect_gt(r2, 0)
  expect_gt(r2, r1)
})

test_that("larger optimality gap produces lower voi of survey scheme", {
  # initialize rng
  set.seed(501)
  RandomFields::RFoptions(seed = 501)
  # data
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 8, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- 5
  # manually define costs
  site_data$management_cost <- 0.1
  site_data$survey_cost <- c(4.8, rep(0.01, 7))
  # create surveys
  site_data$survey1 <- FALSE
  site_data$survey1[1] <- TRUE
  site_data$survey2 <- FALSE
  site_data$survey2[c(3, 6)] <- TRUE
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  site_env_columns <- c("e1", "e2", "e3")
  # prepare xgboost inputs
  xgb_n_folds <- rep(5, n_f)
  xgb_parameters <-
    list(list(objective = "binary:logistic", nrounds = 8, eta = 0.1))[rep(1, 2)]
  # calculations
  r1 <- evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_survey_scheme_column = "survey1",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_alpha_column = "alpha",
    feature_gamma_column = "gamma",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 100,
    seed = 1)
  r2 <- evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_survey_scheme_column = "survey2",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_alpha_column = "alpha",
    feature_gamma_column = "gamma",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1)
  # tests
  expect_true(is.finite(r1))
  expect_true(is.finite(r2))
  expect_gt(r1, 0)
  expect_gt(r2, 0)
  expect_gt(r2, r1)
})

test_that("different voi when xgboost models trained with different weights", {
  # initialize rng
  set.seed(501)
  RandomFields::RFoptions(seed = 501)
  # data
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 12, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$w1 <- runif(nrow(site_data)) + 1
  site_data$w2 <- runif(nrow(site_data)) + 1
  # create survey scheme
  site_data$survey_scheme <- FALSE
  site_data$survey_scheme[which(is.na(site_data$f1))[1:2]] <- TRUE
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  site_env_columns <- c("e1", "e2", "e3")
  site_wgt_columns <- c("w1", "w2")
  # prepare xgboost inputs
  xgb_n_folds <- rep(5, n_f)
  xgb_parameters <-
    list(list(objective = "binary:logistic", nrounds = 8, eta = 0.1))[rep(1, 2)]
  # calculations
  r1 <- evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_survey_scheme_column = "survey_scheme",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_alpha_column = "alpha",
    feature_gamma_column = "gamma",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1)
  r2 <- evdsi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_env_vars_columns = site_env_columns,
    site_weight_columns = site_wgt_columns,
    site_survey_scheme_column = "survey_scheme",
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_alpha_column = "alpha",
    feature_gamma_column = "gamma",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1)
  # tests
  expect_true(is.finite(r1))
  expect_true(is.finite(r2))
  expect_gt(r1, 0)
  expect_gt(r2, 0)
  expect_gte(abs(r2 - r1), 1e-3)
})

test_that("identical outputs given identical inputs", {
  # initialize rng
  set.seed(501)
  RandomFields::RFoptions(seed = 501)
  # data
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 8, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- 5
  # manually define costs
  site_data$management_cost <- 0.1
  site_data$survey_cost <- c(4.8, rep(0.01, 7))
  # create surveys
  site_data$survey <- FALSE
  site_data$survey[c(3, 6)] <- TRUE
  # prepare data
  site_occ_columns <- c("f1", "f2")
  site_prb_columns <- c("p1", "p2")
  site_env_columns <- c("e1", "e2", "e3")
  # prepare xgboost inputs
  xgb_n_folds <- rep(5, n_f)
  xgb_parameters <-
    list(list(objective = "binary:logistic", nrounds = 8, eta = 0.1))[rep(1, 2)]
  # calculations
  r <- vapply(seq_len(20), FUN.VALUE = numeric(1), function(i) {
    evdsi(
      site_data = site_data,
      feature_data = feature_data,
      site_occupancy_columns = site_occ_columns,
      site_probability_columns = site_prb_columns,
      site_env_vars_columns = site_env_columns,
      site_survey_scheme_column = "survey",
      site_management_cost_column = "management_cost",
      site_survey_cost_column = "survey_cost",
      feature_survey_column = "survey",
      feature_survey_sensitivity_column = "survey_sensitivity",
      feature_survey_specificity_column = "survey_specificity",
      feature_model_sensitivity_column = "model_sensitivity",
      feature_model_specificity_column = "model_specificity",
      feature_alpha_column = "alpha",
      feature_gamma_column = "gamma",
      total_budget = total_budget,
      xgb_parameters = lapply(xgb_parameters, append, list(nrounds = 8)),
      n_approx_obj_fun_points = 1000,
      xgb_n_folds = rep(5, n_f),
      optimality_gap = 0,
      seed = 1)
  })
  # tests
  expect_true(all(is.finite(r)))
  expect_true(all(r > 0))
  expect_lte(abs(diff(range(r))), 1e-10)
})

test_that(
  "optimal_survey_scheme/expected_value_of_survey_scheme same results", {
})
