context("approx sensible values")

test_that("lower voi when most of budget spent on surveys", {
  # initialize rng
  RandomFields::RFoptions(seed = 700)
  set.seed(500)
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
    list(list(objective = "binary:logistic", scale_pos_weight = 1.5,
              nrounds = 8, eta = 0.1))[rep(1, 2)]
  # set approximation values
  n_reps <- 10
  n_states_per_rep <- ceiling(n_states(8, n_f) * 0.9)
  # calculations
  r1 <- approx_evdsi(
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = lapply(xgb_parameters, append, list(nrounds = 8)),
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  r2 <- approx_evdsi(
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = lapply(xgb_parameters, append, list(nrounds = 8)),
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  # tests
  expect_true(all(is.finite(r1)))
  expect_true(all(is.finite(r2)))
  expect_true(all(r1 > 0))
  expect_true(all(r2 > 0))
  expect_true(all(r2 > r1))
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
    list(list(objective = "binary:logistic", scale_pos_weight = 1.5,
              nrounds = 8, eta = 0.1))[rep(1, 2)]
  # set approximation values
  n_reps <- 10
  n_states_per_rep <- ceiling(n_states(8, n_f) * 0.9)
  # calculations
  r1 <- approx_evdsi(
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 5,
    seed = 1,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  r2 <- approx_evdsi(
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  # tests
  expect_true(all(is.finite(r1)))
  expect_true(all(is.finite(r2)))
  expect_true(all(r1 > 0))
  expect_true(all(r2 > 0))
  expect_true(all(r2 > r1))
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
  site_data$w1 <- (runif(nrow(site_data)) * 100) + 1
  site_data$w2 <- (runif(nrow(site_data)) * 100) + 1
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
    list(list(objective = "binary:logistic", scale_pos_weight = 1.5,
              nrounds = 8, eta = 0.1))[rep(1, 2)]
  # set approximation values
  n_reps <- 10
  n_states_per_rep <- min(ceiling(n_states(12, n_f) * 0.9), 1e+4)
  # calculations
  r1 <- approx_evdsi(
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  r2 <- approx_evdsi(
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 1000,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  # tests
  expect_true(all(is.finite(r1)))
  expect_true(all(is.finite(r2)))
  expect_true(all(r1 > 0))
  expect_true(all(r2 > 0))
  expect_gte(max(abs(r2 - r1)), 1e-3)
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
    list(list(objective = "binary:logistic", scale_pos_weight = 1.5,
              nrounds = 8, eta = 0.1))[rep(1, 2)]
  # set approximation values
  n_reps <- 10
  n_states_per_rep <- min(ceiling(n_states(12, n_f) * 0.9), 1e+4)
  # calculations
  r <- lapply(seq_len(20), function(i) {
    approx_evdsi(
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
      feature_preweight_column = "preweight",
      feature_postweight_column = "postweight",
      feature_target_column = "target",
      total_budget = total_budget,
      xgb_parameters = xgb_parameters,
      n_approx_obj_fun_points = 1000,
      xgb_n_folds = rep(5, n_f),
      optimality_gap = 0,
      seed = 1,
      n_approx_replicates = n_reps,
      n_approx_states_per_replicate = n_states_per_rep,
      method_approx_state = "uniform_without_replacement")
  })
  # tests
  for (i in seq_along(r)) {
    expect_true(all(is.finite(r[[i]])))
    expect_true(all(r[[i]] > 0))
    expect_lte(max(abs(r[[1]] - r[[i]])), 1e-10)
  }
})

test_that("current = optimal info = perfect info, when all pu selected", {
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(6),
      y = x,
      f1 = c(0, 1, 1, NA, NA, NA),
      f2 = c(0, 1, 0, NA, NA, NA),
      p1 = c(0.99, 0.99, 0.99, 0.99, 0.99, 0.6),
      p2 = c(0.05, 0.99, 0.05, 0.05, 0.99, 0.4),
      e1 = rnorm(6),
      e2 = rnorm(6),
      survey_cost = c(1, 1, 1, 5, 100000, 1),
      management_cost = c(10, 10, 10, 10, 10, 2),
      locked_in = FALSE),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:2],
    survey = rep(TRUE, 2),
    survey_sensitivity = rep(0.95, 2),
    survey_specificity = rep(0.9, 2),
    model_sensitivity = rep(0.8, 2),
    model_specificity = rep(0.85, 2),
    preweight = 100,
    postweight = 5,
    target = 2)
  xgb_parameters <- list(list(nrounds = 3, eta = 0.3, scale_pos_weight = 1.5,
                              objective = "binary:logistic"))[rep(1, 2)]
  # set approximation values
  n_reps <- 10
  n_states_per_rep <- 20
  # calculate expected values
  evd_current <- approx_evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2"),
    site_probability_columns = c("p1", "p2"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = 100,
    site_management_locked_in_column = "locked_in",
    n_approx_obj_fun_points = 1000,
    optimality_gap = 0,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  evd_ss <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2"),
    site_probability_columns = c("p1", "p2"),
    site_env_vars_columns = c("e1", "e2"),
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = 100,
    survey_budget = 10,
    xgb_parameters = xgb_parameters,
    site_management_locked_in_column = "locked_in",
    n_approx_obj_fun_points = 1000,
    optimality_gap = 0,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  evd_perfect <- approx_evdpi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2"),
    site_probability_columns = c("p1", "p2"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = 100,
    site_management_locked_in_column = "locked_in",
    n_approx_obj_fun_points = 1000,
    optimality_gap = 0,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  # tests
  expect_equal(evd_current, evd_perfect)
  expect_equal(attr(evd_ss, "ev")[1, ], evd_current)
  expect_equal(attr(evd_ss, "ev")[1, ], evd_perfect)
})

test_that("current < optimal info < perfect info, some pu selected", {
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(5),
      y = x,
      f1 = c(0, 1, 1, NA, NA),
      f2 = c(0, 1, 0, NA, NA),
      p1 = c(0.99, 0.99, 0.99, 0.99, 0.6),
      p2 = c(0.05, 0.99, 0.05, 0.05, 0.4),
      e1 = rnorm(5),
      e2 = rnorm(5),
      survey_cost = c(1, 1, 1, 1, 100000),
      management_cost = c(10, 10, 10, 10, 10),
      locked_in = FALSE),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:2],
    survey = rep(TRUE, 2),
    survey_sensitivity = rep(0.95, 2),
    survey_specificity = rep(0.9, 2),
    model_sensitivity = rep(0.8, 2),
    model_specificity = rep(0.85, 2),
    preweight = 100,
    postweight = 5,
    target = 2)
  xgb_parameters <- list(list(nrounds = 3, eta = 0.3, scale_pos_weight = 1.5,
                              objective = "binary:logistic"))[rep(1, 2)]
  budget <- 25
  gap <- 1e-4
  # set approximation values
  n_reps <- 3
  n_states_per_rep <- 1000
  # calculate expected values
  evd_current <- approx_evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2"),
    site_probability_columns = c("p1", "p2"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = budget,
    site_management_locked_in_column = "locked_in",
    n_approx_obj_fun_points = 100,
    optimality_gap = gap,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  evd_ss <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2"),
    site_probability_columns = c("p1", "p2"),
    site_env_vars_columns = c("e1", "e2"),
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = budget,
    survey_budget = 10,
    xgb_parameters = xgb_parameters,
    site_management_locked_in_column = "locked_in",
    n_approx_obj_fun_points = 100,
    optimality_gap = gap,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  evd_perfect <- approx_evdpi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2"),
    site_probability_columns = c("p1", "p2"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = budget,
    site_management_locked_in_column = "locked_in",
    n_approx_obj_fun_points = 100,
    optimality_gap = gap,
    n_approx_replicates = n_reps,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  # tests
  expect_true(all(evd_current < evd_perfect))
  expect_true(all(attr(evd_ss, "ev")[1, ] > evd_current))
  expect_true(all(attr(evd_ss, "ev")[1, ] < evd_perfect))
})

test_that("approx = exact given all states (evdci)", {
  # data
 RandomFields::RFoptions(seed = 700)
  set.seed(500)
  n_f <- 1
  site_data <- simulate_site_data(n_sites = 8, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- sum(site_data$management_cost * 0.8)
  n_states_per_rep <- n_states(nrow(site_data), nrow(feature_data))
  # prepare data
  site_occ_columns <- paste0("f", seq_len(n_f))
  site_prb_columns <- paste0("p", seq_len(n_f))
  site_env_columns <- c("e1", "e2", "e3")
  # calculations
  r1 <- approx_evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    n_approx_obj_fun_points = 1000,
    optimality_gap = 0,
    seed = 1,
    n_approx_replicates = 5,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_states = "uniform_without_replacement")
  r2 <- evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    n_approx_obj_fun_points = 1000,
    optimality_gap = 0)
  expect_lte(max(abs(r1 - r2)), 1e-14)
})

test_that("approx = exact given all states (evdpi)", {
  # data
  RandomFields::RFoptions(seed = 800)
  set.seed(500)
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 5, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- sum(site_data$management_cost) * 0.7
  n_states_per_rep <- n_states(nrow(site_data), nrow(feature_data))
  # prepare data
  site_occ_columns <- paste0("f", seq_len(n_f))
  site_prb_columns <- paste0("p", seq_len(n_f))
  site_env_columns <- c("e1", "e2", "e3")
  # calculations
  r1 <- approx_evdpi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    n_approx_obj_fun_points = 100,
    optimality_gap = 0,
    seed = 1,
    n_approx_replicates = 1,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_states = "uniform_without_replacement")
  r2 <- evdpi(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = site_occ_columns,
    site_probability_columns = site_prb_columns,
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    n_approx_obj_fun_points = 100,
    optimality_gap = 0)
  expect_lte(max(abs(r1 - r2)), 1e-10)
})

test_that("approx = exact, when all states used (evdsi)", {
  # data
  RandomFields::RFoptions(seed = 700)
  set.seed(500)
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 8, n_features = n_f, 0.5)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:2]] <- TRUE
  xgb_n_folds <- rep(5, n_f)
  xgb_parameters <-
    list(list(objective = "binary:logistic", scale_pos_weight = 1.5,
              nrounds = 3))[rep(1, n_f)]
  n_states_per_rep <- n_states(nrow(site_data), nrow(feature_data))
  # prepare data
  site_occ_columns <- paste0("f", seq_len(n_f))
  site_prb_columns <- paste0("p", seq_len(n_f))
  site_env_columns <- c("e1", "e2", "e3")
  # calculations
  r1 <- approx_evdsi(
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 100,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1,
    n_approx_replicates = 5,
    n_approx_states_per_replicate = n_states_per_rep,
    method_approx_state = "uniform_without_replacement")
  r2 <- evdsi(
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_parameters = xgb_parameters,
    n_approx_obj_fun_points = 100,
    xgb_n_folds = rep(5, n_f),
    optimality_gap = 0,
    seed = 1)
  expect_lte(max(abs(r1 - r2)), 1e-11)
})
