context("approx sensible values")

test_that("lower voi when most of budget spent on surveys", {
  # initialize rng
  set.seed(501)
  RandomFields::RFoptions(seed = 501)
  # data
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 15, n_features = n_f, 0.1)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  feature_data$target <- c(2, 2)
  total_budget <- 5
  # manually define costs
  site_data$management_cost <- 0.1
  site_data$survey_cost <- c(4.8, rep(0.01, nrow(site_data) - 1))
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
  xgb_tuning_parameters <-
    list(objective = "binary:logistic", lambda = c(0.01, 0.1, 0.5))
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
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_tuning_parameters = xgb_tuning_parameters,
    xgb_n_folds = rep(2, n_f),
    seed = 1)
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
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_tuning_parameters = xgb_tuning_parameters,
    xgb_n_folds = rep(2, n_f),
    seed = 1)
  # tests
  expect_true(all(is.finite(r1)))
  expect_true(all(is.finite(r2)))
  expect_true(all(r1 >= 0))
  expect_true(all(r2 >= 0))
  expect_true(all(r1 <= 1))
  expect_true(all(r2 <= 1))
  expect_true(all(r2 >= r1))
})

test_that("current == optimal info, when all pu selected", {
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(7),
      y = x,
      f1 = c(0, 1, 1, 0, NA, NA, NA),
      f2 = c(0, 1, 0, 1, NA, NA, NA),
      p1 = c(0.99, 0.99, 0.99, 0.05, 0.99, 0.99, 0.6),
      p2 = c(0.05, 0.99, 0.05, 0.99, 0.05, 0.99, 0.4),
      e1 = rnorm(7),
      e2 = rnorm(7),
      survey_cost = c(1, 1, 1, 1, 5, 100000, 1),
      management_cost = c(10, 10, 10, 10, 10, 10, 2),
      locked_in = FALSE),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:2],
    survey = rep(TRUE, 2),
    survey_sensitivity = rep(0.95, 2),
    survey_specificity = rep(0.9, 2),
    model_sensitivity = rep(0.8, 2),
    model_specificity = rep(0.85, 2),
    target = c(1, 1))
  xgb_tuning_parameters <-
    list(objective = "binary:logistic", lambda = c(0.01, 0.1, 0.5))
  # calculate expected values
  evd_current <- evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2"),
    site_probability_columns = c("p1", "p2"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = 100,
    site_management_locked_in_column = "locked_in")
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
    feature_target_column = "target",
    total_budget = 100,
    survey_budget = 10,
    xgb_n_fold = rep(2, nrow(feature_data)),
    xgb_tuning_parameters = xgb_tuning_parameters,
    site_management_locked_in_column = "locked_in")
  # tests
  expect_equal(evd_current, max(attr(evd_ss, "ev")))
})

test_that("current < optimal info, some pu selected", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(7),
      y = x,
      f1 = c(0, 1, 1, 0, NA, 1, NA),
      f2 = c(0, 1, 0, 1, NA, 0, NA),
      f3 = c(1, 0, 1, 0, NA, 1, NA),
      p1 = c(0.51, 0.99, 0.99, 0.05, 0.5, 0.99, 0.5),
      p2 = c(0.51, 0.99, 0.05, 0.99, 0.5, 0.05, 0.5),
      p3 = c(0.51, 0.05, 0.99, 0.99, 0.5, 0.99, 0.5),
      e1 = runif(7),
      e2 = runif(7),
      survey_cost = c(1, 1, 1, 1, 5, 100000, 8),
      management_cost = c(10, 10, 10, 10, 10, 10, 10),
      locked_in = FALSE),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = rep(TRUE, 3),
    survey_sensitivity = rep(0.95, 3),
    survey_specificity = rep(0.9, 3),
    model_sensitivity = rep(0.8, 3),
    model_specificity = rep(0.85, 3),
    target = c(3, 3, 3))
  xgb_parameters <-
    list(eta = c(0.1, 0.3, 0.5),
         lambda = c(0.01, 0.1, 0.5),
         objective = "binary:logistic")
  total_budget <- 57
  survey_budget <- 5
  # calculate expected values
  evd_current <- evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2", "f3"),
    site_probability_columns = c("p1", "p2", "p3"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    site_management_locked_in_column = "locked_in")
  evd_ss <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2", "f3"),
    site_probability_columns = c("p1", "p2", "p3"),
    site_env_vars_columns = c("e1", "e2"),
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    survey_budget = survey_budget,
    xgb_n_folds = rep(2, 3),
    xgb_tuning_parameters = xgb_parameters,
    site_management_locked_in_column = "locked_in")
  # tests
  expect_gt(max(attr(evd_ss, "ev")), evd_current)
})

test_that("locking out planning units lowers voi", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(7),
      y = x,
      f1 = c(0, 1, 1, 0, NA, 1, NA),
      f2 = c(0, 1, 0, 1, NA, 0, NA),
      f3 = c(1, 0, 1, 0, NA, 1, NA),
      p1 = c(0.51, 0.99, 0.99, 0.05, 0.5, 0.99, 0.5),
      p2 = c(0.51, 0.99, 0.05, 0.99, 0.5, 0.05, 0.5),
      p3 = c(0.51, 0.05, 0.99, 0.99, 0.5, 0.99, 0.5),
      e1 = runif(7),
      e2 = runif(7),
      survey_cost = c(1, 1, 1, 1, 5, 100000, 8),
      management_cost = c(10, 10, 10, 10, 10, 10, 10),
      locked_in = FALSE,
      locked_out = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = rep(TRUE, 3),
    survey_sensitivity = rep(0.95, 3),
    survey_specificity = rep(0.9, 3),
    model_sensitivity = rep(0.8, 3),
    model_specificity = rep(0.85, 3),
    target = c(3, 3, 3))
  xgb_parameters <-
    list(eta = c(0.1, 0.3, 0.5),
         lambda = c(0.01, 0.1, 0.5),
         objective = "binary:logistic")
  total_budget <- 1000
  survey_budget <- 5
  # calculate expected values
  ## evd current
  evd_ci1 <- evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2", "f3"),
    site_probability_columns = c("p1", "p2", "p3"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    site_management_locked_in_column = "locked_in")
  evd_ci2 <- evdci(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2", "f3"),
    site_probability_columns = c("p1", "p2", "p3"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    site_management_locked_in_column = "locked_in",
    site_management_locked_out_column = "locked_out")
  ## evdsi
  evd_opt1 <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2", "f3"),
    site_probability_columns = c("p1", "p2", "p3"),
    site_env_vars_columns = c("e1", "e2"),
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    survey_budget = survey_budget,
    xgb_n_folds = rep(2, 3),
    xgb_tuning_parameters = xgb_parameters,
    site_management_locked_in_column = "locked_in")
  evd_opt2 <- approx_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2", "f3"),
    site_probability_columns = c("p1", "p2", "p3"),
    site_env_vars_columns = c("e1", "e2"),
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = total_budget,
    survey_budget = survey_budget,
    xgb_n_folds = rep(2, 3),
    xgb_tuning_parameters = xgb_parameters,
    site_management_locked_in_column = "locked_in",
    site_management_locked_out_column = "locked_out")
  # tests
  expect_lt(evd_ci2, evd_ci1)
  expect_lt(max(attr(evd_opt2, "ev")), min(attr(evd_opt1, "ev")))
})

test_that("approx = exact, when all states used (evdsi)", {
  # data
  RandomFields::RFoptions(seed = 700)
  set.seed(500)
  n_f <- 2
  site_data <- simulate_site_data(n_sites = 15, n_features = n_f, 0.1)
  feature_data <- simulate_feature_data(n_features = n_f, 0.5)
  feature_data$target <- c(1, 1)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(is.na(site_data$f1))[1:2]] <- TRUE
  xgb_n_folds <- rep(5, n_f)
  xgb_parameters <-
    list(eta = c(0.1, 0.3, 0.5),
         lambda = c(0.01, 0.1, 0.5),
         objective = "binary:logistic")
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
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_tuning_parameters = xgb_parameters,
    xgb_n_folds = rep(5, n_f),
    seed = 1,
    n_approx_replicates = 5,
    n_approx_outcomes_per_replicate = n_states_per_rep,
    method_approx_outcomes = "uniform_without_replacement")
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
    feature_target_column = "target",
    total_budget = total_budget,
    xgb_tuning_parameters = xgb_parameters,
    xgb_n_folds = rep(5, n_f),
    seed = 1)
  expect_lte(max(abs(r1 - r2)), 1e-11)
})
