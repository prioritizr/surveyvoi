context("rcpp_approx_expected_value_of_decision_given_survey_scheme")

test_that("single species", {
  skip_if_not_installed("RandomFields")
  skip_on_os("windows")
  # data
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_f <- 1
  n_sites <- 20
  n_vars <- 3
  n_reps <- 2
  n_outcomes_per_rep <- 2
  ## simulate data
  site_data <- simulate_site_data(n_sites, n_f, 0.5, n_vars)
  feature_data <- simulate_feature_data(n_f)
  feature_data$target <- c(6)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(site_data$n1 < 0.5)[1:4]] <- TRUE
  ## create matrices for data
  pij <- prior_probability_matrix(
    site_data, feature_data,
    paste0("f", seq_len(n_f)),
    paste0("n", seq_len(n_f)), paste0("p", seq_len(n_f)),
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  # calculations
  r1 <- rcpp_approx_expected_value_of_decision_given_survey_scheme(
    pij = pij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    seed = 500)
  r2 <- r_approx_expected_value_of_decision_given_survey_scheme(
    pij = pij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
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
  skip_if_not_installed("RandomFields")
  skip_on_os("windows")
  # data
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_f <- 3
  n_sites <- 20
  n_vars <- 3
  n_reps <- 1
  n_outcomes_per_rep <- 2
  ## simulate data
  site_data <- simulate_site_data(n_sites, n_f, 0.5, n_vars)
  feature_data <- simulate_feature_data(n_f)
  feature_data$target <- c(6)
  feature_data$survey <- c(TRUE, FALSE, TRUE)
  total_budget <- sum(site_data$management_cost * 0.8)
  site_data$survey <- FALSE
  site_data$survey[which(site_data$n1 < 0.5)[1:4]] <- TRUE
  ## create matrices for data
  pij <- prior_probability_matrix(
    site_data, feature_data,
    paste0("f", seq_len(n_f)),
    paste0("n", seq_len(n_f)), paste0("p", seq_len(n_f)),
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  # calculations
  r1 <- rcpp_approx_expected_value_of_decision_given_survey_scheme(
    pij = pij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
    obj_fun_target = feature_data$target,
    total_budget = total_budget,
    n_approx_replicates = n_reps,
    n_approx_outcomes_per_replicate = n_outcomes_per_rep,
    seed = 500)
  r2 <- r_approx_expected_value_of_decision_given_survey_scheme(
    pij = pij,
    survey_features = feature_data$survey,
    survey_sensitivity = feature_data$survey_sensitivity,
    survey_specificity = feature_data$survey_specificity,
    pu_survey_solution = site_data$survey,
    pu_survey_costs  = site_data$survey_cost,
    pu_purchase_costs = site_data$management_cost,
    pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
    pu_purchase_locked_out = rep(FALSE, nrow(site_data)),
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
