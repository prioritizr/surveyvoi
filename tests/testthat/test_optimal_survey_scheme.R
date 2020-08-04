context("optimal_survey_scheme")

test_that("expected results", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(5),
      y = x,
      f1 = c(0, 1, 1, NA, 1),
      f2 = c(0, 1, 0, NA, 0),
      f3 = c(0, 0, 1, NA, 1),
      p1 = c(0.05, 0.99, 0.99, 0.99, 0.99),
      p2 = c(0.05, 0.99, 0.05, 0.05, 0.99),
      p3 = c(0.05, 0.05, 0.99, 0.05, 0.99),
      e1 = runif(5),
      e2 = runif(5),
      survey_cost = c(1, 1, 1, 5, 100000),
      management_cost = c(10, 10, 10, 10, 10),
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
  xgb_parameters <- list(list(nrounds = 3, eta = 0.3, scale_pos_weight = 2,
                              objective = "binary:logistic"))[rep(1, 3)]
  # generate prioritisation
  r <- optimal_survey_scheme(
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
    total_budget = 49,
    survey_budget = 10,
    xgb_parameters = xgb_parameters,
    site_management_locked_in_column = "locked_in")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(site_data))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, FALSE, FALSE, TRUE, FALSE))
})

test_that("consistent results", {
  # seeds
  set.seed(505)
  RandomFields::RFoptions(seed = 505)
  # data
  site_data <- simulate_site_data(
    n_sites = 10, n_features = 1, proportion_of_sites_missing_data = 0.5,
    n_env_vars = 2, output_probabilities = FALSE)
  feature_data <- simulate_feature_data(
    n_features = 1, proportion_of_survey_features = 1)
  total_budget <- 500.128863597055
  survey_budget <- 26.4498218037885
  xgb_params <- list(
    max_depth = 5, eta = 0.3, lambda = 1, subsample = 1, colsample_bytree = 1,
    objective = "binary:logistic")
  xgb_model <- fit_occupancy_models(
    site_data, "f1", c("e1", "e2"), parameters = xgb_params,
    n_random_search_iterations = 1)
  site_data$p1 <- xgb_model$predictions$f1
  # run calculations
  r <- lapply(seq_len(5), function(i) {
    optimal_survey_scheme(
        site_data = site_data,
        feature_data = feature_data,
        site_occupancy_columns = "f1",
        site_probability_columns = "p1",
        site_env_vars_columns = c("e1", "e2"),
        site_survey_cost_column = "survey_cost",
        site_management_cost_column = "management_cost",
        feature_survey_column = "survey",
        feature_survey_sensitivity_column = "survey_sensitivity",
        feature_survey_specificity_column = "survey_specificity",
        feature_model_sensitivity_column = "model_sensitivity",
        feature_model_specificity_column = "model_specificity",
        feature_target_column = "target",
        survey_budget = survey_budget,
        total_budget = total_budget,
        xgb_parameters = xgb_model$parameters)
  })
  # verify that all repeat calculations are identical
  for (i in seq_along(r))
    expect_identical(r[[1]], r[[i]])
})

test_that("consistent results (multiple threads)", {
  # seeds
  set.seed(505)
  RandomFields::RFoptions(seed = 505)
  # data
  site_data <- simulate_site_data(
    n_sites = 10, n_features = 1, proportion_of_sites_missing_data = 0.5,
    n_env_vars = 2, output_probabilities = FALSE)
  feature_data <- simulate_feature_data(
    n_features = 1, proportion_of_survey_features = 1)
  total_budget <- 500.128863597055
  survey_budget <- 26.4498218037885
  xgb_params <- list(
    max_depth = 5, eta = 0.3, lambda = 1, subsample = 1, colsample_bytree = 1,
    objective = "binary:logistic")
  xgb_model <- fit_occupancy_models(
    site_data, "f1", c("e1", "e2"), parameters = xgb_params,
    n_random_search_iterations = 1)
  site_data$p1 <- xgb_model$predictions$f1
  # run calculations
  r <- suppressWarnings({
    lapply(seq_len(5), function(i) {
      optimal_survey_scheme(
        site_data = site_data,
        feature_data = feature_data,
        site_occupancy_columns = "f1",
        site_probability_columns = "p1",
        site_env_vars_columns = c("e1", "e2"),
        site_survey_cost_column = "survey_cost",
        site_management_cost_column = "management_cost",
        feature_survey_column = "survey",
        feature_survey_sensitivity_column = "survey_sensitivity",
        feature_survey_specificity_column = "survey_specificity",
        feature_model_sensitivity_column = "model_sensitivity",
        feature_model_specificity_column = "model_specificity",
        feature_target_column = "target",
        survey_budget = survey_budget,
        total_budget = total_budget,
        xgb_parameters = xgb_model$parameters)
    })
  })
  # verify that all repeat calculations are identical
  for (i in seq_along(r))
    expect_identical(r[[1]], r[[i]])
})

test_that("expected results (sparse)", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(5),
      y = x,
      f1 = c(0, 1, 1, 0, NA),
      f2 = c(0, 1, 0, NA, 1),
      f3 = c(0, 0, 1, NA, NA),
      p1 = c(0.99, 0.99, 0.99, 0.99, 0.99),
      p2 = c(0.05, 0.99, 0.05, 0.05, 0.99),
      p3 = c(0.05, 0.05, 0.05, 0.05, 0.99),
      e1 = rnorm(5),
      e2 = rnorm(5),
      survey_cost = c(1, 1, 1, 5, 100000),
      survey_locked_out = c(TRUE, TRUE, FALSE, FALSE, FALSE),
      management_cost = c(10, 10, 10, 10, 10),
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
  xgb_parameters <- list(list(nrounds = 3, eta = 0.3, scale_pos_weight = 2,
                              objective = "binary:logistic"))[rep(1, 3)]
  # generate prioritisation
  r <- optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = c("f1", "f2", "f3"),
    site_probability_columns = c("p1", "p2", "p3"),
    site_env_vars_columns = c("e1", "e2"),
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    site_survey_locked_out_column = "survey_locked_out",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = 49,
    survey_budget = 10,
    xgb_parameters = xgb_parameters,
    site_management_locked_in_column = "locked_in")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(site_data))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, FALSE, FALSE, TRUE, FALSE))
})
