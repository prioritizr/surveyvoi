context("approx_near_optimal_survey_scheme")

test_that("single species", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(6),
      y = x,
      f1 = c(0, 1, 0, 1, 0, 1),
      n1 = c(1, 1, 1, 1, 0, 1),
      p1 = c(0.51, 0.99, 0.01, 0.99, 0.5, 0.99),
      survey_cost = c(1, 1, 1, 1, 5, 100000),
      management_cost = c(10, 10, 10, 10, 10, 10),
      locked_in = FALSE),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:1],
    survey = rep(TRUE, 1),
    survey_sensitivity = rep(0.99, 1),
    survey_specificity = rep(0.99, 1),
    model_sensitivity = rep(0.5, 1),
    model_specificity = rep(0.5, 1),
    target = 3)
  pm <- matrix(NA, ncol = nrow(site_data), nrow = nrow(feature_data))
  pm[1, ] <- site_data$p1
  # generate prioritisation
  r <- approx_near_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = c("f1"),
    site_n_surveys_columns = c("n1"),
    site_probability_columns = c("p1"),
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = 46,
    survey_budget = 10,
    site_management_locked_in_column = "locked_in",
    prior_matrix = pm)
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(site_data))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
})

test_that("multiple species", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(6),
      y = x,
      f1 = c(0, 1, 1, 0, 0, 1),
      f2 = c(0, 1, 0, 1, 0, 0),
      f3 = c(1, 0, 1, 0, 0, 1),
      n1 = c(1, 1, 1, 1, 0, 1),
      n2 = c(1, 1, 1, 1, 0, 1),
      n3 = c(1, 1, 1, 1, 0, 1),
      p1 = c(0.51, 0.99, 0.99, 0.05, 0.5, 0.99),
      p2 = c(0.51, 0.99, 0.05, 0.99, 0.5, 0.05),
      p3 = c(0.51, 0.05, 0.99, 0.99, 0.5, 0.99),
      survey_cost = c(1, 1, 1, 1, 5, 100000),
      management_cost = c(10, 10, 10, 10, 10, 10),
      locked_in = FALSE),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = rep(TRUE, 3),
    survey_sensitivity = rep(0.99, 3),
    survey_specificity = rep(0.99, 3),
    model_sensitivity = rep(0.8, 3),
    model_specificity = rep(0.85, 3),
    target = c(3, 3, 3))
  pij <-
    t(as.matrix(sf::st_drop_geometry(site_data)[, c("p1", "p2", "p3"),
                                                drop = FALSE]))
  # generate prioritisation
  r <- approx_near_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = c("f1", "f2", "f3"),
    site_n_surveys_columns = c("n1", "n2", "n3"),
    site_probability_columns = c("p1", "p2", "p3"),
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = 56,
    survey_budget = 10,
    site_management_locked_in_column = "locked_in",
    prior_matrix = pij)
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(site_data))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
})

test_that("multiple species (sparse)", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(6),
      y = x,
      f1 = c(0, 1, 1, 0, 0, 0),
      f2 = c(0, 1, 0, 1, 0, 0),
      f3 = c(1, 0, 1, 0, 0, 1),
      n1 = c(1, 1, 1, 1, 0, 0),
      n2 = c(1, 1, 1, 1, 0, 1),
      n3 = c(1, 1, 1, 1, 0, 1),
      p1 = c(0.51, 0.99, 0.99, 0.05, 0.5, 1e-10),
      p2 = c(0.51, 0.99, 0.05, 0.99, 0.5, 0.05),
      p3 = c(0.51, 0.05, 0.99, 0.99, 0.5, 0.99),
      survey_cost = c(1, 1, 1, 1, 5, 100000),
      management_cost = c(10, 10, 10, 10, 10, 10),
      locked_in = FALSE),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = rep(TRUE, 3),
    survey_sensitivity = rep(0.99, 3),
    survey_specificity = rep(0.99, 3),
    model_sensitivity = rep(0.8, 3),
    model_specificity = rep(0.85, 3),
    target = c(3, 3, 3))
  pij <-
    t(as.matrix(sf::st_drop_geometry(site_data)[, c("p1", "p2", "p3"),
                                                drop = FALSE]))
  # generate prioritisation
  r <- approx_near_optimal_survey_scheme(
    site_data = site_data,
    feature_data = feature_data,
    site_detection_columns = c("f1", "f2", "f3"),
    site_n_surveys_columns = c("n1", "n2", "n3"),
    site_probability_columns = c("p1", "p2", "p3"),
    site_management_cost_column = "management_cost",
    site_survey_cost_column = "survey_cost",
    feature_survey_column = "survey",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_target_column = "target",
    total_budget = 56,
    survey_budget = 10,
    site_management_locked_in_column = "locked_in",
    prior_matrix = pij)
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(site_data))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
})

test_that("consistent results", {
  skip_if_not_installed("RandomFields")
  skip_on_os("Windows")
  # seeds
  set.seed(505)
  RandomFields::RFoptions(seed = 505)
  # data
  site_data <- simulate_site_data(
    n_sites = 30, n_features = 1, proportion_of_sites_missing_data = 0.1,
    n_env_vars = 2, output_probabilities = FALSE)
  feature_data <- simulate_feature_data(
    n_features = 1, proportion_of_survey_features = 1)
  total_budget <- 500.128863597055
  survey_budget <- 26.4498218037885
  xgb_parameters <-
    list(eta = c(0.1, 0.3),
         lambda = c(0.01),
         objective = "binary:logistic")
  xgb_model <- fit_xgb_occupancy_models(
    site_data, feature_data, "f1", "n1", c("e1", "e2"),
    "survey_sensitivity", "survey_specificity",
    xgb_tuning_parameters = xgb_parameters,
    xgb_early_stopping_rounds = 5, xgb_n_rounds = 10, n_folds = 2)
  site_data$p1 <- xgb_model$predictions$f1
  # run calculations
  r <- lapply(seq_len(5), function(i) {
    approx_near_optimal_survey_scheme(
        site_data = site_data,
        feature_data = feature_data,
        site_detection_columns = "f1",
        site_n_surveys_columns = "n1",
        site_probability_columns = "p1",
        site_survey_cost_column = "survey_cost",
        site_management_cost_column = "management_cost",
        feature_survey_column = "survey",
        feature_survey_sensitivity_column = "survey_sensitivity",
        feature_survey_specificity_column = "survey_specificity",
        feature_model_sensitivity_column = "model_sensitivity",
        feature_model_specificity_column = "model_specificity",
        feature_target_column = "target",
        survey_budget = survey_budget,
        total_budget = total_budget)
  })
  # verify that all repeat calculations are identical
  for (i in seq_along(r))
    expect_identical(r[[1]], r[[i]])
})

test_that("consistent results (multiple threads)", {
  # skip if using PSOCK cluster and package not installed
  skip_if_not_installed("RandomFields")
  skip_on_os("Windows")
  skip_if(!requireNamespace("surveyvoi") &&
          !identical(.Platform$OS.type, "unix"))
  # seeds
  set.seed(505)
  RandomFields::RFoptions(seed = 505)
  # data
  site_data <- simulate_site_data(
    n_sites = 30, n_features = 1, proportion_of_sites_missing_data = 0.1,
    n_env_vars = 2, output_probabilities = FALSE)
  feature_data <- simulate_feature_data(
    n_features = 1, proportion_of_survey_features = 1)
  total_budget <- 500.128863597055
  survey_budget <- 26.4498218037885
  xgb_parameters <-
    list(eta = c(0.1, 0.3),
         lambda = c(0.01),
         objective = "binary:logistic")
  xgb_model <- fit_xgb_occupancy_models(
    site_data, feature_data, "f1", "n1", c("e1", "e2"),
    "survey_sensitivity", "survey_specificity",
    xgb_tuning_parameters = xgb_parameters,
    xgb_early_stopping_rounds = 5, xgb_n_rounds = 10, n_folds = 2)
  site_data$p1 <- xgb_model$predictions$f1
  # run calculations
  suppressWarnings({
    r <- lapply(seq_len(5), function(i) {
      approx_near_optimal_survey_scheme(
          site_data = site_data,
          feature_data = feature_data,
          site_detection_columns = "f1",
          site_n_surveys_columns = "n1",
          site_probability_columns = "p1",
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
          n_threads = 2)
    })
  })
  # verify that all repeat calculations are identical
  for (i in seq_along(r))
    expect_identical(r[[1]], r[[i]])
})
