load_all()
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(5),
      y = x,
      f1 = c(0, 1, 1, NA, NA),
      f2 = c(0, 1, 0, NA, NA),
      f3 = c(0, 0, 1, NA, NA),
      p1 = c(0.99, 0.99, 0.99, 0.99, 0.99),
      p2 = c(0.05, 0.99, 0.05, 0.05, 0.99),
      p3 = c(0.05, 0.05, 0.05, 0.05, 0.99),
      e1 = rnorm(5),
      e2 = rnorm(5),
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
    preweight = runif(3, 100, 200),
    postweight = runif(3, 5, 20),
    target = c(1, 1, 3))
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = 49,
    survey_budget = 10,
    xgb_parameters = xgb_parameters,
    site_management_locked_in_column = "locked_in",
    optimality_gap = 0)
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(site_data))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, FALSE, FALSE, TRUE, FALSE))
load_all()
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(5),
      y = x,
      f1 = c(0, 1, 1, NA, NA),
      f2 = c(0, 1, 0, NA, NA),
      f3 = c(0, 0, 1, NA, NA),
      p1 = c(0.99, 0.99, 0.99, 0.99, 0.99),
      p2 = c(0.05, 0.99, 0.05, 0.05, 0.99),
      p3 = c(0.05, 0.05, 0.05, 0.05, 0.99),
      e1 = rnorm(5),
      e2 = rnorm(5),
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
    preweight = runif(3, 100, 200),
    postweight = runif(3, 5, 20),
    target = c(1, 1, 3))
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
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = 49,
    survey_budget = 10,
    xgb_parameters = xgb_parameters,
    site_management_locked_in_column = "locked_in",
    optimality_gap = 0)
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(site_data))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, FALSE, FALSE, TRUE, FALSE))
