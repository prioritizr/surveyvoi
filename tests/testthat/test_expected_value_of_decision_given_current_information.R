context("expected_value_of_decision_given_current_information")

test_that("exact calculations", {
  # data
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      management_cost = c(100, 500, 200),
      locked_in = c(FALSE, FALSE, FALSE),
      f1 = c(1, 1, 0),
      f2 = c(0, 1, 0),
      p1 = c(0.99, 0.99, 0.99),
      p2 = c(0.05, 0.99, 0.99)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:2],
    survey = rep(TRUE, 2),
    sensitivity = c(0.5, 0.96),
    specificity = c(0.34, 0.92),
    model_sensitivity = c(0.8, 0.7),
    alpha = abs(rnorm(2)) + 1,
    gamma = runif(2))
  site_occupancy_columns <- c("f1", "f2")
  site_probability_columns <-  c("p1", "p2")
  prior_data <- prior_probability_matrix(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "sensitivity", "specificity", "model_sensitivity")
  # calculations
  set.seed(500)
  r1 <- rcpp_expected_value_of_decision_given_current_info(
    prior_data, site_data$management_cost, site_data$locked_in,
    feature_data$alpha, feature_data$gamma, 1000, 301, 0)
  r2 <- expected_value_of_decision_given_current_information(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "management_cost", "sensitivity", "specificity", "model_sensitivity",
    "alpha", "gamma", 301, "locked_in", NULL, 1000, NULL, 0, 500)
  # tests
  expect_equal(r1, r2)
})

test_that("approximate calculations", {
  # data
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      management_cost = c(100, 500, 200),
      locked_in = c(FALSE, FALSE, FALSE),
      f1 = c(1, 1, 0),
      f2 = c(0, 1, 0),
      p1 = c(0.99, 0.99, 0.99),
      p2 = c(0.05, 0.99, 0.99)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:2],
    survey = rep(TRUE, 2),
    sensitivity = c(0.5, 0.96),
    specificity = c(0.34, 0.92),
    model_sensitivity = c(0.8, 0.7),
    alpha = abs(rnorm(2)) + 1,
    gamma = runif(2))
  site_occupancy_columns <- c("f1", "f2")
  site_probability_columns <-  c("p1", "p2")
  prior_data <- prior_probability_matrix(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "sensitivity", "specificity", "model_sensitivity")
  # calculations
  set.seed(500)
  r1 <- rcpp_approx_expected_value_of_decision_given_current_info_n_states(
    prior_data, site_data$management_cost, site_data$locked_in,
    feature_data$alpha, feature_data$gamma, 1000, 301, 0, 20)
  r2 <- expected_value_of_decision_given_current_information(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "management_cost", "sensitivity", "specificity", "model_sensitivity",
    "alpha", "gamma", 301, "locked_in", NULL, 1000, 20, 0, 500)
  # tests
  expect_equal(r1, r2)
})
