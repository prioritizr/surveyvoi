context("rcpp_approximate_expected_value_of_management_action")

test_that("expected result", {
  # data
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      solution = c(TRUE, FALSE, TRUE),
      f1 = c(1, 1, 1),
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
  states <- c(0, 1, 4, 8, 10)
  # calculations
  r1 <- rcpp_approximate_expected_value_of_management_action(
    site_data$solution, prior_data, feature_data$alpha, feature_data$gamma,
    states)
  r2 <- r_approximate_expected_value_of_management_action(
    site_data$solution, prior_data, feature_data$alpha, feature_data$gamma,
    states)
  # tests
  expect_equal(r1, r2)
})

test_that("consistent result", {
  # data
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      solution = c(TRUE, FALSE, TRUE),
      f1 = c(1, 1, 1),
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
  states <- c(0, 1, 4, 8, 10)
  # calculations
  r1 <- sapply(seq_len(50), function(x) {
    rcpp_approximate_expected_value_of_management_action(
      site_data$solution, prior_data, feature_data$alpha, feature_data$gamma,
      states)
  })
  r2 <- sapply(seq_len(50), function(x) {
    r_approximate_expected_value_of_management_action(
      site_data$solution, prior_data, feature_data$alpha, feature_data$gamma,
      states)
  })
  # tests
  expect_length(unique(c(r1, r2)), 1)
})
