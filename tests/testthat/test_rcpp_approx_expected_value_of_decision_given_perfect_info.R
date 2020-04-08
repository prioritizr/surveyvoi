context("rcpp_approx_expected_value_of_decision_given_perfect_info")

test_that("correct result (fixed states)", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      management_cost = c(100, 500, 200),
      locked_in = c(FALSE, FALSE, FALSE),
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
    model_specificity = c(0.92, 0.9),
    preweight = runif(2, 100, 200),
    postweight = runif(2, 5, 20),
    target = c(1, 1))
  site_data <- site_data[c(1, 2), ]
  site_occupancy_columns <- c("f1", "f2")
  site_probability_columns <-  c("p1", "p2")
  prior_data <- prior_probability_matrix(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "sensitivity", "specificity", "model_sensitivity", "model_specificity")
  states <- c(0, 1, 5, 10, 20, 23)
  # calculations
  r1 <- rcpp_approx_expected_value_of_decision_given_perfect_info_fixed_states(
    prior_data, site_data$management_cost, site_data$locked_in,
    feature_data$preweight, feature_data$postweight, feature_data$target,
    1000, 301, 0, states)
  r2 <- r_approx_expected_value_of_decision_given_perfect_info_fixed_states(
    prior_data, site_data$management_cost, site_data$locked_in,
    feature_data$preweight, feature_data$postweight, feature_data$target,
    1000, 301, 0, states)
  # tests
  expect_equal(r1, r2)
})

test_that("correct result (n states)", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      management_cost = c(100, 500, 200),
      locked_in = c(FALSE, FALSE, FALSE),
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
    model_specificity = c(0.92, 0.9),
    preweight = runif(2, 100, 200),
    postweight = runif(2, 5, 20),
    target = c(1, 1))
  site_occupancy_columns <- c("f1", "f2")
  site_probability_columns <-  c("p1", "p2")
  prior_data <- prior_probability_matrix(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "sensitivity", "specificity", "model_sensitivity", "model_specificity")
  states <- 10
  reps <- 5
  # calculations
  set.seed(100)
  r1 <- rcpp_approx_expected_value_of_decision_given_perfect_info_n_states(
    prior_data, site_data$management_cost, site_data$locked_in,
    feature_data$preweight, feature_data$postweight, feature_data$target,
    1000, 301, 0, reps, states,
    method_approx_states = "weighted_without_replacement")
  set.seed(100)
  r2 <- r_approx_expected_value_of_decision_given_perfect_info_n_states(
    prior_data, site_data$management_cost, site_data$locked_in,
    feature_data$preweight, feature_data$postweight, feature_data$target,
    1000, 301, 0, reps, states)
  # tests
  expect_lte(max(abs(c(r1 - r2))), 1e-15)
})

test_that("expected error", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      management_cost = c(100, 500, 200),
      locked_in = c(FALSE, FALSE, FALSE),
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
    model_specificity = c(0.92, 0.9),
    preweight = runif(2, 100, 200),
    postweight = runif(2, 5, 20),
    target = c(1, 1))
  site_occupancy_columns <- c("f1", "f2")
  site_probability_columns <-  c("p1", "p2")
  prior_data <- prior_probability_matrix(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "sensitivity", "specificity", "model_sensitivity", "model_specificity")
  # tests
  expect_error({
    rcpp_approx_expected_value_of_decision_given_perfect_info_fixed_states(
      prior_data, site_data$management_cost, site_data$locked_in,
      feature_data$preweight, feature_data$postweight, feature_data$target,
      1000, 301, 0, 0)
  })
})
