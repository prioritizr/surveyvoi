context("approx_evdci")

test_that("expected result", {
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
    model_specificity = c(0.92, 0.9),
    alpha = abs(rnorm(2)) + 1,
    gamma = runif(2))
  site_occupancy_columns <- c("f1", "f2")
  site_probability_columns <-  c("p1", "p2")
  prior_data <- prior_probability_matrix(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "sensitivity", "specificity", "model_sensitivity", "model_specificity")
  # calculations
  set.seed(500)
  r1 <- rcpp_approx_expected_value_of_decision_given_current_info_n_states(
    prior_data, site_data$management_cost, site_data$locked_in,
    feature_data$alpha, feature_data$gamma, 1000, 301, 0, 10, 20)
  r2 <- approx_evdci(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "management_cost", "sensitivity", "specificity", "model_sensitivity",
    "model_specificity", "alpha", "gamma", 301, "locked_in", NULL, 1000, 0, 10,
    20, 500)
  # tests
  expect_equal(setNames(r1, c("mean", "se")), r2)
})
