context("rcpp_expected_value_of_action")

test_that("correct result (target = 1)", {
  # data
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      solution = c(TRUE, FALSE, TRUE),
      f1 = c(1, 1, 1),
      f2 = c(0, 1, 0),
      n1 = c(1, 1, 1),
      n2 = c(1, 1, 1),
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
    target = c(1, 1))
  site_det_columns <- c("f1", "f2")
  site_n_columns <- c("n1", "n2")
  site_prb_columns <-  c("p1", "p2")
  prior_data <- prior_probability_matrix(
    site_data, feature_data,
    site_det_columns, site_n_columns, site_prb_columns,
    "sensitivity", "specificity",
    "model_sensitivity", "model_specificity")
  # calculations
  r1 <- rcpp_expected_value_of_action(
    site_data$solution, prior_data, feature_data$target)
  r2 <- r_expected_value_of_action(
    site_data$solution, prior_data, feature_data$target)
  # tests
  expect_equal(r1, r2)
})

test_that("correct result (target = c(2, 2))", {
  # data
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(3),
      y = x,
      solution = c(TRUE, FALSE, TRUE),
      f1 = c(1, 1, 1),
      f2 = c(0, 1, 0),
      n1 = c(1, 1, 1),
      n2 = c(1, 1, 1),
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
    target = c(2, 2))
  site_det_columns <- c("f1", "f2")
  site_n_columns <- c("n1", "n2")
  site_prb_columns <-  c("p1", "p2")
  prior_data <- prior_probability_matrix(
    site_data, feature_data,
    site_det_columns, site_n_columns, site_prb_columns,
    "sensitivity", "specificity",
    "model_sensitivity", "model_specificity")
  # calculations
  r1 <- rcpp_expected_value_of_action(
    site_data$solution, prior_data, feature_data$target)
  r2 <- r_expected_value_of_action(
    site_data$solution, prior_data, feature_data$target)
  # tests
  expect_equal(r1, r2)
})

test_that("correct result (target = c(10, 10))", {
  # data
  set.seed(123)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(20),
      y = x,
      solution = sample(c(TRUE, FALSE), size = 20, replace = TRUE),
      f1 = sample(c(1, 0), size = 20, replace = TRUE),
      f2 = sample(c(1, 0), size = 20, replace = TRUE),
      n1 = sample(c(1, 0), size = 20, replace = TRUE),
      n2 =sample(c(1, 0), size = 20, replace = TRUE),
      p1 = runif(20),
      p2 = runif(20)
    ),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:2],
    survey = rep(TRUE, 2),
    sensitivity = c(0.5, 0.96),
    specificity = c(0.34, 0.92),
    model_sensitivity = c(0.8, 0.7),
    model_specificity = c(0.92, 0.9),
    target = c(10, 10))
  site_det_columns <- c("f1", "f2")
  site_n_columns <- c("n1", "n2")
  site_prb_columns <-  c("p1", "p2")
  prior_data <- prior_probability_matrix(
    site_data, feature_data,
    site_det_columns, site_n_columns, site_prb_columns,
    "sensitivity", "specificity",
    "model_sensitivity", "model_specificity")
  # calculations
  r1 <- rcpp_expected_value_of_action(
    site_data$solution, prior_data, feature_data$target)
  r2 <- r_expected_value_of_action(
    site_data$solution, prior_data, feature_data$target)
  # tests
  expect_equal(r1, r2)
})
