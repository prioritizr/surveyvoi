context("posterior_probability_matrix")

test_that("correct result", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(6),
      y = x,
      solution = c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE),
      management_cost = c(100, 500, 200, 500, 12, 3),
      locked_in = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
      f1 = c(1, 1, 1, -1, -1, -1),
      f2 = c(0, 1, 0, -1, -1, -1),
      f3 = c(0, 0, 0, -1, -1, -1),
      p1 = c(0.99, 0.99, 0.99, 0.5, 0.1, 0.32),
      p2 = c(0.05, 0.99, 0.99, 0.2, 0.88, 0.67),
      p3 = c(0.21, 0.768, 0.98, 0.233, 0.56, 0.123),
      o1 = c(1, 1, 1, 1, 0.7, 0),
      o2 = c(0, 1, 0, -1, -1, -1),
      o3 = c(0, 1, 0, 0, 0.2, 0)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = c(TRUE, FALSE, TRUE),
    sensitivity = c(0.7, 0.96, 0.8),
    specificity = c(0.54, 0.92, 0.6),
    model_sensitivity = c(0.8, 0.7, 0.657),
    model_specificity = c(0.92, 0.9, 0.65))
  site_occ_columns <- c("f1", "f2", "f3")
  site_prb_columns <- c("p1", "p2", "p3")
  site_out_columns <- c("o1", "o2", "o3") # combination of outcome/new models
  # extract data
  rij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_occ_columns])))
  pij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_prb_columns])))
  oij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_out_columns])))
  # calculations
  r1 <- rcpp_posterior_probability_matrix(
    rij, pij, oij,
    site_data$solution,
    feature_data$survey,
    feature_data$sensitivity, feature_data$specificity,
    feature_data$model_sensitivity, feature_data$model_specificity)
  r2 <- r_posterior_probability_matrix(
    rij, pij, oij,
    site_data$solution,
    feature_data$survey,
    feature_data$sensitivity, feature_data$specificity,
    feature_data$model_sensitivity, feature_data$model_specificity)
  dimnames(r1) <- dimnames(pij)
  dimnames(r2) <- dimnames(pij)
  # tests
  ## R and Rcpp give same answers
  expect_equal(r1, r2)
  ## feature 1 posterior
  expect_equal(r1[1, 1:3], pij[1, 1:3])
  expect_gt(r1[1, 4], pij[1, 4])
  expect_gt(r1[1, 5], pij[1, 5])
  expect_lt(r1[1, 6], pij[1, 6])
  ## feature 2 posterior
  expect_equal(r1[2, ], pij[2, ])
  ## feature 3 posterior
  expect_equal(r1[3, 1:3], pij[3, 1:3])
  expect_lt(r1[3, 4], pij[3, 4])
  expect_lt(r1[3, 5], pij[3, 5])
  expect_lt(r1[3, 6], pij[3, 6])
})
