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
      n1 = c(1, 1, 3, 0, 0, 0), # number of surveys for f1
      n2 = c(1, 1, 1, 0, 0, 0), # number of surveys for f2
      n3 = c(2, 1, 1, 0, 0, 0), # number of surveys for f3
      p1 = c(0.99, 0.99, 0.99, 0.5, 0.1, 0.32), # prior for f1
      p2 = c(0.05, 0.99, 0.99, 0.2, 0.88, 0.67), # prior for f2
      p3 = c(0.21, 0.768, 0.98, 0.233, 0.56, 0.123), # prior for f3
      o1 = c(-1, -1, -1,  1, 0.9, 0.3), # outcomes and model probs for f1
      o2 = c(-1, -1, -1, -1, -1, -1), # outcomes and model probs for f2
      o3 = c(-1, -1, -1,  0, 0.8, 0.4)), # outcomes and model probs for f3
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = c(TRUE, FALSE, TRUE),
    sensitivity = c(0.9, 0.96, 0.9),
    specificity = c(0.95, 0.92, 0.95),
    model_sensitivity = c(0.8, 0.7, 0.657),
    model_specificity = c(0.92, 0.9, 0.65))
  site_n_columns <- c("n1", "n2", "n3")
  site_prb_columns <- c("p1", "p2", "p3")
  site_out_columns <- c("o1", "o2", "o3") # combination of outcome/new models
  # extract data
  nij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_n_columns])))
  pij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_prb_columns])))
  oij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_out_columns])))
  # calculations
  r1 <- rcpp_posterior_probability_matrix(
    nij, pij, oij,
    site_data$solution,
    feature_data$survey,
    feature_data$sensitivity, feature_data$specificity,
    feature_data$model_sensitivity, feature_data$model_specificity)
  r2 <- r_posterior_probability_matrix(
    nij, pij, oij,
    site_data$solution,
    feature_data$survey,
    feature_data$sensitivity, feature_data$specificity,
    feature_data$model_sensitivity, feature_data$model_specificity)
  dimnames(r1) <- dimnames(pij)
  dimnames(r2) <- dimnames(pij)
  # tests
  ## expect prior probs in places with existing data
  expect_equivalent(r2[1, 1:3], site_data$p1[1:3])
  expect_equivalent(r2[2, ], site_data$p2)
  expect_equivalent(r2[3, 1:3], site_data$p3[1:3])
  ## expect highest probs in places with surveys that have detection
  expect_gt(r2[1, 4], max(r2[1, c(5, 6)]))
  ## expect lowest probs in places with surveys that have detection
  expect_lt(r2[3, 4], max(r2[3, c(5, 6)]))
  ## expect high probs in places with modelled presence
  expect_gt(r2[1, 5], max(r2[1, 6]))
  ## expect low probs in places with modelled presence
  expect_lt(r2[3, 6], max(r2[3, 5]))
  ## equal results for r and rcpp implementations
  expect_equal(r1, r2)
})
