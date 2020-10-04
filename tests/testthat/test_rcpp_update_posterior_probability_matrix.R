context("rcpp_update_posterior_probability_matrix")

test_that("correct result", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(6),
      y = x,
      survey = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
      p1 = c(0.99, 0.99, 0.99, 0.51, 0.1, 0.32), # prior for f1
      p2 = c(0.05, 0.99, 0.99, 0.2, 0.88, 0.67), # prior for f2
      p3 = c(0.21, 0.768, 0.98, 0.233, 0.56, 0.123), # prior for f3
      o1 = c(-1, -1, -1,  1, 1, 0), # survey outcomes for f1
      o2 = c(-1, -1, -1, -1, -1, -1), # survey outcomes for f2
      o3 = c(-1, -1, -1,  0, 1, 0)), # survey outcomes for f3
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = c(TRUE, FALSE, TRUE),
    sensitivity = c(0.9, 0.96, 0.52),
    specificity = c(0.95, 0.92, 0.51))
  site_pij_columns <- c("p1", "p2", "p3")
  site_oij_columns <- c("o1", "o2", "o3")
  # extract data
  pij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_pij_columns])))
  oij <- t(as.matrix(sf::st_drop_geometry(site_data[, site_oij_columns])))
  # calculations
  r1 <- rcpp_update_posterior_probability_matrix(
    pij, oij,
    feature_data$survey,
    feature_data$sensitivity, feature_data$specificity,
    site_data$survey)
  r2 <- r_update_posterior_probability_matrix(
    pij, oij,
    feature_data$survey,
    feature_data$sensitivity, feature_data$specificity,
    site_data$survey)
  dimnames(r1) <- dimnames(pij)
  dimnames(r2) <- dimnames(pij)
  # tests
  ## equal results for r and rcpp implementations
  expect_equal(r1, r2)
  ## expect same posterior probs in places with no surveys
  expect_equivalent(r2[1, 1:3], pij[1, 1:3])
  expect_equivalent(r2[2, ], pij[2, ])
  expect_equivalent(r2[3, 1:3], pij[3, 1:3])
  ## expect different posterior probs in places with surveys
  expect_true(all(r2[1, 4:6] != pij[1, 4:6]))
  expect_true(all(r2[3, 4:6] != pij[3, 4:6]))
  ## expect highest probs in places with surveys that have detection
  expect_gt(r2[1, 4], max(r2[1, c(5, 6)]))
  ## expect lowest probs in places with surveys that have detection
  expect_lt(r2[3, 4], max(r2[3, c(5, 6)]))
})
