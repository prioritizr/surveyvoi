context("prior_probability_matrix")

test_that("expected result (survey performance > model performance)", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq(1, 6),
      y = x,
      f1 = c(1, 1, 1, 0, 0, 0),
      f2 = c(0, 1, 0, 0, 0, 0),
      f3 = c(0, 0, 1, 0, 0, 0),
      n1 = c(3, 3, 3, 0, 0, 0),
      n2 = c(1, 1, 1, 0, 0, 0),
      n3 = c(1, 1, 1, 0, 0, 0),
      p1 = c(0.99, 0.99, 0.99, 0.5, 0.1, 0.32),
      p2 = c(0.05, 0.99, 0.99, 0.2, 0.88, 0.67),
      p3 = c(0.21, 0.768, 0.98, 0.233, 0.56, 0.123)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    sensitivity = c(0.7, 0.96, 0.8),
    specificity = c(0.54, 0.92, 0.6),
    model_sensitivity = c(0.8, 0.7, 0.657),
    model_specificity = c(0.92, 0.9, 0.65))
  site_det_columns <- c("f1", "f2", "f3")
  site_prb_columns <- c("p1", "p2", "p3")
  site_n_columns <- c("n1", "n2", "n3")
  # calculations
  r <- prior_probability_matrix(
    site_data, feature_data, site_det_columns, site_n_columns, site_prb_columns,
    "sensitivity", "specificity", "model_sensitivity", "model_specificity")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), length(site_det_columns))
  expect_equal(ncol(r), nrow(site_data))
  expect_true(all(c(r) <= 1))
  expect_true(all(c(r) >= 0))
})


test_that("expected result (survey performance < model performance)", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq(1, 6),
      y = x,
      f1 = c(1, 1, 1, 0, 0, 0),
      f2 = c(0, 1, 0, 0, 0, 0),
      f3 = c(0, 0, 1, 0, 0, 0),
      n1 = c(1, 1, 1, 0, 0, 0),
      n2 = c(1, 5, 1, 0, 0, 0),
      n3 = c(1, 1, 1, 0, 0, 0),
      p1 = c(0.99, 0.99, 0.99, 0.5, 0.1, 0.32),
      p2 = c(0.05, 0.99, 0.99, 0.2, 0.88, 0.67),
      p3 = c(0.21, 0.768, 0.98, 0.233, 0.56, 0.123)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    sensitivity = c(0.7, 0.96, 0.8),
    specificity = c(0.54, 0.92, 0.6),
    model_sensitivity = c(0.99, 0.99, 0.99),
    model_specificity = c(0.99, 0.99, 0.99))
  site_det_columns <- c("f1", "f2", "f3")
  site_prb_columns <- c("p1", "p2", "p3")
  site_n_columns <- c("n1", "n2", "n3")
  # calculations
  pij <-
    t(as.matrix(sf::st_drop_geometry(site_data)[, c("p1", "p2", "p3"),
                                                drop = FALSE]))
  nij <-
    t(as.matrix(sf::st_drop_geometry(site_data)[, c("n1", "n2", "n3"),
                                                drop = FALSE]))
  r <- prior_probability_matrix(
    site_data, feature_data, site_det_columns, site_n_columns, site_prb_columns,
    "sensitivity", "specificity", "model_sensitivity", "model_specificity")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), length(site_det_columns))
  expect_equal(ncol(r), nrow(site_data))
  expect_lte(max(abs(r[pij >= 0.5 & nij == 1] - 0.99)), 1e-15)
  expect_lte(max(abs(r[pij < 0.5 & nij == 1] - 0.01)), 1e-15)
})

test_that("expected result (no survey data)", {
  # data
  set.seed(500)
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq(1, 6),
      y = x,
      zeros = 0,
      f1 = c(1, 1, 1, 0, 0, 0),
      f2 = c(0, 1, 0, 0, 0, 0),
      f3 = c(0, 0, 1, 0, 0, 0),
      n1 = c(1, 1, 1, 0, 0, 0),
      n2 = c(1, 1, 1, 0, 0, 0),
      n3 = c(1, 1, 1, 0, 0, 0),
      p1 = c(0.99, 0.99, 0.99, 0.5, 0.1, 0.32),
      p2 = c(0.05, 0.99, 0.99, 0.2, 0.88, 0.67),
      p3 = c(0.21, 0.768, 0.98, 0.233, 0.56, 0.123)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    sensitivity = c(0.7, 0.96, 0.8),
    specificity = c(0.54, 0.92, 0.6),
    model_sensitivity = c(0.99, 0.99, 0.99),
    model_specificity = c(0.99, 0.99, 0.99))
  site_det_columns <- c("f1", "f2", "f3")
  site_prb_columns <- c("p1", "p2", "p3")
  site_n_columns <- c("n1", "n2", "n3")
  # calculations
  ## note that surveys are always worse then models, so we should use
  ## models when specifying zero available surveys or even when survey data
  ## are present
  r1 <- prior_probability_matrix(
    site_data, feature_data,
    rep("zeros", length(site_det_columns)),
    rep("zeros", length(site_n_columns)),
    site_prb_columns,
    "sensitivity", "specificity", "model_sensitivity", "model_specificity")
  r2 <- prior_probability_matrix(
    site_data, feature_data, site_det_columns, site_n_columns, site_prb_columns,
    "sensitivity", "specificity", "model_sensitivity", "model_specificity")
  # tests
  expect_is(r1, "matrix")
  expect_is(r2, "matrix")
  rownames(r1) <- site_det_columns
  expect_equal(r1, r2)
})
