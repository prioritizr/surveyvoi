context("prior_probability_matrix")

test_that("correct result", {
  # data
  site_data <- sf::st_as_sf(
    tibble::tibble(
      x = seq_len(5),
      y = x,
      f1 = c(1, 1, 1, NA, NA),
      f2 = c(0, 1, 0, NA, NA),
      f3 = c(0, 0, 0, NA, NA),
      p1 = c(0.99, 0.99, 0.99, 0.99, 0.99),
      p2 = c(0.05, 0.99, 0.99, 0.05, 0.99),
      p3 = c(0.05, 0.05, 0.05, 0.05, 0.99)),
    coords = c("x", "y"))
  feature_data <- tibble::tibble(
    name = letters[1:3],
    survey = rep(TRUE, 3),
    sensitivity = c(0.5, 0.96, 0.97),
    specificity = c(0.34, 0.92, 0.98),
    model_sensitivity = c(0.8, 0.7, 0.6),
    alpha = abs(rnorm(3)) + 1,
    gamma = runif(3))
  site_occupancy_columns <- c("f1", "f2", "f3")
  site_probability_columns <-  c("p1", "p2", "p3")
  # calculations
  pij <- prior_probability_matrix(
    site_data, feature_data, site_occupancy_columns, site_probability_columns,
    "sensitivity", "specificity", "model_sensitivity")
  # tests
  expect_is(pij, "matrix")
  expect_true(all(is.finite(pij)))
  expect_lte(max(pij), 1)
  expect_gte(min(pij), 0)
  correct <- matrix(0, ncol = nrow(site_data), nrow = nrow(feature_data))
  for (f in seq_len(nrow(feature_data))) {
    for (j in seq_len(nrow(site_data))) {
      h <- site_data[[site_occupancy_columns[f]]][j]
      m <- site_data[[site_probability_columns[f]]][j]
      if (!is.na(h)) {
        if (h >= 0.5) {
          correct[f, j] <- feature_data$sensitivity[f]
        } else {
          correct[f, j] <- 1 - feature_data$specificity[f]
        }
      } else {
        correct[f, j] <- feature_data$model_sensitivity[f] * m
      }
    }
  }
  expect_equivalent(pij, correct)
})
