context("rcpp_model_performance")

test_that("formula (Foody example)", {
  # create data
  conf_matrix <- matrix(c(120, 180, 220, 480), ncol = 2, nrow = 2)
  rse <- 0.75
  rsp <- 0.75
  # calculations
  r0 <- c(0.7, 0.7)
  r1 <- r_formula_sensitivity_and_specificity(conf_matrix, rse, rsp)
  r2 <- rcpp_formula_sensitivity_and_specificity(conf_matrix, rse, rsp)
  # tests
  expect_equal(r1, r0)
  expect_equal(r2, r0)
})

test_that("maxlik (Foody example)", {
  skip_on_os("mac")
  # create data
  conf_matrix <- matrix(c(120, 180, 220, 480), ncol = 2, nrow = 2)
  rse <- 0.75
  rsp <- 0.75
  # calculations
  r0 <- c(0.7, 0.7)
  r1 <- r_maxlik_sensitivity_and_specificity(conf_matrix, rse, rsp)
  r2 <- rcpp_maxlik_sensitivity_and_specificity(conf_matrix, rse, rsp)
  # tests
  expect_equal(r1, r0)
  expect_lte(max(abs(r2 - r0)), 1e-5)
})

test_that("formula (complex example)", {
  # create data
  conf_matrix <- matrix(c(1005, 134, 147, 1668), ncol = 2, nrow = 2)
  rse <- 0.98
  rsp <- 0.85
  # calculations
  r0 <- c(1.1958615, 0.9294045)
  r1 <- r_formula_sensitivity_and_specificity(conf_matrix, rse, rsp)
  r2 <- rcpp_formula_sensitivity_and_specificity(conf_matrix, rse, rsp)
  # tests
  expect_lte(max(abs(r1 - r0)), 1e-4)
  expect_lte(max(abs(r2 - r0)), 1e-4)
})

test_that("maxlik (complex example)", {
  skip_on_os("mac")
  # create data
  conf_matrix <- matrix(c(1005, 134, 147, 1668), ncol = 2, nrow = 2)
  rse <- 0.98
  rsp <- 0.85
  # calculations
  r0 <- c(1.0000000, 0.9234635)
  r1 <- r_maxlik_sensitivity_and_specificity(conf_matrix, rse, rsp)
  r2 <- rcpp_maxlik_sensitivity_and_specificity(conf_matrix, rse, rsp)
  # tests
  expect_lte(max(abs(r1 - r0)), 1e-4)
  expect_lte(max(abs(r2 - r0)), 1e-4)
})

test_that("model_performance (Foody example)", {
  # create data
  conf_matrix <- matrix(c(120, 180, 220, 480), ncol = 2, nrow = 2)
  rse <- 0.75
  rsp <- 0.75
  y <- c(rep(1, sum(conf_matrix[, 1])), rep(0, sum(conf_matrix[, 2])))
  yhat <- c(rep(1, conf_matrix[1, 1]), rep(0, conf_matrix[2, 1]),
            rep(0, conf_matrix[2, 2]), rep(1, conf_matrix[1, 2]))
  w <- rep(1, length(y))
  # calculations
  r0 <- c(0.4, 0.7, 0.7)
  r1 <- r_model_performance(y, yhat, w, rse, rsp)
  r2 <- rcpp_model_performance(y, yhat, w, rse, rsp)
  # tests
  expect_lte(max(abs(r0 - r1)), 1e-4)
  expect_lte(max(abs(r0 - r2)), 1e-4)
})

test_that("model_performance (complex example)", {
  # create data
  conf_matrix <- matrix(c(1005, 134, 147, 1668), ncol = 2, nrow = 2)
  rse <- 0.98
  rsp <- 0.85
  y <- c(rep(1, sum(conf_matrix[, 1])), rep(0, sum(conf_matrix[, 2])))
  yhat <- c(rep(1, conf_matrix[1, 1]), rep(0, conf_matrix[2, 1]),
            rep(0, conf_matrix[2, 2]), rep(1, conf_matrix[1, 2]))
  w <- rep(1, length(y))
  # calculations
  r0 <- c(0.9234635, 1.0000000, 0.9234635)
  r1 <- r_model_performance(y, yhat, w, rse, rsp)
  r2 <- rcpp_model_performance(y, yhat, w, rse, rsp)
  # tests
  expect_lte(max(abs(r0 - r1)), 1e-4)
  expect_lte(max(abs(r0 - r2)), 1e-4)
})
