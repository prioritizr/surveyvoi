context("functions")

test_that("rcpp_log_sum", {
  expect_equal(rcpp_log_sum(c(log(3), log(4))), log(3 + 4))
})
