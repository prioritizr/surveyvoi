context("rcpp_k_nth_state")

test_that("correct results", {
  p <- matrix(runif(6), nrow = 3, ncol = 2)
  max_n <- rcpp_n_states(length(p))
  k <- rcpp_sample_k_nth_states(4, p)
  expect_is(k, "numeric")
  expect_length(k, 4)
  expect_lte(max(k), max_n)
  expect_gte(min(k), 0)
})

test_that("expected results", {
  p <- matrix(round(runif(6), 2), nrow = 3, ncol = 2)
  k <- rcpp_sample_k_nth_states(100000, p)
  states <- lapply(k, rcpp_nth_state, p)
  p2 <- matrix(rowMeans(sapply(states, c)), nrow = nrow(p), ncol = ncol(p))
  # tests
  expect_equal(p, round(p2, 2))
})
