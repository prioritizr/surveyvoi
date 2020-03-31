context("rcpp sample states")

test_that("rcpp_sample_n_uniform_states_with_replacement (expected results)", {
  p <- matrix(runif(6), nrow = 3, ncol = 2)
  max_n <- rcpp_n_states(length(p))
  k <- rcpp_sample_n_uniform_states_with_replacement(4, p)
  expect_is(k, "numeric")
  expect_length(k, 4)
  expect_lte(max(k), max_n)
  expect_gte(min(k), 0)
})

test_that("rcpp_sample_n_uniform_states_without_replacement (k < n_states)", {
  p <- matrix(runif(6), nrow = 3, ncol = 2)
  k <- 4
  n <- rcpp_n_states(length(p))
  s <- character(20)
  for (i in seq_len(20)) {
    o <- rcpp_sample_n_uniform_states_without_replacement(k, p)
    s[i] <- paste(o, collapse = "")
    expect_is(o, "numeric")
    expect_length(o, k)
    expect_lte(max(o), n)
    expect_gte(min(o), 0)
    expect_equal(anyDuplicated(o), 0)
  }
  expect_equal(anyDuplicated(s), 0)
})

test_that("rcpp_sample_n_uniform_states_without_replacement (k = n_states)", {
  p <- matrix(runif(6), nrow = 3, ncol = 2)
  k <- rcpp_n_states(length(p))
  s <- character(20)
  o <- rcpp_sample_n_uniform_states_without_replacement(k + 1, p)
  expect_is(o, "numeric")
  expect_equal(sort(o), seq(0, k))
})

test_that("rcpp_sample_n_weighted_states_with_replacement (expected results)", {
  p <- matrix(runif(6), nrow = 3, ncol = 2)
  max_n <- rcpp_n_states(length(p))
  k <- rcpp_sample_n_weighted_states_with_replacement(4, p)
  expect_is(k, "numeric")
  expect_length(k, 4)
  expect_lte(max(k), max_n)
  expect_gte(min(k), 0)
})

test_that("rcpp_sample_n_weighted_states_with_replacement (correct results)", {
  p <- matrix(round(runif(6), 2), nrow = 3, ncol = 2)
  k <- rcpp_sample_n_weighted_states_with_replacement(100000, p)
  states <- lapply(k, rcpp_nth_state, p)
  p2 <- matrix(rowMeans(sapply(states, c)), nrow = nrow(p), ncol = ncol(p))
  # tests
  expect_equal(p, round(p2, 2))
})
