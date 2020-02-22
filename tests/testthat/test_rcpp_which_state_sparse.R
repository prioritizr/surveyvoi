context("rcpp_which_state_sparse")

test_that("expected result", {
  # initialize variables
  m <- matrix(2, nrow = 2, ncol = 2)
  idx <- c(1, 3)
  # 0'th state
  expect_equal(0, rcpp_which_state_sparse(matrix(c(0, 2, 0, 2), nrow = 2,
                                                 ncol = 2), idx))
  # 1st state
  expect_equal(1, rcpp_which_state_sparse(matrix(c(0, 2, 1, 2), ncol = 2), idx))
  # 2nd state
  expect_equal(3, rcpp_which_state_sparse(matrix(c(1, 2, 1, 2), ncol = 2), idx))
})

test_that("correct full iteration", {
  # generate all possible states, we know this is 16 ahead of time
  m <- matrix(0, nrow = 2, ncol = 2)
  counter <- 0
  idx <- seq_len(4)
  out <- c()
  for (i in seq_len(15)) {
    expect_equal(
      rcpp_which_state_sparse(rcpp_nth_state_sparse(i, m, idx), idx), i)
  }
})
