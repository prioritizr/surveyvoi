context("rcpp_nth_state_sparse")

test_that("correct results", {
  # initialize variables
  m <- matrix(2, nrow = 2, ncol = 2)
  idx <- c(1, 3)
  # 0'th state
  expect_equal(rcpp_nth_state_sparse(0, m, idx),
               matrix(c(0, 2, 0, 2), nrow = 2, ncol = 2))
  # 1st state
  expect_equal(rcpp_nth_state_sparse(1, m, idx),
               matrix(c(0, 2, 1, 2), ncol = 2))
  # 3rd state
  expect_equal(rcpp_nth_state_sparse(3, m, idx),
               matrix(c(1, 2, 1, 2), ncol = 2))
})

test_that("correct full iteration", {
  # generate all possible states, we know this is 16 ahead of time
  m <- matrix(0, nrow = 2, ncol = 2)
  counter <- 0
  idx <- seq_len(4)
  out <- c()
  while (sum(m) != 4) {
    m <- rcpp_nth_state_sparse(counter, m, idx)
    out[counter + 1] <- paste(c(m), collapse = "")
    counter <- counter + 1
  }
  # tests
  expect_equal(sum(duplicated(out)), 0)
  expect_length(out, 16)
})
