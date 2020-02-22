context("which_state")

test_that("expected result", {
  # initialize variables
  m <- matrix(2, nrow = 2, ncol = 2)
  idx <- c(1, 3)
  # 0'th state
  expect_equal(0, which_state(matrix(c(0, 2, 0, 2), nrow = 2, ncol = 2), idx))
  # 1st state
  expect_equal(1, which_state(matrix(c(0, 2, 1, 2), ncol = 2), idx))
  # 2nd state
  expect_equal(3, which_state(matrix(c(1, 2, 1, 2), ncol = 2), idx))
})

test_that("correct full iteration", {
  # generate all possible states, we know this is 16 ahead of time
  m <- matrix(0, nrow = 2, ncol = 2)
  counter <- 0
  idx <- seq_len(4)
  out <- c()
  for (i in seq_len(15)) {
    expect_equal(which_state(nth_state(i, m, idx), idx), i)
  }
})
