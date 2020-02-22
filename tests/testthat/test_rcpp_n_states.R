context("rcpp_n_states")

# define helper function
correct_n_states <- function(x) {
  sum(vapply(seq_len(x), function(i) length(combn(x, i, simplify = FALSE)),
             numeric(1)))
}

test_that("correct result", {
  # 1 planning unit and 1 feature
  expect_equal(rcpp_n_states(1), correct_n_states(1))
  # 1 planning unit and 2 features
  expect_equal(rcpp_n_states(2), correct_n_states(2))
  # 2 planning unit and 2 features
  expect_equal(rcpp_n_states(4), correct_n_states(4))
  # 3 planning units and 2 features
  expect_equal(rcpp_n_states(6), correct_n_states(6))
  # verify several numbers
  for (i in seq_len(20))
    expect_equal(rcpp_n_states(i), correct_n_states(i))
})
