context("rcpp_conservation_value")

test_that("rcpp_conservation_value_amount", {
  # data
  set.seed(123)
  n_pu <- 10
  n_f <- 3
  target <- rep(1, n_f)
  rij <- matrix(runif(n_pu * n_f), ncol = n_pu, nrow = n_f)
  status <- matrix(0, ncol = n_pu)
  sols <- lapply(seq_len(n_states(n_pu - 1, 1)), function(i)
    rcpp_nth_state(i, status))
  # calculations
  r1 <- sapply(sols, function(x) {
    r_conservation_value(rij[, which(c(x) > 0.5), drop = FALSE], target)
  })
  r2 <- sapply(sols, function(x) {
    rcpp_expected_value_of_action(as.logical(c(x)), rij, target)
  })
  # tests
  expect_equal(r1, r2)
})

test_that("rcpp_conservation_value_states", {
  # data
  set.seed(500)
  target <- rep(1, 3)
  m <- matrix(0, ncol = 5, nrow = 3)
  states <- lapply(seq_len(n_states(5, 3)), function(i) rcpp_nth_state(i, m))
  # calculations
  r1 <- sapply(states, r_conservation_value, target = target)
  r2 <- sapply(states, rcpp_expected_value_of_action,
               solution = rep(TRUE, ncol(m)), target = target)
  # tests
  expect_equal(r1, r2)
})
