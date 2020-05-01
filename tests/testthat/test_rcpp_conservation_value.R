context("rcpp_conservation_value")

test_that("rcpp_conservation_value_amount", {
  # data
  set.seed(500)
  prew <- 200
  postw <- 5
  target <- 30
  total <- 100
  amount <- seq(0, 100)
  # calculations
  r1 <- sapply(amount, r_conservation_value_amount, preweight = prew,
    postweight = postw, target = target, total = total)
  r2 <- sapply(amount, rcpp_conservation_value_amount, preweight = prew,
    postweight = postw, target = target, total = total)
  # tests
  expect_equal(r1, r2)
})

test_that("rcpp_conservation_value_states", {
  # data
  set.seed(500)
  prew <- runif(3, 100, 200)
  postw <- runif(3, 5, 10)
  target <- ceiling(runif(3, 5, 20))
  m <- matrix(0, ncol = 5, nrow = 3)
  states <- lapply(seq_len(n_states(5, 3)), function(i) rcpp_nth_state(i, m))
  # calculations
  r1 <- sapply(states, r_conservation_value_state, preweight = prew,
    postweight = postw, target = target, total = ncol(states[[1]]))
  r2 <- sapply(states, rcpp_conservation_value_state, preweight = prew,
    postweight = postw, target = target, total = ncol(states[[1]]))
  # tests
  expect_equal(r1, r2)
})
