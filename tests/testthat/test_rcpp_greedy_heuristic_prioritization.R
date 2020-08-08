context("rcpp_greedy_heuristic_prioritization")

test_that("expected result", {
  # data
  set.seed(123)
  n_pu <- 5
  n_f <- 2
  target <- rep(1, n_f)
  rij <- matrix(runif(n_pu * n_f), ncol = n_pu, nrow = n_f)
  pu_costs <- runif(n_pu)
  pu_locked_in <- sample(c(0, 1), n_pu, replace = TRUE, prob = c(0.8, 0.2))
  pu_locked_out <- sample(c(0, 1), n_pu, replace = TRUE, prob = c(0.8, 0.2))
  # pu_costs[sample(which(pu_locked_in < 0.5), 5)] <- 5000
  pu_locked_in[as.logical(pu_locked_out)] <- 0
  budget <- sum(pu_costs) * 0.7
  # results
  r1 <- r_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  r2 <- rcpp_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  # tests
  expect_equal(r1$x, r2$x)
  expect_equal(r1$objval, r2$objval)
})

test_that("correct result", {
  # data
  set.seed(123)
  n_pu <- 20
  n_f <- 5
  target <- rep(3, n_f)
  rij <- matrix(runif(n_pu * n_f), ncol = n_pu, nrow = n_f)
  pu_costs <- runif(n_pu)
  pu_locked_in <- sample(c(0, 1), n_pu, replace = TRUE, prob = c(0.8, 0.2))
  pu_locked_out <- sample(c(0, 1), n_pu, replace = TRUE, prob = c(0.8, 0.2))
  pu_costs[sample(which(pu_locked_in < 0.5), 2)] <- 5000
  pu_locked_in[as.logical(pu_locked_out)] <- 0
  budget <- ceiling(sum(pu_costs) * 0.7)
  # results
  r1 <- r_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  r2 <- rcpp_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  r3 <- brute_force_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  # tests
  expect_equal(r1$x, r3$x)
  expect_equal(r2$x, r3$x)
  expect_lte(abs(r1$objval - r2$objval), 1e-5)
  expect_lte(abs(r1$objval - r3$objval), 1e-5)
})

test_that("expected result (one feature has all 0s)", {
  # data
  set.seed(123)
  n_pu <- 15
  n_f <- 2
  target <- rep(10, n_f)
  rij <- matrix(0, ncol = n_pu, nrow = n_f)
  rij[1, ] <- 1
  pu_costs <- runif(n_pu)
  pu_locked_in <- sample(c(0, 1), n_pu, replace = TRUE, prob = c(0.8, 0.2))
  pu_locked_out <- rep(0, n_pu)
  budget <- sum(pu_costs) * 1.1
  # results
  r1 <- r_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  r2 <- rcpp_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  r3 <- brute_force_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  # tests
  expect_equal(r1$x, r2$x)
  expect_true(all(r1$x))
  expect_true(all(r2$x))
  expect_lte(abs(r1$objval - r2$objval), 1e-5)
  expect_lte(abs(r1$objval - r3$objval), 1e-4)
})

test_that("complex example", {
  # data
  set.seed(123)
  target <- c(2, 2)
  rij <- matrix(c(
    0.8887591, 0.8642433, 0.1467236, 0.9039545, 0.8795562,  0.002713299,
    0.7849601, 0.14347625,
    0.9015369, 0.1301346, 0.0178074, 0.5321106, 0.8479811, 0.097441299,
    0.7894374, 0.02387501),
    nrow = 2, byrow = TRUE)
  budget <- 391.4425
  pu_costs <- c(64.00377, 66.53150, 67.59650, 54.29766, 61.80682, 57.46875,
                55.49012, 68.57491)
  pu_locked_in <- rep(0, ncol(rij))
  pu_locked_out <- rep(0, ncol(rij))
  # results
  r1 <- r_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  r2 <- rcpp_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  r3 <- brute_force_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  # tests
  expect_equal(r1$x, r2$x)
  expect_equal(r1$x, r3$x)
  expect_lte(abs(r1$objval - r2$objval), 1e-4)
  expect_lte(abs(r1$objval - r3$objval), 1e-4)
})

test_that("highly variable costs", {
  # data
  set.seed(123)
  n_sites <- 12
  n_f <- 2
  target <- c(2, n_f)
  rij <- matrix(runif(n_sites * n_f), ncol = n_sites, nrow = n_f)
  pu_costs <- runif(n_sites)
  budget <- sum(pu_costs) * 0.3
  pu_costs[c(3:7)] <- 1e+5
  pu_locked_in <- rep(0, ncol(rij))
  pu_locked_out <- rep(0, ncol(rij))
  # results
  r1 <- r_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  r2 <- rcpp_greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  r3 <- brute_force_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget)
  # tests
  expect_equal(r1$x, r2$x)
  expect_lte(abs(r1$objval - r2$objval), 1e-4)
  expect_gte(mean(r1$x == r3$x), 0.6) # expect similar solutions
  expect_lte((r3$objval - r1$objval) / r3$objval, 0.15) # expect <15% optimality
})
