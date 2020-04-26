context("prioritization")

test_that("expected result", {
  # data
  set.seed(123)
  n_pu <- 100
  n_f <- 20
  prew <- 200
  postw <- 30
  target <- 10
  n_approx_points <- 1000
  rij <- matrix(runif(n_pu * n_f), ncol = n_pu, nrow = n_f)
  pu_costs <- runif(n_pu)
  pu_locked_in <- sample(c(0, 1), n_pu, replace = TRUE, prob = c(0.8, 0.2))
  budget <- sum(pu_costs) * 0.7
  gap <- 0.0
  # results
  r1 <- r_prioritization(
    rij, pu_costs, pu_locked_in, rep(prew, n_f), rep(postw, n_f),
    rep(target, n_f), n_approx_points, budget, gap, "r1.lp")
  r2 <- rcpp_prioritization(
    rij, pu_costs, pu_locked_in, rep(prew, n_f), rep(postw, n_f),
    rep(target, n_f), n_approx_points, budget, gap, "r2.lp")
  # tests
  objval <-
    r_conservation_value_state(
      matrix(r1$x, ncol = n_pu, nrow = n_f, byrow = TRUE) * rij,
      rep(prew, n_f), rep(postw, n_f), rep(target, n_f), ncol(rij))
  expect_lte(abs(r1$objval - objval), 1e-4)
  expect_lte(abs(r2$objval - objval), 1e-4)
  expect_equal(r1$x, r2$x)
  # clean up
  unlink("r1.lp")
  unlink("r2.lp")
})

test_that("correct result", {
  # data
  set.seed(123)
  n_pu <- 8
  n_f <- 5
  prew <- 200
  postw <- 0.1
  target <- 10
  n_approx_points <- 1000
  rij <- matrix(runif(n_pu * n_f), ncol = n_pu, nrow = n_f)
  pu_costs <- runif(n_pu)
  pu_locked_in <- sample(c(0, 1), n_pu, replace = TRUE, prob = c(0.8, 0.2))
  budget <- ceiling(sum(pu_costs) * 0.7)
  gap <- 0.0
  # results
  r1 <- r_prioritization(
    rij, pu_costs, pu_locked_in, rep(prew, n_f), rep(postw, n_f),
    rep(target, n_f), n_approx_points, budget, gap, "")
  r2 <- rcpp_prioritization(
    rij, pu_costs, pu_locked_in, rep(prew, n_f), rep(postw, n_f),
    rep(target, n_f), n_approx_points, budget, gap, "")
  r3 <- brute_force_prioritization(
    rij, rep(prew, n_f), rep(postw, n_f), rep(target, n_f), pu_costs,
    pu_locked_in, budget)
  # tests
  expect_equal(r1$x, r3$x)
  expect_equal(r2$x, r3$x)
  expect_lte(abs(r1$objval - r3$objval), 1e-4)
  expect_lte(abs(r2$objval - r3$objval), 1e-4)
})

test_that("expected result (one feature has all zeros)", {
  # data
  set.seed(123)
  n_pu <- 5
  n_f <- 2
  prew <- 200
  postw <- 0.1
  target <- 10
  n_approx_points <- 1000
  rij <- matrix(0, ncol = n_pu, nrow = n_f)
  rij[1, ] <- 1
  pu_costs <- runif(n_pu)
  pu_locked_in <- sample(c(0, 1), n_pu, replace = TRUE, prob = c(0.8, 0.2))
  budget <- sum(pu_costs) * 1.1
  gap <- 0.0
  # results
  r1 <- r_prioritization(
    rij, pu_costs, pu_locked_in, rep(prew, n_f), rep(postw, n_f),
    rep(target, n_f), n_approx_points, budget, gap, "r1.lp")
  r2 <- rcpp_prioritization(
    rij, pu_costs, pu_locked_in, rep(prew, n_f), rep(postw, n_f),
    rep(target, n_f), n_approx_points, budget, gap, "r2.lp")
  # tests
  objval <-
    r_conservation_value_state(
      matrix(r1$x, ncol = n_pu, nrow = n_f, byrow = TRUE) * rij,
      rep(prew, n_f), rep(postw, n_f), rep(target, n_f), ncol(rij))
  expect_lte(abs(r1$objval - objval), 1e-4)
  expect_lte(abs(r2$objval - objval), 1e-4)
  expect_equal(r1$x, r2$x)
  expect_true(all(r1$x))
  expect_true(all(r2$x))
  # clean up
 unlink("r1.lp")
 unlink("r2.lp")
})

test_that("complex example", {
  # data
  set.seed(123)
  preweight <- c(2.348694, 2.453980) * 200
  postweight <- c(0.007071096, 0.799745903)
  target <- c(2.5, 2.5)
  rij <- matrix(c(
    0.8887591, 0.8642433, 0.1467236, 0.9039545, 0.8795562,  0.002713299,
      0.7849601, 0.14347625,
    0.9015369, 0.1301346, 0.0178074, 0.5321106, 0.8479811, 0.097441299,
      0.7894374, 0.02387501),
    nrow = 2, byrow = TRUE)
  budget <- 391.4425
  pu_costs <- c(64.00377, 66.53150, 67.59650, 54.29766, 61.80682, 57.46875,
                55.49012, 68.57491)
  gap <- 0
  pu_locked_in <- rep(0, ncol(rij))
  # results
  r1 <- r_prioritization(
    rij, pu_costs, pu_locked_in, preweight, postweight, target, 1000, budget,
    gap, "r1.lp")
  r2 <- rcpp_prioritization(
    rij, pu_costs, pu_locked_in, preweight, postweight, target, 1000, budget,
    gap, "r2.lp")
  r3 <- brute_force_prioritization(
    rij, preweight, postweight, target, pu_costs, pu_locked_in, budget)
  # tests
  expect_lte(abs(r1$objval - r3$objval), 1e-4)
  expect_lte(abs(r2$objval - r3$objval), 1e-4)
  expect_equal(r1$x, r3$x)
  expect_equal(r2$x, r3$x)
  # clean up
  unlink("r1.lp")
  unlink("r2.lp")
})
