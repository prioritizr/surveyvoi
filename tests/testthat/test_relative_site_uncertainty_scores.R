context("relative_site_uncertainty_scores")

# https://stat.ethz.ch/pipermail/r-help/2008-July/167112.html
shannon_entropy <- function(x, names) {
  shen <- function(p) {
    if (min(p) < 0 || sum(p) <= 0)
      return(NA)
    p.norm <- p[p>0]/sum(p)
    -sum(log2(p.norm)*p.norm)
  }
  x <- as.matrix(st_drop_geometry(x)[, names])
  for (i in seq_len(length(x)))
    x[i] <- shen(c(x[i], 1 - x[i]))
  rowSums(x)
}

test_that("correct_result", {
  # data
  x <- tibble::tibble(x = rnorm(5), y = rnorm(5))
  x <- sf::st_as_sf(x, coords = c("x", "y"))
  for (i in seq_len(3)) x[[paste0("p", i)]] <- runif(5)
  # calculation
  r1 <- relative_site_uncertainty_scores(x, paste0("p", seq_len(3)))
  r2 <- shannon_entropy(x, paste0("p", seq_len(3)))
  # tests
  expect_is(r1, "numeric")
  expect_length(r1, nrow(x))
  expect_true(all(is.finite(r1)))
  expect_equal(r1, r2)
})

test_that("sensible results (only 0.5, 0, and 1)", {
  # data
  x <- tibble::tibble(x = rnorm(5), y = rnorm(5))
  x <- sf::st_as_sf(x, coords = c("x", "y"))
  pij <- matrix(NA, nrow = 3, ncol = 5)
  pij[1, ] <- c(0.5, 0, 1, 0, 1)
  pij[2, ] <- c(0.5, 0.5, 1, 0, 1)
  pij[3, ] <- c(0.5, 0.5, 0.5, 0, 1)
  pij <- t(pij)
  for (i in seq_len(3)) x[[paste0("p", i)]] <- pij[, i]
  # calculation
  r1 <- relative_site_uncertainty_scores(x, c("p1", "p2", "p3"))
  r2 <- shannon_entropy(x, c("p1", "p2", "p3"))
  # tests
  expect_gt(r1[1], r2[2])
  expect_gt(r1[1], r2[3])
  expect_gt(r1[2], r2[3])
  expect_equal(r1, r2)
})

test_that("sensible results (variable numbers)", {
  # data
  x <- tibble::tibble(x = rnorm(4), y = rnorm(4))
  x <- sf::st_as_sf(x, coords = c("x", "y"))
  pij <- matrix(NA, nrow = 3, ncol = 4)
  pij[, 1] <- runif(3, 0.01, 0.1)
  pij[, 2] <- runif(3, 0.01, 0.1)
  pij[, 3] <- c(0.5, 0.4, 0.9)
  pij[, 4] <- runif(3, 0.9, 0.99)
  pij <- t(pij)
  for (i in seq_len(3)) x[[paste0("p", i)]] <- pij[, i]
  # calculation
  r1 <- relative_site_uncertainty_scores(x, c("p1", "p2", "p3"))
  r2 <- shannon_entropy(x, c("p1", "p2", "p3"))
  # tests
  expect_equal(which.max(r1), 3)
  expect_equal(r1, r2)
})
