context("relative_site_richness_scores")

test_that("correct_result", {
  # data
  x <- tibble::tibble(x = rnorm(5), y = rnorm(5))
  x <- sf::st_as_sf(x, coords = c("x", "y"))
  for (i in seq_len(3)) x[[paste0("p", i)]] <- runif(5)
  # calculation
  r <- relative_site_richness_scores(x, paste0("p", seq_len(3)))
  # tests
  expect_is(r, "numeric")
  expect_length(r, nrow(x))
  expect_true(all(is.finite(r)))
})

test_that("sensible results (only 0.5, 0, and 1)", {
  # data
  x <- tibble::tibble(x = rnorm(5), y = rnorm(5))
  x <- sf::st_as_sf(x, coords = c("x", "y"))
  pij <- matrix(NA, nrow = 3, ncol = 5)
  pij[1, ] <- c(0.5, 1, 1, 1, 1)
  pij[2, ] <- c(0.5, 0.5, 1, 1, 1)
  pij[3, ] <- c(0.5, 0.5, 0.5, 1, 1)
  pij <- t(pij)
  for (i in seq_len(3)) x[[paste0("p", i)]] <- pij[, i]
  # calculation
  r <- relative_site_richness_scores(x, c("p1", "p2", "p3"))
  # tests
  expect_equal(r[4], r[5])
  expect_gt(r[4], r[3])
  expect_gt(r[4], r[2])
  expect_gt(r[4], r[1])
  expect_gt(r[3], r[2])
  expect_gt(r[3], r[1])
  expect_gt(r[2], r[1])
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
  r <- relative_site_richness_scores(x, c("p1", "p2", "p3"))
  # tests
  expect_equal(which.max(r), 4)
})
