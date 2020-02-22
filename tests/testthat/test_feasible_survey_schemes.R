context("feasible_survey_schemes")

test_that("equal costs (no locked sites)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # calculations
  r <- feasible_survey_schemes(x, "cost", 4)
  # tests
  expect_is(r, "matrix")
  expect_is(r[1], "logical")
  expect_equal(ncol(r), nrow(x))
  expect_equal(nrow(r), 16)
  expect_equal(anyDuplicated(apply(r, 1, paste, collapse = "")), 0)
})

test_that("equal costs (locked in sites)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   cost = rep(1, 4),
                   locked_in = c(TRUE, TRUE, FALSE, TRUE)),
    coords = c("x", "y"))
  # calculations
  r <- feasible_survey_schemes(x, "cost", 4, "locked_in")
  r <- r[order(apply(r, 1, paste, collapse = "")), , drop = FALSE]
  r2 <- matrix(c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
               ncol = nrow(x))
  # tests
  expect_equal(r, r2)
})

test_that("equal costs (locked out sites)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   cost = rep(1, 4),
                   locked_out = c(TRUE, FALSE, FALSE, TRUE)),
    coords = c("x", "y"))
  # calculations
  r <- feasible_survey_schemes(x, "cost", 4, locked_out = "locked_out")
  r <- r[order(apply(r, 1, paste, collapse = "")), , drop = FALSE]
  r2 <- matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE,
                 TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
               ncol = nrow(x))
  # tests
  expect_equal(r, r2)
})


test_that("equal costs (all sites fixed)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   cost = rep(1, 4),
                   locked_out = c(TRUE, FALSE, FALSE, TRUE),
                   locked_in = c(FALSE, TRUE, TRUE, FALSE)),
    coords = c("x", "y"))
  # calculations
  r <- feasible_survey_schemes(x, "cost", 4, locked_in = "locked_in",
                               locked_out = "locked_out")
  r <- r[order(apply(r, 1, paste, collapse = "")), , drop = FALSE]
  r2 <- matrix(c(FALSE, TRUE, TRUE, FALSE), ncol = nrow(x))
  # tests
  expect_equal(r, r2)
})

test_that("equal costs (one site fixed)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   cost = rep(1, 4),
                   locked_in = c(FALSE, TRUE, FALSE, FALSE),
                   locked_out = c(TRUE, FALSE, FALSE, TRUE)),
    coords = c("x", "y"))
  # calculations
  r <- feasible_survey_schemes(x, "cost", 4, locked_in = "locked_in",
                               locked_out = "locked_out")
  r <- r[order(apply(r, 1, paste, collapse = "")), , drop = FALSE]
  r2 <- matrix(c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
               ncol = nrow(x))
  # tests
  expect_equal(r, r2)
})

test_that("variable costs (no locked sites)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   cost = c(100, 200, 0.2, 1)),
    coords = c("x", "y"))
  # calculations
  r <- feasible_survey_schemes(x, "cost", 4)
  r <- r[order(apply(r, 1, paste, collapse = "")), , drop = FALSE]
  r2 <- matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                 TRUE, TRUE, TRUE, FALSE, TRUE), ncol = nrow(x))
  r2 <- rbind(FALSE, r2)
  # tests
  expect_equal(r, r2)
})

test_that("variable costs (locked in sites)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   cost = c(100, 200, 0.2, 1)),
    coords = c("x", "y"))
  # calculations
  r <- feasible_survey_schemes(x, "cost", 4)
  r <- r[order(apply(r, 1, paste, collapse = "")), , drop = FALSE]
  r2 <- matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                 TRUE, TRUE, TRUE, FALSE, TRUE), ncol = nrow(x))
  r2 <- rbind(FALSE, r2)
  # tests
  expect_equal(r, r2)
})

test_that("variable costs (locked out sites)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   cost = c(100, 200, 0.2, 1)),
    coords = c("x", "y"))
  # calculations
  r <- feasible_survey_schemes(x, "cost", 4)
  r <- r[order(apply(r, 1, paste, collapse = "")), , drop = FALSE]
  r2 <- matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                 TRUE, TRUE, TRUE, FALSE, TRUE), ncol = nrow(x))
  r2 <- rbind(FALSE, r2)
  # tests
  expect_equal(r, r2)
})
