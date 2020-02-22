context("env_div_survey_scheme")

test_that("single solution", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   v1 = c(0.1, 0.2, 0.3, 10),
                   v2 = c(0.1, 0.2, 0.3, 10),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- env_div_survey_scheme(x, "cost", 2, c("v1", "v2"), "euclidean")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("multiple solutions", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   v1 = c(0.1, 0.2, 0.3, 10),
                   v2 = c(0.1, 0.2, 0.3, 10),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- env_div_survey_scheme(x, "cost", c(2, 3), c("v1", "v2"), "euclidean")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 2)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(c(r[2, ]), c(TRUE, FALSE, TRUE, TRUE))
})

test_that("variable costs", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   v1 = c(0.1, 0.2, 0.5, 10),
                   v2 = c(0.1, 0.2, 0.5, 10),
                   locked_in = rep(FALSE, 4),
                   cost = c(1, 100, 1, 1)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- env_div_survey_scheme(x, "cost", 2, c("v1", "v2"), "euclidean")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(TRUE, FALSE, FALSE, TRUE))
})

test_that("locked in", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   v1 = c(0.1, 0.2, 0.5, 10),
                   v2 = c(0.1, 0.2, 0.5, 10),
                   locked_in = c(TRUE, FALSE, TRUE, FALSE),
                   cost = c(1, 1, 1, 1)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- env_div_survey_scheme(x, "cost", 2, c("v1", "v2"), "euclidean",
                             "locked_in")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), x$locked_in)
})

test_that("locked out", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   v1 = c(0.1, 0.2, 0.5, 10),
                   v2 = c(0.1, 0.2, 0.5, 10),
                   locked_out = c(FALSE, FALSE, FALSE, TRUE),
                   cost = c(1, 1, 1, 1)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- env_div_survey_scheme(x, "cost", 2, c("v1", "v2"), "euclidean",
                             NULL, "locked_out")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, TRUE, FALSE))
})
