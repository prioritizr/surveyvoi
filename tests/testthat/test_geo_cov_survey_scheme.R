context("geo_cov_survey_scheme")

test_that("single solution", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = c(0.1, 0.2, 0.3, 10),
                   y = c(0.1, 0.2, 0.3, 10),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- geo_cov_survey_scheme(x, "cost", 2)
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
    tibble::tibble(x = c(0.1, 0.2, 0.3, 10),
                   y = c(0.1, 0.2, 0.3, 10),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- geo_cov_survey_scheme(x, "cost", c(2, 3))
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
    tibble::tibble(x = c(0.1, 0.2, 0.5, 10),
                   y = c(0.1, 0.2, 0.5, 10),
                   locked_in = rep(FALSE, 4),
                   cost = c(1, 100, 1, 1)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- geo_cov_survey_scheme(x, "cost", 2)
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
    tibble::tibble(x = c(0.1, 0.2, 0.5, 10),
                   y = c(0.1, 0.2, 0.5, 10),
                   locked_in = c(TRUE, FALSE, TRUE, FALSE),
                   cost = c(1, 1, 1, 1)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- geo_cov_survey_scheme(x, "cost", 2, "locked_in")
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
    tibble::tibble(x = c(0.1, 0.2, 0.5, 10),
                   y = c(0.1, 0.2, 0.5, 10),
                   locked_out = c(FALSE, FALSE, FALSE, TRUE),
                   cost = c(1, 1, 1, 1)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- geo_cov_survey_scheme(x, "cost", 2, NULL, "locked_out")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, TRUE, FALSE))
})

test_that("locked out, exclude = TRUE", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = c(0.1, 0.2, 0.5, 10),
                   y = c(0.1, 0.2, 0.5, 10),
                   locked_out = c(FALSE, FALSE, FALSE, TRUE),
                   cost = c(1, 1, 1, 1)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- geo_cov_survey_scheme(x, "cost", 2, NULL, "locked_out",
                             exclude_locked_out = TRUE)
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, TRUE, FALSE))
})
