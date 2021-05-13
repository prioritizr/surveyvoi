context("weighted_survey_scheme")

test_that("single solution (gurobi)", {
  skip_if_not_installed("gurobi")
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(5, 10, 8, 1),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", 2, "w", solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, TRUE, FALSE))
})

test_that("multiple solutions (gurobi)", {
  skip_if_not_installed("gurobi")
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(5, 10, 8, 1),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", c(1, 2, 3, 4), "w", solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 4)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(r[1, ], c(FALSE, TRUE, FALSE, FALSE))
  expect_equal(r[2, ], c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(r[3, ], c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(r[4, ], c(TRUE, TRUE, TRUE, TRUE))
})

test_that("variable costs (gurobi)", {
  skip_if_not_installed("gurobi")
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(5, 10, 8, 1),
                   locked_in = rep(FALSE, 4),
                   cost = c(1, 100, 1, 0.1)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", c(2, 3), "w", solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 2)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(r[1, ], c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(r[2, ], c(TRUE, FALSE, TRUE, TRUE))
})

test_that("locked in (gurobi)", {
  skip_if_not_installed("gurobi")
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(0.01, 10, 8, 1),
                   locked_in = c(TRUE, FALSE, FALSE, FALSE),
                   cost = c(1, 1, 1, 1)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", c(2, 3), "w", "locked_in", solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 2)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(r[1, ], c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(r[2, ], c(TRUE, TRUE, TRUE, FALSE))
})

test_that("locked_out (gurobi)", {
  skip_if_not_installed("gurobi")
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(5, 10, 8, 1),
                   locked_out = c(TRUE, FALSE, FALSE, TRUE),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", 2, "w", NULL, "locked_out", solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, TRUE, FALSE))
})

test_that("single solution (Rsymphony)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(5, 10, 8, 1),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", 2, "w", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, TRUE, FALSE))
})

test_that("multiple solutions (Rsymphony)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(5, 10, 8, 1),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", c(1, 2, 3, 4), "w", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 4)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(r[1, ], c(FALSE, TRUE, FALSE, FALSE))
  expect_equal(r[2, ], c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(r[3, ], c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(r[4, ], c(TRUE, TRUE, TRUE, TRUE))
})

test_that("variable costs (Rsymphony)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(5, 10, 8, 1),
                   locked_in = rep(FALSE, 4),
                   cost = c(1, 100, 1, 0.1)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", c(2, 3), "w", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 2)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(r[1, ], c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(r[2, ], c(TRUE, FALSE, TRUE, TRUE))
})

test_that("locked in (Rsymphony)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(0.01, 10, 8, 1),
                   locked_in = c(TRUE, FALSE, FALSE, FALSE),
                   cost = c(1, 1, 1, 1)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", c(2, 3), "w", "locked_in", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 2)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(r[1, ], c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(r[2, ], c(TRUE, TRUE, TRUE, FALSE))
})

test_that("locked_out (Rsymphony)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   w = c(5, 10, 8, 1),
                   locked_out = c(TRUE, FALSE, FALSE, TRUE),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate scheme
  r <- weighted_survey_scheme(
    x, "cost", 2, "w", NULL, "locked_out", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, TRUE, FALSE))
})
