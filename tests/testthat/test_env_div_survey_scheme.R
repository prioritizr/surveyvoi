context("env_div_survey_scheme")

test_that("single solution (gurobi)", {
  skip_if_not_installed("gurobi")
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
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean", solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("multiple solutions (gurobi)", {
  skip_if_not_installed("gurobi")
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(5),
                   y = rnorm(5),
                   v1 = c(0.1, 0.21, 0.22, 0.23, 10),
                   v2 = c(0.1, 0.21, 0.22, 0.23, 10),
                   locked_in = rep(FALSE, 5),
                   cost = rep(1, 5)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- env_div_survey_scheme(
    x, "cost", c(2, 3), c("v1", "v2"), "euclidean", solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 2)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, FALSE, TRUE, FALSE, TRUE))
  expect_equal(c(r[2, ]), c(TRUE, FALSE, TRUE, FALSE, TRUE))
})

test_that("variable costs (gurobi)", {
  skip_if_not_installed("gurobi")
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
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean", solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(TRUE, FALSE, FALSE, TRUE))
})

test_that("locked in (gurobi)", {
  skip_if_not_installed("gurobi")
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
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean", "locked_in", solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), x$locked_in)
})

test_that("locked out (gurobi)", {
  skip_if_not_installed("gurobi")
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   v1 = c(0.1, 0.21, 0.5, 10),
                   v2 = c(0.1, 0.21, 0.5, 10),
                   locked_out = c(FALSE, FALSE, FALSE, TRUE),
                   cost = c(1, 1, 1, 1)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean", NULL, "locked_out",
    solver = "gurobi")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(TRUE, FALSE, TRUE, FALSE))
})

test_that("locked out, exclude = TRUE (gurobi)", {
  skip_if_not_installed("gurobi")
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
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean", NULL, "locked_out",
    exclude_locked_out = TRUE, solver = "gurobi")
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
                   v1 = c(0.1, 0.2, 0.3, 10),
                   v2 = c(0.1, 0.2, 0.3, 10),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("multiple solutions (Rsymphony)", {
  # data
  x <- sf::st_as_sf(
    tibble::tibble(x = rnorm(4),
                   y = rnorm(4),
                   v1 = c(0.1, 0.2, 0.21, 10),
                   v2 = c(0.1, 0.2, 0.21, 10),
                   locked_in = rep(FALSE, 4),
                   cost = rep(1, 4)),
    coords = c("x", "y"))
  # generate prioritisation
  r <- env_div_survey_scheme(
    x, "cost", c(2, 3), c("v1", "v2"), "euclidean", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 2)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r[1, ]), c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(c(r[2, ]), c(TRUE, FALSE, TRUE, TRUE))
})

test_that("variable costs (Rsymphony)", {
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
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(TRUE, FALSE, FALSE, TRUE))
})

test_that("locked in (Rsymphony)", {
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
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean",
    "locked_in", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), x$locked_in)
})

test_that("locked out (Rsymphony)", {
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
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean",
   NULL, "locked_out", solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, TRUE, FALSE))
})

test_that("locked out, exclude = TRUE (Rsymphony)", {
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
  r <- env_div_survey_scheme(
    x, "cost", 2, c("v1", "v2"), "euclidean",
    NULL, "locked_out", exclude_locked_out = TRUE,
    solver = "Rsymphony")
  # tests
  expect_is(r, "matrix")
  expect_equal(nrow(r), 1)
  expect_equal(ncol(r), nrow(x))
  expect_is(r[1], "logical")
  expect_equal(c(r), c(FALSE, TRUE, TRUE, FALSE))
})
