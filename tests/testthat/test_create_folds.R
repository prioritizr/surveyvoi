context("create_folds")

test_that("expected result", {
  # data
  y <- c(rep(0, 5), rep(1, 5))
  # calculations
  f <- create_folds(y, 3)
  # tests
  expect_is(f, "list")
  expect_length(f, 2)
  expect_equal(names(f), c("train", "test"))
  expect_equivalent(lengths(f), rep(3, 2))
  for (i in seq_along(f$train)) {
    expect_gte(sum(y[f$train[[i]]] > 0.5), 1)
    expect_gte(sum(y[f$train[[i]]] < 0.5), 1)
  }
  for (i in seq_along(f$test)) {
    expect_gte(sum(y[f$test[[i]]] > 0.5), 1)
    expect_gte(sum(y[f$test[[i]]] < 0.5), 1)
  }
})

test_that("few presences and absences", {
  # data
  y <- c(0, 1)
  # calculations
  f <- create_folds(y, 10)
  # tests
  expect_is(f, "list")
  expect_length(f, 2)
  expect_equal(names(f), c("train", "test"))
  expect_equivalent(lengths(f), rep(10, 2))
  for (i in seq_along(f$train)) {
    expect_gte(sum(y[f$train[[i]]] > 0.5), 1)
    expect_gte(sum(y[f$train[[i]]] < 0.5), 1)
  }
  for (i in seq_along(f$test)) {
    expect_gte(sum(y[f$test[[i]]] > 0.5), 1)
    expect_gte(sum(y[f$test[[i]]] < 0.5), 1)
  }
})

test_that("missing values", {
  # data
  y <- c(0, 1, NA)
  # calculations
  f <- create_folds(y, 10, na.fail = FALSE)
  # tests
  expect_is(f, "list")
  expect_length(f, 2)
  expect_equal(names(f), c("train", "test"))
  expect_equivalent(lengths(f), rep(10, 2))
  for (i in seq_along(f$train)) {
    expect_gte(sum(y[f$train[[i]]] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(y[f$train[[i]]] < 0.5, na.rm = TRUE), 1)
  }
  for (i in seq_along(f$test)) {
    expect_gte(sum(y[f$test[[i]]] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(y[f$test[[i]]] < 0.5, na.rm = TRUE), 1)
  }
  expect_true(3 %in% unlist(f))
})

test_that("custom indices", {
  # data
  x <- c(0, 1, NA)
  index <- c(14, 100, 300)
  # calculations
  f <- create_folds(x, 3, index = index, na.fail = FALSE)
  # tests
  expect_is(f, "list")
  expect_length(f, 2)
  expect_equal(names(f), c("train", "test"))
  expect_equivalent(lengths(f), rep(3, 2))
  for (i in seq_along(f$train)) {
    idx <- match(f$train[[i]], index)
    expect_gte(sum(x[idx] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(x[idx] < 0.5, na.rm = TRUE), 1)
  }
  for (i in seq_along(f$test)) {
    idx <- match(f$test[[i]], index)
    expect_gte(sum(x[idx] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(x[idx] < 0.5, na.rm = TRUE), 1)
  }
  expect_true(setequal(index, unlist(f)))
})
