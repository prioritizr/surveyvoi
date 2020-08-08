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
  for (i in seq_along(f$train)) {
    expect_gte(sum(y[f$train[[i]]] > 0.5), 1)
    expect_gte(sum(y[f$train[[i]]] < 0.5), 1)
  }
  for (i in seq_along(f$test)) {
    expect_gte(sum(y[f$test[[i]]] > 0.5), 1)
    expect_gte(sum(y[f$test[[i]]] < 0.5), 1)
  }
})

test_that("a single missing value", {
  # data
  y <- c(rep(0, 10), rep(1, 5), rep(NA_real_, 1))
  # calculations
  f <- create_folds(y, 4, na.fail = FALSE)
  # tests
  expect_is(f, "list")
  expect_length(f, 2)
  expect_equal(names(f), c("train", "test"))
  for (i in seq_along(f$train)) {
    expect_gte(sum(y[f$train[[i]]] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(y[f$train[[i]]] < 0.5, na.rm = TRUE), 1)
  }
  for (i in seq_along(f$test)) {
    expect_gte(sum(y[f$test[[i]]] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(y[f$test[[i]]] < 0.5, na.rm = TRUE), 1)
  }
  expect_true(which(is.na(y)) %in% unlist(f))
})

test_that("lots of missing values", {
  # data
  y <- c(rep(0, 10), rep(1, 5), rep(NA_real_, 20))
  # calculations
  f <- create_folds(y, 3, na.fail = FALSE)
  # tests
  expect_is(f, "list")
  expect_length(f, 2)
  expect_equal(names(f), c("train", "test"))
  for (i in seq_along(f$train)) {
    expect_gte(sum(y[f$train[[i]]] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(y[f$train[[i]]] < 0.5, na.rm = TRUE), 1)
  }
  for (i in seq_along(f$test)) {
    expect_gte(sum(y[f$test[[i]]] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(y[f$test[[i]]] < 0.5, na.rm = TRUE), 1)
  }
  expect_true(all(which(is.na(y)) %in% unlist(f)))
})

test_that("custom indices", {
  # data
  y <- c(rep(0, 10), rep(1, 5), rep(NA_real_, 20))
  index <- sample.int(1000, length(y))
  # calculations
  f <- create_folds(y, 3, index = index, na.fail = FALSE)
  # tests
  expect_is(f, "list")
  expect_length(f, 2)
  expect_equal(names(f), c("train", "test"))
  for (i in seq_along(f$train)) {
    idx <- match(f$train[[i]], index)
    expect_gte(sum(y[idx] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(y[idx] < 0.5, na.rm = TRUE), 1)
  }
  for (i in seq_along(f$test)) {
    idx <- match(f$test[[i]], index)
    expect_gte(sum(y[idx] > 0.5, na.rm = TRUE), 1)
    expect_gte(sum(y[idx] < 0.5, na.rm = TRUE), 1)
  }
  expect_true(setequal(index, unlist(f)))
})
