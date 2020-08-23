context("create_site_folds")

test_that("expected result", {
  # data
  n_det <- as.integer(c(0, 1, 2, 3, 0))
  n_non_det <- as.integer(c(3, 2, 1, 3, 0))
  n_total <- n_det + n_non_det
  prop_detected <- n_det / n_total
  prop_detected[!is.finite(prop_detected)] <- 0
  id <- sample(seq_len(5) + 100)
  n <- 2
  # calculations
  f <- create_site_folds(prop_detected, n_total, n, id)
  # tests
  expect_is(f, "list")
  expect_length(f, 2)
  expect_equal(names(f), c("train", "test"))
  for (i in seq_along(f$train)) {
    expect_gte(sum(prop_detected[id %in% f$train[[i]]] > 0), 1)
    expect_gte(sum(prop_detected[id %in% f$train[[i]]] < 1), 1)
  }
  for (i in seq_along(f$test)) {
    expect_gte(sum(prop_detected[id %in% f$test[[i]]] > 0), 1)
    expect_gte(sum(prop_detected[id %in% f$test[[i]]] < 1), 1)
  }
})

test_that("lots of sites with no surveys", {
  # data
  n_det <- as.integer(c(0, 1, 2, 3, 0, 0, 0, 0))
  n_non_det <- as.integer(c(3, 2, 1, 3, 0, 0, 0, 0))
  id <- sample(seq_along(n_det) + 100)
  n_total <- n_det + n_non_det
  prop_detected <- n_det / n_total
  prop_detected[!is.finite(prop_detected)] <- 0
  n <- 3
  # calculations
  f <- create_site_folds(prop_detected, n_total, n, id)
  # tests
  expect_is(f, "list")
  expect_length(f, 2)
  expect_equal(names(f), c("train", "test"))
  for (i in seq_along(f$train)) {
    expect_gte(sum(prop_detected[id %in% f$train[[i]]] > 0), 1)
    expect_gte(sum(prop_detected[id %in% f$train[[i]]] < 1), 1)
  }
  for (i in seq_along(f$test)) {
    expect_gte(sum(prop_detected[id %in% f$test[[i]]] > 0), 1)
    expect_gte(sum(prop_detected[id %in% f$test[[i]]] < 1), 1)
  }
  expect_true(
    all(id %in% unique(unlist(f, recursive = TRUE, use.names = FALSE))))
})
