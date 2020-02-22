context("tune_occupancy_models")

test_that("single species", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## simulate data
  n_f <- 1
  n_vars <- 3
  x <- simulate_site_data(10000, n_f, 0.5, n_env_vars = n_vars)
  ## set parameters for modelling
  tuning_parameters <- list(max_depth = seq(1, 10, 1),
                            eta = seq(0.1, 0.5, 0.1),
                            nrounds = seq(100, 1000, 100),
                            lambda = 10 ^ seq(-1.0, 0.0, 0.25),
                            subsample = seq(0.5, 1.0, 0.1),
                            colsample_bytree = seq(0.4, 1.0, 0.1),
                            objective = "binary:logistic")
  # main calculations
  r <- tune_occupancy_models(
    x, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
    n_folds = rep(5, n_f), n_random_search_iterations = 10,
    early_stopping_rounds = 5, parameters = tuning_parameters, n_threads = 1)
  # tests
  expect_is(r, "list")
  for (i in seq_along(r)) {
    expect_equal(names(r[[i]]), names(tuning_parameters))
    expect_true(all(is.finite(unlist(r[[i]])) |
                    is.character(unlist(r[[i]]))))
  }
})

test_that("multiple species", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## simulate data
  n_f <- 5
  n_vars <- 3
  x <- simulate_site_data(10000, n_f, 0.5, n_env_vars = n_vars)
  ## set parameters for modelling
  tuning_parameters <- list(max_depth = seq(1, 10, 1),
                            eta = seq(0.1, 0.5, 0.1),
                            nrounds = seq(100, 1000, 100),
                            lambda = 10 ^ seq(-1.0, 0.0, 0.25),
                            subsample = seq(0.5, 1.0, 0.1),
                            colsample_bytree = seq(0.4, 1.0, 0.1),
                            objective = "binary:logistic")
  # main calculations
  r <- tune_occupancy_models(
    x, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
    n_folds = rep(5, n_f), n_random_search_iterations = 10,
    early_stopping_rounds = 5, parameters = tuning_parameters, n_threads = 1)
  # tests
  expect_is(r, "list")
  for (i in seq_along(r)) {
    expect_equal(names(r[[i]]), names(tuning_parameters))
    expect_true(all(is.finite(unlist(r[[i]])) |
                    is.character(unlist(r[[i]]))))
  }
})

test_that("weights species", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## simulate data
  n_f <- 1
  n_vars <- 3
  x <- simulate_site_data(10000, n_f, 0.5, n_env_vars = n_vars)
  x$weight <- runif(nrow(x), 0.1, 1e+4)
  ## set parameters for modelling
  tuning_parameters <- list(max_depth = seq(1, 10, 1),
                            eta = seq(0.1, 0.5, 0.1),
                            nrounds = seq(100, 1000, 100),
                            lambda = 10 ^ seq(-1.0, 0.0, 0.25),
                            subsample = seq(0.5, 1.0, 0.1),
                            colsample_bytree = seq(0.4, 1.0, 0.1),
                            objective = "reg:logistic")
  # main calculations
  set.seed(123)
  r1 <- tune_occupancy_models(
    x, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
    n_folds = rep(5, n_f), n_random_search_iterations = 30,
    early_stopping_rounds = 5, parameters = tuning_parameters, n_threads = 1)
  set.seed(123)
  r2 <- tune_occupancy_models(
    x, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
    n_folds = rep(5, n_f), n_random_search_iterations = 30,
    early_stopping_rounds = 5, parameters = tuning_parameters, n_threads = 1)
  set.seed(123)
  r3 <- tune_occupancy_models(
    x, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
    n_folds = rep(5, n_f), n_random_search_iterations = 30,
    early_stopping_rounds = 5, parameters = tuning_parameters,
    site_weight_columns = "weight", n_threads = 1)
  # tests
  expect_identical(r1, r2)
  expect_false(identical(r1, r3))
})

test_that("species with few presences", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## simulate data
  n_f <- 1
  n_vars <- 3
  x <- simulate_site_data(10000, n_f, 0.5, n_env_vars = n_vars)
  ## set parameters for modelling
  tuning_parameters <- list(max_depth = seq(1, 10, 1),
                            eta = seq(0.1, 0.5, 0.1),
                            nrounds = seq(100, 1000, 100),
                            lambda = 10 ^ seq(-1.0, 0.0, 0.25),
                            subsample = seq(0.5, 1.0, 0.1),
                            colsample_bytree = seq(0.4, 1.0, 0.1),
                            objective = "binary:logistic")
  # manually encoded a single presence
  x$f1[!is.na(x$f1)] <- 0
  x$f1[sample(which(!is.na(x$f1)), 2)] <- 1
  # main calculations
  r <- tune_occupancy_models(
    x, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
    n_folds = rep(5, n_f), n_random_search_iterations = 10,
    early_stopping_rounds = 5, parameters = tuning_parameters, n_threads = 1)
  # tests
  expect_is(r, "list")
  for (i in seq_along(r)) {
    expect_equal(names(r[[i]]), names(tuning_parameters))
    expect_true(all(is.finite(unlist(r[[i]])) |
                    is.character(unlist(r[[i]]))))
  }
})

test_that("species with few absences", {
  # data
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## simulate data
  n_f <- 1
  n_vars <- 3
  x <- simulate_site_data(10000, n_f, 0.5, n_env_vars = n_vars)
  ## set parameters for modelling
  tuning_parameters <- list(max_depth = seq(1, 10, 1),
                            eta = seq(0.1, 0.5, 0.1),
                            nrounds = seq(100, 1000, 100),
                            lambda = 10 ^ seq(-1.0, 0.0, 0.25),
                            subsample = seq(0.5, 1.0, 0.1),
                            colsample_bytree = seq(0.4, 1.0, 0.1),
                            objective = "binary:logistic")
  # manually encoded a single absence
  x$f1[!is.na(x$f1)] <- 1
  x$f1[sample(which(!is.na(x$f1)), 2)] <- 0
  # main calculations
  r <- tune_occupancy_models(
    x, paste0("f", seq_len(n_f)), paste0("e", seq_len(n_vars)),
    n_folds = rep(5, n_f), n_random_search_iterations = 10,
    early_stopping_rounds = 5, parameters = tuning_parameters, n_threads = 1)
  # tests
  expect_is(r, "list")
  for (i in seq_along(r)) {
    expect_equal(names(r[[i]]), names(tuning_parameters))
    expect_true(all(is.finite(unlist(r[[i]])) |
                    is.character(unlist(r[[i]]))))
  }
})
