context("fit_xgboost_models_and_assess_performance")

test_that("equal weights", {
  # data
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_pu <- 100
  n_f <- 1
  n_vars <- 2
  n_folds <- 5
  ## simulate data
  x <- simulate_site_data(n_pu, n_f, 0.2, n_env_vars = n_vars)
  x <- sf::st_drop_geometry(x)
  survey_features <- rep(TRUE, n_f)
  # weights
  wij <- matrix(1, ncol = n_pu, nrow = n_f)
  ## specify presence/absences
  rij <- t(as.matrix(x[, paste0("f", seq_len(n_f))], ncol = n_f))
  rij[is.na(rij)] <- -1
  ## specify environmental data
  pu_env_data <- as.matrix(x[, paste0("e", seq_len(n_vars))])
  ## model fitting parameters
  xgb_folds <- lapply(paste0("f", seq_len(n_f)), function(f) {
    non_na_idx <- which(!is.na(x[[f]]))
    create_folds(x[[f]][non_na_idx], index = non_na_idx, n_folds)
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  xgb_nrounds <- rep(10, n_f)
  xgb_early_stopping_rounds <- rep(5, n_f)
  tuning_parameters <-
    expand.grid(eta = c(0.1, 0.5, 1.0),
                lambda = c(0.001, 0.01, 0.05),
                objective = "binary:logistic",
                seed = "123")
  tuning_parameters <- as.matrix(tuning_parameters)
  # run calculations
  r1 <- rcpp_fit_xgboost_models_and_assess_performance(
    rij, wij, pu_env_data, survey_features,
    colnames(tuning_parameters), tuning_parameters,
    xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds)
  r2 <- r_fit_xgboost_models_and_assess_performance(
    rij, wij, pu_env_data, survey_features,
    tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds)
  # tests
  expect_equal(r1$sens, r2$sens)
  expect_equal(r1$spec, r2$spec)
})

test_that("variable weights", {
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_pu <- 100
  n_f <- 3
  n_vars <- 2
  n_folds <- 5
  ## simulate data
  x <- simulate_site_data(n_pu, n_f, 0.2, n_env_vars = n_vars)
  x <- sf::st_drop_geometry(x)
  survey_features <- c(TRUE, FALSE, rep(TRUE, n_f - 2))
  ## specify presence/absences
  rij <- t(as.matrix(x[, paste0("f", seq_len(n_f))], ncol = n_f))
  rij[is.na(rij)] <- -1
  # weights
  wij <- matrix(round(runif(length(rij)), 6),
                ncol = ncol(rij), nrow = nrow(rij))
  ## specify environmental data
  pu_env_data <- as.matrix(x[, paste0("e", seq_len(n_vars))])
  ## model fitting parameters
  xgb_folds <- lapply(paste0("f", seq_len(n_f)), function(f) {
    non_na_idx <- which(!is.na(x[[f]]))
    create_folds(x[[f]][non_na_idx], index = non_na_idx, n_folds)
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  xgb_nrounds <- rep(100, n_f)
  xgb_early_stopping_rounds <- rep(20, n_f)
  tuning_parameters <-
    expand.grid(eta = c(0.1, 0.5, 1.0),
                lambda = c(0.001, 0.01, 0.05),
                objective = "binary:logistic",
                seed = "123")
  tuning_parameters <- as.matrix(tuning_parameters)
  # run calculations
  r1 <- rcpp_fit_xgboost_models_and_assess_performance(
    rij, wij, pu_env_data, survey_features,
    colnames(tuning_parameters), tuning_parameters,
    xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds)
  r2 <- r_fit_xgboost_models_and_assess_performance(
    rij, wij, pu_env_data, survey_features,
    tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds)
  # tests
  expect_lte(max(abs(r1$sens - r2$sens)), 1e-6)
  expect_lte(max(abs(r1$spec - r2$spec)), 1e-6)
})
