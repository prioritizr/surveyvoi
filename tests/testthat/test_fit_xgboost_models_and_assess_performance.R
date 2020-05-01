context("fit_xgboost_models_and_assess_performance")

test_that("equal weights", {
  # data
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_pu <- 1000
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
    na_idx <- which(is.na(x[[f]]))
    o <- create_folds(x[[f]], n_folds, na.fail = FALSE)
    o$train <- lapply(o$train, function(z) z[!z %in% na_idx])
    o$test <- lapply(o$test, function(z) z[!z %in% na_idx])
    o
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  xgb_nrounds <- rep(10, n_f)
  tuning_parameters <-
    list(list(objective = "binary:logistic", scale_pos_weight = "2",
              seed = "123"))[rep(1, n_f)]
  # run calculations
  r1 <- rcpp_fit_xgboost_models_and_assess_performance(
    rij, wij, pu_env_data, survey_features, tuning_parameters, xgb_nrounds,
    xgb_train_folds, xgb_test_folds)
  r2 <- r_fit_xgboost_models_and_assess_performance(
    rij, wij, pu_env_data, survey_features, tuning_parameters, xgb_nrounds,
    xgb_train_folds, xgb_test_folds)
  # tests
  expect_equal(r1$sens, r2$sens)
  expect_equal(r1$spec, r2$spec)
})

test_that("variable weights", {
  # data
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_pu <- 10
  n_f <- 3
  n_vars <- 3
  n_folds <- 5
  ## simulate data
  x <- simulate_site_data(n_pu, n_f, 0.2, n_env_vars = n_vars)
  x <- sf::st_drop_geometry(x)
  survey_features <- rep(TRUE, n_f)
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
    na_idx <- which(is.na(x[[f]]))
    o <- create_folds(x[[f]], n_folds, na.fail = FALSE)
    o$train <- lapply(o$train, function(z) z[!z %in% na_idx])
    o$test <- lapply(o$test, function(z) z[!z %in% na_idx])
    o
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  xgb_nrounds <- rep(10, n_f)
  tuning_parameters <-
    list(list(objective = "binary:logistic", scale_pos_weight = "2",
              seed = "123"))[rep(1, n_f)]
  # run calculations
  r1 <- rcpp_fit_xgboost_models_and_assess_performance(
    rij, wij, pu_env_data, survey_features, tuning_parameters, xgb_nrounds,
    xgb_train_folds, xgb_test_folds)
  r2 <- r_fit_xgboost_models_and_assess_performance(
    rij, wij, pu_env_data, survey_features, tuning_parameters, xgb_nrounds,
    xgb_train_folds, xgb_test_folds)
  # tests
  expect_lte(max(abs(r1$sens - r2$sens)), 1e-7)
  expect_lte(max(abs(r1$spec - r2$spec)), 1e-7)
})
