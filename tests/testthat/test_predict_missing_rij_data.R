context("predict_missing_rij_data")

test_that("equal weights", {
  # data
  ## set seed
  set.seed(401)
  ## set constants for simulating data
  n_total_f <- 5
  n_f <- ceiling(n_total_f * 0.5)
  n_pu <- 20
  n_vars <- 3
  n_folds <- 5
  ## simulate planning units that need model predictions
  pu_model_prediction_idx <- sample.int(n_pu, ceiling(n_pu * 0.5))
  pu_model_prediction_idx <- sort(pu_model_prediction_idx)
  ## simulate planning units that have been surveyed
  pu_survey_solution_idx <- setdiff(seq_len(n_pu), pu_model_prediction_idx)
  pu_survey_solution_idx <- sample(pu_survey_solution_idx,
                                   ceiling(length(pu_survey_solution_idx) *
                                           0.4))
  pu_survey_solution_idx <- sort(pu_survey_solution_idx)
  ## simulate features that have uncertain distributions
  features <- sort(sample.int(n_total_f, n_f))
  survey_features <- replace(rep(FALSE, n_total_f), features, TRUE)
  survey_features_idx <- which(replace(rep(FALSE, n_total_f), features, TRUE))
  ## simulate matrix of predictor variables for modelling,
  ## including an intercept term
  x <- matrix(rnorm(n_vars * n_pu), nrow = n_pu)
  x[, 1] <- 1
  ## simulate coefficients
  beta <- matrix(rnorm(n_vars * n_total_f), ncol = n_vars, nrow = n_total_f)
  ## simulate rij data
  rij <- matrix(NA, ncol = n_pu, nrow = n_total_f)
  for (i in seq_len(n_total_f)) {
    rij[i, ] <- (x %*% beta[i, ])
    rij[i, ] <- exp(rij[i, ]) / (1 + exp(rij[i, ]))
  }
  rij[] <- rbinom(prod(dim(rij)), 1, rij[])
  true_rij <- rij
  rij[, pu_model_prediction_idx] <- -1
  ## set constant weights
  wij <- rij
  wij[] <- 1
  ## set xgboost modelling information
  xgb_parameters <-
    list(list(seed = "0", scale_pos_weight = "2",
              objective = "binary:logistic"))[rep(1, n_total_f)]
  xgb_nrounds <- rep(10, n_total_f)
  xgb_folds <- lapply(survey_features_idx, function(i) {
    na_idx <- which(rij[i, ] < -0.5)
    o <- create_folds(rij[i, ], n_folds, na.fail = FALSE)
    o$train <- lapply(o$train, function(z) z[!z %in% na_idx])
    o$test <- lapply(o$test, function(z) z[!z %in% na_idx])
    o
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  # results
  r1 <- r_predict_missing_rij_data(
    rij, wij, x, features, pu_model_prediction_idx,
    xgb_parameters, xgb_nrounds, xgb_train_folds, xgb_test_folds)
  r2 <- rcpp_predict_missing_rij_data(
    rij, wij, x, survey_features,
    pu_model_prediction_idx - 1, xgb_parameters, xgb_nrounds, xgb_train_folds,
    xgb_test_folds)
  # test
  expect_lte(max(abs(r1 - r2)), 1e-7)
})

test_that("variable weights", {
  # data
  ## set seed
  set.seed(401)
  ## set constants for simulating data
  n_total_f <- 5
  n_f <- ceiling(n_total_f * 0.5)
  n_pu <- 20
  n_vars <- 3
  n_folds <- 5
  ## simulate planning units that need model predictions
  pu_model_prediction_idx <- sample.int(n_pu, ceiling(n_pu * 0.5))
  pu_model_prediction_idx <- sort(pu_model_prediction_idx)
  ## simulate planning units that have been surveyed
  pu_survey_solution_idx <- setdiff(seq_len(n_pu), pu_model_prediction_idx)
  pu_survey_solution_idx <- sample(pu_survey_solution_idx,
                                   ceiling(length(pu_survey_solution_idx) *
                                           0.4))
  pu_survey_solution_idx <- sort(pu_survey_solution_idx)
  ## simulate features that have uncertain distributions
  features <- sort(sample.int(n_total_f, n_f))
  survey_features <- replace(rep(FALSE, n_total_f), features, TRUE)
  survey_features_idx <- which(replace(rep(FALSE, n_total_f), features, TRUE))
  ## simulate matrix of predictor variables for modelling,
  ## including an intercept term
  x <- matrix(rnorm(n_vars * n_pu), nrow = n_pu)
  x[, 1] <- 1
  ## simulate coefficients
  beta <- matrix(rnorm(n_vars * n_total_f), ncol = n_vars, nrow = n_total_f)
  ## simulate rij data
  rij <- matrix(NA, ncol = n_pu, nrow = n_total_f)
  for (i in seq_len(n_total_f)) {
    rij[i, ] <- (x %*% beta[i, ])
    rij[i, ] <- exp(rij[i, ]) / (1 + exp(rij[i, ]))
  }
  rij[] <- rbinom(prod(dim(rij)), 1, rij[])
  true_rij <- rij
  rij[, pu_model_prediction_idx] <- -1
  ## set constant weights
  wij <- rij
  wij[] <- round(runif(length(wij)), 6)
  ## set xgboost modelling information
  xgb_parameters <-
    list(list(seed = "0", scale_pos_weight = "2",
              objective = "binary:logistic"))[rep(1, n_total_f)]
  xgb_nrounds <- rep(10, n_total_f)
  xgb_folds <- lapply(survey_features_idx, function(i) {
    na_idx <- which(rij[i, ] < -0.5)
    o <- create_folds(rij[i, ], n_folds, na.fail = FALSE)
    o$train <- lapply(o$train, function(z) z[!z %in% na_idx])
    o$test <- lapply(o$test, function(z) z[!z %in% na_idx])
    o
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  # results
  r1 <- r_predict_missing_rij_data(
    rij, wij, x, features, pu_model_prediction_idx,
    xgb_parameters, xgb_nrounds, xgb_train_folds, xgb_test_folds)
  r2 <- rcpp_predict_missing_rij_data(
    rij, wij, x, survey_features,
    pu_model_prediction_idx - 1, xgb_parameters, xgb_nrounds, xgb_train_folds,
    xgb_test_folds)
  # test
  expect_lte(max(abs(r1 - r2)), 1e-7)
})
