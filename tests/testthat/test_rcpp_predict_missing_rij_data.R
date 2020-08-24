context("predict_missing_rij_data")

test_that("single species", {
  # data
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_pu <- 30
  n_f <- 1
  n_vars <- 2
  n_folds <- 2
  ## simulate data
  x <- simulate_site_data(n_pu, n_f, 0.2, n_env_vars = n_vars)
  x <- sf::st_drop_geometry(x)
  y <- simulate_feature_data(n_f)
  survey_features <- y$survey
  survey_sensitivity <- y$survey_sensitivity
  survey_specificity <- y$survey_specificity
  ## create matrices for data
  dij <- t(as.matrix(x[, paste0("f", seq_len(n_f))], ncol = n_f))
  nij <- t(as.matrix(x[, paste0("n", seq_len(n_f))], ncol = n_f))
  pij <- prior_probability_matrix(x, y, paste0("f", seq_len(n_f)),
    paste0("n", seq_len(n_f)), paste0("p", seq_len(n_f)),
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  ## specify environmental data
  pu_env_data <- as.matrix(x[, paste0("e", seq_len(n_vars))])
  ## specify planning units for predictions
  pu_model_prediction_idx <- lapply(seq_len(n_f), function(i) {
    which(nij[i, ] < 0.5)
  })
  ## model fitting parameters
  xgb_folds <- lapply(seq_len(n_f), function(f) {
    fn <- paste0("f", f)
    nn <- paste0("n", f)
    has_data_idx <- which(x[[nn]] > 0)
    create_site_folds(x[[fn]][has_data_idx], x[[nn]][has_data_idx],
                      index = has_data_idx, n_folds)
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  xgb_nrounds <- rep(10, n_f)
  xgb_early_stopping_rounds <- rep(10, n_f)
  tuning_parameters <-
    expand.grid(eta = c(0.1),
                lambda = c(0.001),
                objective = "binary:logistic",
                seed = "123")
  tuning_parameters <- as.matrix(tuning_parameters)
  # results
  r1 <- rcpp_predict_missing_rij_data(
    dij, nij, pij, pu_env_data,
    survey_features, survey_sensitivity, survey_specificity,
    colnames(tuning_parameters), tuning_parameters,
    xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds, pu_model_prediction_idx)
  r2 <- r_predict_missing_rij_data(
    dij, nij, pij, pu_env_data,
    survey_features, survey_sensitivity, survey_specificity,
    tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds, pu_model_prediction_idx)
  ## tests
  ## expect rcpp and r implementations to give same result
  expect_equal(r1, r2)
  ## expect values for pu's that shouldn't be updated to remain the same
  for (i in seq_len(n_f)) {
    expect_equal(
      r1[i, -1 * pu_model_prediction_idx[[i]]],
      pij[i, -1 * pu_model_prediction_idx[[i]]])
  }
  ## expect values for species that shouldn't be updated to remain the same
  for (i in seq_len(n_f)) {
    if (!survey_features[i]) {
      expect_equal(r1[i, ], pij[i, ])
    }
  }
  ## expect model predictions to be semi-accurate
  for (i in seq_len(n_f)) {
    if (!survey_features[i]) {
      expect_lte(max(abs(
        r1[i, pu_model_prediction_idx[[i]]] -
        x[[paste0("p", i)]][pu_model_prediction_idx[[i]]])), 0.5)
    }
  }
})

test_that("multiple species", {
  # data
  ## set seeds
  set.seed(123)
  RandomFields::RFoptions(seed = 123)
  ## set constants
  n_pu <- 30
  n_f <- 3
  n_vars <- 2
  n_folds <- 2
  ## simulate data
  x <- simulate_site_data(n_pu, n_f, 0.2, n_env_vars = n_vars)
  x <- sf::st_drop_geometry(x)
  y <- simulate_feature_data(n_f)
  y$survey[2] <- FALSE
  survey_features <- y$survey
  survey_sensitivity <- y$survey_sensitivity
  survey_specificity <- y$survey_specificity
  ## create matrices for data
  dij <- t(as.matrix(x[, paste0("f", seq_len(n_f))], ncol = n_f))
  nij <- t(as.matrix(x[, paste0("n", seq_len(n_f))], ncol = n_f))
  pij <- prior_probability_matrix(x, y, paste0("f", seq_len(n_f)),
    paste0("n", seq_len(n_f)), paste0("p", seq_len(n_f)),
    "survey_sensitivity", "survey_specificity",
    "model_sensitivity", "model_specificity")
  ## specify environmental data
  pu_env_data <- as.matrix(x[, paste0("e", seq_len(n_vars))])
  ## specify planning units for predictions
  pu_model_prediction_idx <- lapply(seq_len(n_f), function(i) {
    which(nij[i, ] < 0.5)
  })
  ## model fitting parameters
  xgb_folds <- lapply(seq_len(n_f), function(f) {
    fn <- paste0("f", f)
    nn <- paste0("n", f)
    has_data_idx <- which(x[[nn]] > 0)
    create_site_folds(x[[fn]][has_data_idx], x[[nn]][has_data_idx],
                      index = has_data_idx, n_folds)
  })
  xgb_train_folds <- lapply(xgb_folds, `[[`, "train")
  xgb_test_folds <- lapply(xgb_folds, `[[`, "test")
  xgb_nrounds <- rep(10, n_f)
  xgb_early_stopping_rounds <- rep(10, n_f)
  tuning_parameters <-
    expand.grid(eta = c(0.1),
                lambda = c(0.001),
                objective = "binary:logistic",
                seed = "123")
  tuning_parameters <- as.matrix(tuning_parameters)
  # results
  r1 <- rcpp_predict_missing_rij_data(
    dij, nij, pij, pu_env_data,
    survey_features, survey_sensitivity, survey_specificity,
    colnames(tuning_parameters), tuning_parameters,
    xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds, pu_model_prediction_idx)
  r2 <- r_predict_missing_rij_data(
    dij, nij, pij, pu_env_data,
    survey_features, survey_sensitivity, survey_specificity,
    tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds, pu_model_prediction_idx)
  ## tests
  ## expect rcpp and r implementations to give same result
  expect_lte(abs(max(r1 - r2)), 1e-7)
  ## expect values for pu's that shouldn't be updated to remain the same
  for (i in seq_len(n_f)) {
    expect_equal(
      r1[i, -1 * pu_model_prediction_idx[[i]]],
      pij[i, -1 * pu_model_prediction_idx[[i]]])
  }
  ## expect values for species that shouldn't be updated to remain the same
  for (i in seq_len(n_f)) {
    if (!survey_features[i]) {
      expect_equal(r1[i, ], pij[i, ])
    }
  }
  ## expect model predictions to be semi-accurate
  for (i in seq_len(n_f)) {
    if (!survey_features[i]) {
      expect_lte(max(abs(
        r1[i, pu_model_prediction_idx[[i]]] -
        x[[paste0("p", i)]][pu_model_prediction_idx[[i]]])), 0.5)
    }
  }
})
