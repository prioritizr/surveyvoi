weighted_sensitivity <- function(y, yhat, weights) {
  total_positive <- sum(weights * (y >= 0.5))
  true_positive <- sum(weights * ((y >= 0.5) & (yhat >= 0.5)))
  if (total_positive < 1e-10) return(as.numeric(all(yhat == y)))
  true_positive / total_positive
}

weighted_specificity <- function(y, yhat, weights) {
  total_negative <- sum(weights * (y < 0.5))
  true_negative <- sum(weights * ((y < 0.5) & (yhat < 0.5)))
  if (total_negative < 1e-10) return(as.numeric(all(yhat == y)))
  true_negative / total_negative
}

r_xgboost <- function(y, x_train, x_weight, predict_x, xgb_parameters,
                              xgb_nrounds) {
  spw <- round(sum(y < 0.5) / sum(y > 0.5), 6) # C++ code uses 1e-6 precision
  withr::with_seed(as.numeric(xgb_parameters$seed), {
  model <- xgboost::xgboost(
    data = x_train, label = y, weight = round(x_weight, 6),
    nrounds = xgb_nrounds,
    objective = xgb_parameters$objective,
    scale_pos_weight = spw,
    verbose = FALSE)
  })
  predict(model, predict_x)
}

r_fit_xgboost_models_and_assess_performance <- function(
  rij, wij, pu_env_data, survey_features,
  tuning_parameters, xgb_nrounds, xgb_train_folds, xgb_test_folds) {
  ## make predictions
  yhat <- lapply(seq_len(nrow(rij)), function(i) {
    lapply(seq_along(xgb_train_folds[[i]]), function(j) {
      r_xgboost(rij[i, xgb_train_folds[[i]][[j]]],
                pu_env_data[xgb_train_folds[[i]][[j]], , drop = FALSE],
                wij[i, xgb_train_folds[[i]][[j]]],
                pu_env_data,
                tuning_parameters[[i]], xgb_nrounds[i])
    })
  })
  ## calculate sensitivity
  sens <- sapply(seq_len(nrow(rij)), function(i) {
    mean(sapply(seq_along(xgb_train_folds[[i]]), function(j) {
      weighted_sensitivity(rij[i, xgb_test_folds[[i]][[j]]],
                           yhat[[i]][[j]][xgb_test_folds[[i]][[j]]],
                           wij[i, xgb_test_folds[[i]][[j]]])
    }))
  })
  sens <- matrix(sens, ncol = 1)

  ## calculate specificity
  spec <- sapply(seq_len(nrow(rij)), function(i) {
    mean(sapply(seq_along(xgb_train_folds[[i]]), function(j) {
      weighted_specificity(rij[i, xgb_test_folds[[i]][[j]]],
                           yhat[[i]][[j]][xgb_test_folds[[i]][[j]]],
                           wij[i, xgb_test_folds[[i]][[j]]])
    }))
  })
  spec <- matrix(spec, ncol = 1)

  ## clamp values to (1e-10) and (1 - 1e-10) to arise from probabilities
  ## that are exactly zero and one
  spec[] <- pmax(spec[], 1e-10)
  sens[] <- pmax(sens[], 1e-10)
  spec[] <- pmin(spec[], 1 - 1e-10)
  sens[] <- pmin(sens[], 1 - 1e-10)

  ## return result
  list(sens = sens, spec = spec)
}
