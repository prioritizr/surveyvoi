r_xgboost <- function(y, x_train, x_weights, predict_x, xgb_parameters,
                              xgb_nrounds) {
  spw <- round(sum(y < 0.5) / sum(y > 0.5), 6) # C++ code uses 1e-6 precision
  set.seed(as.numeric(xgb_parameters$seed))
  model <- xgboost::xgboost(
    data = x_train, label = y, nrounds = xgb_nrounds, weight = x_weights,
    objective = xgb_parameters$objective,
    scale_pos_weight = spw,
    verbose = FALSE)
  set.seed(as.numeric(xgb_parameters$seed))
  predict(model, predict_x)
}

r_predict_missing_rij_data <- function(
  rij, wij, x, features, pu_model_prediction_idx,
  xgb_parameters, xgb_nrounds, xgb_train_folds, xgb_test_folds) {
  # main
  for (i in seq_along(features)) {
    # generate predictions
    yhat <- sapply(seq_along(xgb_train_folds[[i]]), function(j) {
      r_xgboost(
        rij[features[i], xgb_train_folds[[i]][[j]]],
        x[xgb_train_folds[[i]][[j]], , drop = FALSE],
        wij[features[i], xgb_train_folds[[i]][[j]]],
        x, xgb_parameters[[features[i]]], xgb_nrounds[features[i]])
    })
    # assign values
    rij[features[i], pu_model_prediction_idx] <-
      rowMeans(yhat)[pu_model_prediction_idx]
  }
  # exports
  rij
}
