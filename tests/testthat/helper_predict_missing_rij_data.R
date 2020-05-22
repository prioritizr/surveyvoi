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
    rij[features[i], pu_model_prediction_idx[[features[i]]]] <-
      rowMeans(yhat)[pu_model_prediction_idx[[features[i]]]]
  }
  # exports
  rij
}
