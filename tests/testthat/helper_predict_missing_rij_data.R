r_predict_missing_rij_data <- function(
  rij, wij, x, survey_features,
  tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
  xgb_train_folds, xgb_test_folds, pu_model_prediction_idx) {
  # fit models
  m <- r_fit_xgboost_models_and_assess_performance(
    rij, wij, x, survey_features,
    tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds)
  # prepare output
  out <- rij
  for (i in seq_along(survey_features)) {
    if (survey_features[i]) {
      out[i, pu_model_prediction_idx[[i]]] <-
        m$pred[pu_model_prediction_idx[[i]], i]
    }
  }
  # return result
  attr(out, "dimnames") <- NULL
  out
}
