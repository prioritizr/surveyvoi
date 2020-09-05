r_predict_missing_rij_data <- function(
  dij, nij, pij, pu_env_data,
  survey_features, survey_sensitivity, survey_specificity,
  tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
  xgb_train_folds, xgb_test_folds, pu_model_prediction) {
  # fit models
  m <- r_fit_xgboost_models_and_assess_performance(
    dij, nij, pij, pu_env_data,
    survey_features, survey_sensitivity, survey_specificity,
    tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
    xgb_train_folds, xgb_test_folds)
  # prepare output
  out <- pij
  survey_features_idx <- which(survey_features)
  for (i in seq_along(survey_features_idx)) {
    out[survey_features_idx[i], pu_model_prediction] <-
      m$pred[pu_model_prediction, i]
  }
  # return result
  attr(out, "dimnames") <- NULL
  out
}
