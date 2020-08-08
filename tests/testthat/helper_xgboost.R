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

weighted_tss <- function(y, yhat, weights) {
  weighted_sensitivity(y, yhat, weights) +
  weighted_specificity(y, yhat, weights) - 1
}

r_xgboost <- function(
  y_train, x_train, w_train, y_test, x_test, w_test, x_predict,
  xgb_parameters, xgb_nrounds, xgb_early_nrounds) {
  # calculate parameters for model fitting
  spw <- round(sum(y_train < 0.5) / sum(y_train > 0.5), 6)
  seed <- as.numeric(xgb_parameters$seed)
  # prepare data for model fitting
  dtrain <- xgboost::xgb.DMatrix(
    x_train, missing = NA,
    info = list(label = y_train, weight =round(w_train, 6)))
  # prepare data for model evaluation
  dtest <- xgboost::xgb.DMatrix(
    x_test, missing = NA,
    info = list(label = y_test, weight = round(w_test, 6)))
  # prepare xgboost call
  args <- list(data = dtrain, verbose = FALSE, scale_pos_weight = spw,
               watchlist = list(test = dtest), eval_metric = feval_tss,
               maximize = TRUE, nrounds = xgb_nrounds,
               early_stopping_rounds = xgb_early_nrounds)
  args <- append(args, xgb_parameters)
  args$seed <- NULL
  # fit model
  withr::with_seed(seed, {
    model <- do.call(what = xgboost::xgb.train, args)
  })
  # calculate performance statistics
  yhat_test <- predict(model, x_test, ntreelimit = model$best_ntreelimit)
  sens <- weighted_sensitivity(y_test, yhat_test, w_test)
  spec <- weighted_specificity(y_test, yhat_test, w_test)
  tss <- weighted_tss(y_test, yhat_test, w_test)
  # generate predictions
  yhat_predict <- predict(model, x_predict, ntreelimit = model$best_ntreelimit)
  # return result
  list(sens = sens, spec = spec, tss = tss, yhat = yhat_predict)
}

r_fit_xgboost_models_and_assess_performance <- function(
  rij, wij, pu_env_data, survey_features,
  tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
  xgb_train_folds, xgb_test_folds) {
  ## fit models
  models <- lapply(which(survey_features), function(i) {
    lapply(seq_len(nrow(tuning_parameters)), function(p) {
      curr_p <- as.list(tuning_parameters[p, , drop = FALSE])
      names(curr_p) <- colnames(tuning_parameters)
      lapply(seq_along(xgb_train_folds[[i]]), function(k) {
        r_xgboost(
          y_train = rij[i, xgb_train_folds[[i]][[k]]],
          x_train = pu_env_data[xgb_train_folds[[i]][[k]], , drop = FALSE],
          w_train = wij[i, xgb_train_folds[[i]][[k]]],
          y_test = rij[i, xgb_test_folds[[i]][[k]]],
          x_test = pu_env_data[xgb_test_folds[[i]][[k]], , drop = FALSE],
          w_test = wij[i, xgb_test_folds[[i]][[k]]],
          x_predict = pu_env_data,
          xgb_parameters = curr_p,
          xgb_nrounds = xgb_nrounds[i],
          xgb_early_nrounds =  xgb_early_stopping_rounds[i])
      })
    })
  })
  ## calculate average tss per feature
  tss <- lapply(models, function(x) {
    colSums(sapply(x, function(z) {sapply(z, function(w) w[["tss"]])}))
  })
  ## find the best model per species
  best_model_idx <- vapply(tss, which.max, numeric(1))
  ## calculate average sensitivity of best model per feature
  sens <- sapply(seq_along(best_model_idx), function(i) {
    mean(vapply(models[[i]][[best_model_idx[i]]], `[[`, numeric(1), "sens"))
  })
  ## calculate average specificity of best model per feature
  spec <- sapply(seq_along(best_model_idx), function(i) {
    mean(vapply(models[[i]][[best_model_idx[i]]], `[[`, numeric(1), "spec"))
  })
  ## calculate average probability predictions
  pred <- sapply(seq_along(best_model_idx), function(i) {
    rowMeans(vapply(models[[i]][[best_model_idx[i]]], `[[`,
             numeric(nrow(pu_env_data)), "yhat"))
  })
  ## clamp values to (1e-10) and (1 - 1e-10) to arise from probabilities
  ## that are exactly zero and one
  spec <- pmax(spec, 1e-10)
  sens <- pmax(sens, 1e-10)
  spec <- pmin(spec, 1 - 1e-10)
  sens <- pmin(sens, 1 - 1e-10)
  ## return result
  list(sens = sens, spec = spec, pred = pred)
}
