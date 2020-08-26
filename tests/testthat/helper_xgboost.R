r_xgboost <- function(
  y_train, x_train, w_train, y_test, x_test, w_test, x_predict,
  survey_sensitivity, survey_specificity,
  xgb_parameters, xgb_nrounds, xgb_early_nrounds) {
  # calculate parameters for model fitting
  spw <- round(sum(y_train < 0.5) / sum(y_train > 0.5), 6)
  seed <- as.numeric(xgb_parameters$seed)
  # round data to account for different precision between floats and doubles
  # y_train <- round(y_train, 5)
  # x_train <- round(x_train, 5)
  # w_train <- round(w_train * 100, 5)
  w_train <- w_train * 100
  # y_test <- round(y_test, 5)
  # x_test <- round(x_test, 5)
  # w_test <- round(w_test, 5)
  # prepare data for model fitting
  dtrain <- xgboost::xgb.DMatrix(
    x_train, missing = NA,
    info = list(label = y_train, weight = w_train))
  # prepare data for model evaluation
  dtest <- xgboost::xgb.DMatrix(
    x_test, missing = NA,
    info = list(label = y_test, weight = w_test))
  # prepare evaluation function
  curr_feval_tss <- feval_tss
  environment(curr_feval_tss)$sens <- survey_sensitivity
  environment(curr_feval_tss)$spec <- survey_specificity
  # prepare xgboost call
  args <- list(data = dtrain, verbose = FALSE, scale_pos_weight = spw,
               watchlist = list(test = dtest), eval_metric = curr_feval_tss,
               maximize = TRUE, nrounds = xgb_nrounds,
               early_stopping_rounds = xgb_early_nrounds, nthread = 1)
  args <- append(args, xgb_parameters)
  args$seed <- NULL
  # fit model
  withr::with_seed(seed, {
    model <- do.call(what = xgboost::xgb.train, args)
  })
  # generate predictions
  yhat_train <- predict(model, x_train, ntreelimit = model$best_iteration)
  yhat_test <- predict(model, x_test, ntreelimit = model$best_iteration)
  yhat_predict <- predict(model, x_predict, ntreelimit = model$best_iteration)
  # calculate performance statistics
  perf <- rcpp_model_performance(
    y_test, yhat_test, w_test, survey_sensitivity, survey_specificity)
  # return result
  list(sens = perf[[2]], spec = perf[[3]], tss = perf[[1]], yhat = yhat_predict)
}

r_fit_xgboost_models_and_assess_performance <- function(
  dij, nij, pij, pu_env_data,
  survey_features, survey_sensitivity, survey_specificity,
  tuning_parameters, xgb_nrounds, xgb_early_stopping_rounds,
  xgb_train_folds, xgb_test_folds) {
  ## prepare data
  data <- lapply(which(survey_features), function(i) {
    lapply(seq_along(xgb_train_folds[[i]]), function(k) {
      list(
        y_train = c(rep(1, length(xgb_train_folds[[i]][[k]])),
                    rep(0, length(xgb_train_folds[[i]][[k]]))),
        w_train = c(pij[i, xgb_train_folds[[i]][[k]]],
                    1 - pij[i, xgb_train_folds[[i]][[k]]]),
        x_train = pu_env_data[rep(xgb_train_folds[[i]][[k]], 2), ,
                              drop = FALSE],
        y_test = c(rep(1, length(xgb_test_folds[[i]][[k]])),
                      rep(0, length(xgb_test_folds[[i]][[k]]))),
        w_test = c(dij[i, xgb_test_folds[[i]][[k]]],
                   1 - dij[i, xgb_test_folds[[i]][[k]]]),
        x_test = pu_env_data[rep(xgb_test_folds[[i]][[k]], 2), ,
                              drop = FALSE])
    })
  })

  ## fit models
  models <- lapply(seq_along(which(survey_features)), function(i) {
    lapply(seq_len(nrow(tuning_parameters)), function(p) {
      curr_p <- as.list(tuning_parameters[p, , drop = FALSE])
      names(curr_p) <- colnames(tuning_parameters)
      ii <- which(survey_features)[i]
      lapply(seq_along(xgb_train_folds[[ii]]),
             function(k) {
        r_xgboost(
          y_train = data[[i]][[k]]$y_train,
          x_train = data[[i]][[k]]$x_train,
          w_train = data[[i]][[k]]$w_train,
          y_test = data[[i]][[k]]$y_test,
          x_test = data[[i]][[k]]$x_test,
          w_test = data[[i]][[k]]$w_test,
          x_predict = pu_env_data,
          survey_sensitivity = survey_sensitivity[ii],
          survey_specificity = survey_specificity[ii],
          xgb_parameters = curr_p,
          xgb_nrounds = xgb_nrounds[ii],
          xgb_early_nrounds = xgb_early_stopping_rounds[ii])
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
