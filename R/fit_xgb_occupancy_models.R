#' @include internal.R
NULL

#' Fit boosted regression tree models to predict occupancy
#'
#' Estimate probability of occupancy for a set of features in a set of
#' planning units. Models are fitted using gradient boosted trees (via
#' [xgboost::xgb.train()]).
#'
#' @param site_data [sf::sf()] object with site data.
#'
#' @param feature_data [base::data.frame()] object with feature data.
#'
#' @param site_detection_columns `character` names of `numeric`
#'   columns in the argument to `site_data` that contain the proportion of
#'   surveys conducted within each site that detected each feature.
#'   Each column should correspond to a different feature, and contain
#'   a proportion value (between zero and one). If a site has
#'   not previously been surveyed, a value of zero should be used.
#'
#' @param site_n_surveys_columns `character` names of `numeric`
#'   columns in the argument to `site_data` that contain the total
#'   number of surveys conducted for each each feature within each site.
#'   Each column should correspond to a different feature, and contain
#'   a non-negative integer number (e.g. 0, 1, 2, 3). If a site has
#'   not previously been surveyed, a value of zero should be used.
#'
#' @param site_env_vars_columns `character` names of columns in the
#'   argument to `site_data` that contain environmental information
#'   for fitting updated occupancy models based on possible survey outcomes.
#'   Each column should correspond to a different environmental variable,
#'   and contain `numeric`, `factor`, or `character` data.
#'   No missing (`NA`) values are permitted in these columns.
#'
#' @param feature_survey_sensitivity_column `character` name of the
#'   column in the argument to `feature_data` that contains
#'   probability of future surveys correctly detecting a presence of each
#'   feature in a given site (i.e. the sensitivity of the survey methodology).
#'   This column should have `numeric` values that are between zero and
#'   one. No missing (`NA`) values are permitted in this column.
#'
#' @param feature_survey_specificity_column `character` name of the
#'   column in the argument to `feature_data` that contains
#'   probability of future surveys correctly detecting an absence of each
#'   feature in a given site (i.e. the specificity of the survey methodology).
#'   This column should have `numeric` values that are between zero and
#'   one. No missing (`NA`) values are permitted in this column.
#'
#' @param xgb_tuning_parameters `list` object containing the candidate
#'  parameter values for fitting models. Valid parameters include:
#'  `"max_depth"`, `"eta"`, `"lambda"`,
#'  `"min_child_weight"`, `"subsample"`, `"colsample_by_tree"`,
#'  `"objective"`. See documentation for the `params` argument in
#'  [xgboost::xgb.train()] for more information.
#'
#' @param xgb_early_stopping_rounds `numeric` model rounds for parameter
#'   tuning. See [xgboost::xgboost()] for more information.
#'   Defaults to 10 for each feature.
#'
#' @param xgb_n_rounds `numeric` model rounds for model fitting
#'   See [xgboost::xgboost()] for more information.
#'   Defaults to 100 for each feature.
#'
#' @param n_folds `numeric` number of folds to split the training
#'   data into when fitting models for each feature.
#'   Defaults to 5 for each feature.
#'
#' @param n_threads `integer` number of threads to use for parameter
#'   tuning. Defaults to 1.
#'
#' @param seed `integer` initial random number generator state for model
#'   fitting. Defaults to 500.
#'
#' @param verbose `logical` indicating if information should be
#'   printed during computations. Defaults to `FALSE`.
#'
#' @details
#'  This function (i) prepares the data for model fitting, (ii) calibrates
#'  the tuning parameters for model fitting, (iii) generate predictions using
#'  the best found tuning parameters, and (iv) assess the performance of the
#'  best supported models. These analyses are performed separately for each
#'  feature. For a given feature:
#'
#'  \enumerate{
#'
#'  \item The data are prepared for model fitting by partitioning the data using
#'  k-fold cross-validation (set via argument to `n_folds`). The
#'  training and evaluation folds are constructed
#'  in such a manner as to ensure that each training and evaluation
#'  fold contains at least one presence and one absence observation.
#'
#'  \item A grid search method is used to tune the model parameters. The
#'  candidate values for each parameter (specified via `parameters`) are
#'  used to generate a full set of parameter combinations, and these
#'  parameter combinations are subsequently used for tuning the models.
#'  To account for unbalanced datasets, the
#'  `scale_pos_weight` [xgboost::xgboost()] parameter
#'  is calculated as the mean value across each of the training folds
#'  (i.e. number of absence divided by number of presences per feature).
#'  For a given parameter combination, models are fit using k-fold cross-
#'  validation (via [xgboost::xgb.cv()]) -- using the previously
#'  mentioned training and evaluation folds -- and the True Skill Statistic
#'  (TSS) calculated using the data held out from each fold is
#'  used to quantify the performance (i.e. `"test_tss_mean"` column in
#'  output). These models are also fitted using the
#'  `early_stopping_rounds` parameter to reduce time-spent
#'  tuning models. If relevant, they are also fitted using the supplied weights
#'  (per by the argument to `site_weights_data`). After exploring the
#'  full set of parameter combinations, the best parameter combination is
#'  identified, and the associated parameter values and models are stored for
#'  later use.
#'
#'  \item The cross-validation models associated with the best parameter
#'   combination are used to generate predict the average probability that the
#'   feature occupies each site. These predictions include sites that have
#'   been surveyed before, and also sites that have not been surveyed before.
#'
#'  \item The performance of the cross-validation models is evaluated.
#'  Specifically, the TSS, sensitivity, and specificity statistics are
#'  calculated (if relevant, weighted by the argument to
#'  `site_weights_data`). These performance values are calculated using
#'  the models' training and evaluation folds.
#'
#' }
#'
#' @return `list` object containing:
#' \describe{
#'
#' \item{parameters}{`list` of `list` objects containing the best
#' tuning parameters for each feature.}
#'
#' \item{predictions}{[tibble::tibble()] object containing
#'  predictions for each feature.}
#'
#' \item{performance}{[tibble::tibble()] object containing the
#'  performance of the best models for each feature. It contains the following
#'  columns:
#'
#'  \describe{
#'  \item{feature}{name of the feature.}
#'  \item{train_tss_mean}{
#'    mean TSS statistic for models calculated using training data in
#'    cross-validation.}
#'  \item{train_tss_std}{
#'    standard deviation in TSS statistics for models calculated using training
#'    data in cross-validation.}
#'  \item{train_sensitivity_mean}{
#'    mean sensitivity statistic for models calculated using training data in
#'    cross-validation.}
#'  \item{train_sensitivity_std}{
#'    standard deviation in sensitivity statistics for models calculated using
#'    training data in cross-validation.}
#'  \item{train_specificity_mean}{
#'    mean specificity statistic for models calculated using training data in
#'    cross-validation.}
#'  \item{train_specificity_std}{
#'    standard deviation in specificity statistics for models calculated using
#'    training data in cross-validation.}
#'  \item{test_tss_mean}{
#'    mean TSS statistic for models calculated using test data in
#'    cross-validation.}
#'  \item{test_tss_std}{
#'    standard deviation in TSS statistics for models calculated using test
#'    data in cross-validation.}
#'  \item{test_sensitivity_mean}{
#'    mean sensitivity statistic for models calculated using test data in
#'    cross-validation.}
#'  \item{test_sensitivity_std}{
#'    standard deviation in sensitivity statistics for models calculated using
#'    test data in cross-validation.}
#'  \item{test_specificity_mean}{
#'    mean specificity statistic for models calculated using test data in
#'    cross-validation.}
#'  \item{test_specificity_std}{
#'    standard deviation in specificity statistics for models calculated using
#'    test data in cross-validation.}
#'  }}
#'
#' }
#'
#' @examples
#' # set seeds for reproducibility
#' library(RandomFields)
#' set.seed(123)
#' RFoptions(seed = 123)
#'
#' # simulate data for 200 sites, 2 features, and 3 environmental variables
#' site_data <- simulate_site_data(n_sites = 30, n_features = 2, prop = 0.1)
#' feature_data <- simulate_feature_data(n_features = 2, prop = 1)
#'
#' # create list of possible tuning parameters for modelling
#' parameters <- list(eta = seq(0.1, 0.5, length.out = 3),
#'                    lambda = 10 ^ seq(-1.0, 0.0, length.out = 3),
#'                    objective = "binary:logistic")
#'
#' \dontrun{
#' # fit models
#' # note that we use 10 random search iterations here so that the example
#' # finishes quickly, you would probably want something like 1000+
#' results <- fit_xgb_occupancy_models(
#'    site_data, feature_data,
#'    c("f1", "f2"), c("n1", "n2"), c("e1", "e2", "e3"),
#'    "survey_sensitivity", "survey_specificity",
#'    n_folds = rep(5, 2), xgb_early_stopping_rounds = rep(100, 2),
#'    xgb_tuning_parameters = parameters, n_threads = 1)
#'
#' # print best found model parameters
#' print(results$parameters)
#'
#' # print model predictions
#' print(results$predictions)
#'
#' # print model performance
#' print(results$performance, width = Inf)
#' }
#' @export
fit_xgb_occupancy_models <- function(
  site_data, feature_data,
  site_detection_columns, site_n_surveys_columns,
  site_env_vars_columns,
  feature_survey_sensitivity_column, feature_survey_specificity_column,
  xgb_tuning_parameters,
  xgb_early_stopping_rounds = rep(20, length(site_detection_columns)),
  xgb_n_rounds = rep(100, length(site_detection_columns)),
  n_folds = rep(5, length(site_detection_columns)),
  n_threads = 1, seed = 500, verbose = FALSE) {
  # assert that arguments are valid
  assertthat::assert_that(
    ## site data
    inherits(site_data, "sf"), nrow(site_data) > 0, ncol(site_data) > 0,
    ## feature data
    inherits(feature_data, "data.frame"), nrow(feature_data) > 0,
    ncol(feature_data) > 0,
    ## site_detection_columns
    is.character(site_detection_columns),
    length(site_detection_columns) > 0,
    assertthat::noNA(site_detection_columns),
    all(assertthat::has_name(site_data, site_detection_columns)),
    length(site_detection_columns) == nrow(feature_data),
    ## site_n_surveys_columns
    is.character(site_n_surveys_columns),
    length(site_n_surveys_columns) > 0,
    assertthat::noNA(site_n_surveys_columns),
    all(assertthat::has_name(site_data, site_n_surveys_columns)),
    length(site_n_surveys_columns) == nrow(feature_data),
    ## site_env_vars_columns
    is.character(site_env_vars_columns),
    assertthat::noNA(site_env_vars_columns),
    all(assertthat::has_name(site_data, site_env_vars_columns)),
    ## feature_survey_sensitivity_column
    assertthat::is.string(feature_survey_sensitivity_column),
    assertthat::noNA(feature_survey_sensitivity_column),
    all(assertthat::has_name(feature_data, feature_survey_sensitivity_column)),
    is.numeric(feature_data[[feature_survey_sensitivity_column]]),
    ## feature_survey_specificity_column
    assertthat::is.string(feature_survey_specificity_column),
    assertthat::noNA(feature_survey_specificity_column),
    all(assertthat::has_name(feature_data, feature_survey_specificity_column)),
    is.numeric(feature_data[[feature_survey_specificity_column]]),
    ## xgb_tuning_parameters
    is.list(xgb_tuning_parameters),
    ## n_folds
    is.numeric(n_folds),
    assertthat::noNA(n_folds),
    all(xgb_n_rounds > 0),
    length(n_folds) == nrow(feature_data),
    ## xgb_n_rounds
    is.numeric(xgb_n_rounds),
    assertthat::noNA(xgb_n_rounds),
    length(xgb_n_rounds) == nrow(feature_data),
    all(xgb_n_rounds > 0),
    ## xgb_early_stopping_rounds
    is.numeric(xgb_early_stopping_rounds),
    assertthat::noNA(xgb_early_stopping_rounds),
    all(xgb_early_stopping_rounds > 0),
    length(xgb_early_stopping_rounds) == nrow(feature_data),
    ## seed
    assertthat::is.count(seed),
    assertthat::noNA(seed),
    ## n_threads
    assertthat::is.count(n_threads), assertthat::noNA(n_threads))
  ## validate survey data
  validate_site_detection_data(site_data, site_detection_columns)
  validate_site_n_surveys_data(site_data, site_n_surveys_columns)
  ## validate tuning parameters
  validate_xgboost_tuning_parameters(xgb_tuning_parameters)

  # drop geometry
  site_data <- sf::st_drop_geometry(site_data)

  # convert env data to model matrix format
  site_env_data <-
    as.matrix(stats::model.matrix(~ . - 1,
                                  data = site_data[, site_env_vars_columns]))

  # prepare folds
  f <- lapply(seq_len(nrow(feature_data)), function(i) {
    n_surveys <- site_data[[site_n_surveys_columns[i]]]
    idx <- seq_along(n_surveys)
    create_site_folds(
      prop_detected = site_data[[site_detection_columns[i]]][n_surveys > 0],
      n_total = n_surveys[n_surveys > 0],
      n = n_folds[i], idx[n_surveys > 0], seed = seed)
  })

  # prepare data
  d <- lapply(seq_len(nrow(feature_data)), function(i) {
    lapply(seq_len(n_folds[i]), function(k) {
      # extract data
      sens <- feature_data[[feature_survey_sensitivity_column]][i]
      spec <- feature_data[[feature_survey_specificity_column]][i]
      n_surveys <- site_data[[site_n_surveys_columns[i]]]
      n_det <- round(site_data[[site_detection_columns[i]]] * n_surveys)
      n_nondet <-
        round((1 - site_data[[site_detection_columns[i]]]) * n_surveys)
      site_data <- tibble::tibble(n_det = n_det, n_nondet = n_nondet,
                                  n_surveys = n_surveys, idx = seq_along(n_det))
      # prepare fitting fold data
      ## here we will prepare the data for model fitting by following the
      ## direct method described in: https://doi.org/10.1214/15-AOAS812
      ## note that this paper reports pretty poor performance for this
      ## method, but that is because they do not have observation-level
      ## training data (unlike here, where we have data for each site).
      ## This method is also conceptually similar to the EM algorithm outlined
      ## by Magder and Hughes
      ## (https://doi.org/10.1093/oxfordjournals.aje.a009251)
      ## except that it involves a single iteration, as opposed to iterating
      ## until parameter convergence.
      ## essentially, this approach involves using weights to account for
      ## imperfect detection during model fitting
      ## (similar in spirit to: https://doi.org/10.1002/env.2446)
      train_data <- site_data[f[[i]]$train[[k]], , drop = FALSE]
      prior_prob_pres <- vapply(
        seq_len(nrow(train_data)), FUN.VALUE = numeric(1), function(j) {
        prior_probability_of_occupancy(
          train_data$n_det[j], train_data$n_nondet[j], sens, spec, 0.5)
      })
      y_fit <- c(rep(1, nrow(train_data)), rep(0, nrow(train_data)))
      x_fit <- site_env_data[c(train_data$idx, train_data$idx), ,
                               drop = FALSE]
      ## note: multiply training weights by 100 to avoid floating point issues
      w_fit <- c(prior_prob_pres, 1 - prior_prob_pres) * 100
      # prepare train fold data (for performance stats)
      y_train <- c(rep(1, nrow(train_data)), rep(0, nrow(train_data)))
      x_train <- site_env_data[c(train_data$idx, train_data$idx), ,
                               drop = FALSE]
      w_train <- c(train_data$n_det / train_data$n_surveys,
                  train_data$n_nondet / train_data$n_surveys)
      # prepare test fold data (for model evaluation + performance stats)
      test_data <- site_data[f[[i]]$test[[k]], , drop = FALSE]
      y_test <- c(rep(1, nrow(test_data)), rep(0, nrow(test_data)))
      x_test <- site_env_data[c(test_data$idx, test_data$idx), ,
                               drop = FALSE]
      w_test <- c(test_data$n_det / test_data$n_surveys,
                  test_data$n_nondet / test_data$n_surveys)
      # return data
      list(fit = list(y = y_fit, x = x_fit, w = w_fit),
           train = list(y = y_train, x = x_train, w = w_train),
           test = list(y = y_test, x = x_test, w = w_test))
    })
  })

  # tune and fit models
  m <- lapply(seq_len(nrow(feature_data)), function(i) {
    tune_model(
      data = d[[i]], folds = f[[i]],
      survey_sensitivity =
        feature_data[[feature_survey_sensitivity_column]][i],
      survey_specificity =
        feature_data[[feature_survey_specificity_column]][i],
      parameters = xgb_tuning_parameters,
      early_stopping_rounds = xgb_early_stopping_rounds[i],
      n_rounds = xgb_n_rounds[i], n_folds = n_folds[i],
      n_threads = n_threads, verbose = verbose, seed = seed)
  })

  # assess models
  perf <- plyr::ldply(seq_len(nrow(feature_data)), function(i) {
    out <- plyr::ldply(seq_len(n_folds[i]), function(k) {
      ## extract fold training and test data
      m_k <- m[[i]]$models[[k]]
      nround_k <- m[[i]]$models[[k]]$best_iteration
      x_train_k <- d[[i]][[k]]$train$x
      x_test_k <- d[[i]][[k]]$test$x
      y_train_k <- d[[i]][[k]]$train$y
      y_test_k <- d[[i]][[k]]$test$y
      w_train_k <- d[[i]][[k]]$train$w
      w_test_k <- d[[i]][[k]]$test$w
      survey_sens <- feature_data[[feature_survey_sensitivity_column]][[i]]
      survey_spec <- feature_data[[feature_survey_specificity_column]][[i]]
      ## make predictions
      p_train_k <- c(withr::with_package("xgboost",
        stats::predict(m_k, x_train_k, ntreelimit = nround_k)))
      p_test_k <- c(withr::with_package("xgboost",
        stats::predict(m_k, x_test_k, ntreelimit = nround_k)))
      ## calculate performance
      perf_train_k <- rcpp_model_performance(
        y_train_k, p_train_k, w_train_k, survey_sens, survey_spec)
      perf_test_k <- rcpp_model_performance(
        y_test_k, p_test_k, w_test_k, survey_sens, survey_spec)
      ## calculate performance
      data.frame(
        train_tss = perf_train_k[[1]],
        train_sensitivity = perf_train_k[[2]],
        train_specificity = perf_train_k[[3]],
        test_tss = perf_test_k[[1]],
        test_sensitivity = perf_test_k[[2]],
        test_specificity = perf_test_k[[3]])
    })
    data.frame(feature = site_detection_columns[i],
               train_tss_mean = mean(out$train_tss),
               train_tss_std = stats::sd(out$train_tss),
               train_sensitivity_mean = mean(out$train_sensitivity),
               train_sensitivity_std = stats::sd(out$train_sensitivity),
               train_specificity_mean = mean(out$train_specificity),
               train_specificity_std = stats::sd(out$train_specificity),
               test_tss_mean = mean(out$test_tss),
               test_tss_std = stats::sd(out$test_tss),
               test_sensitivity_mean = mean(out$test_sensitivity),
               test_sensitivity_std = stats::sd(out$test_sensitivity),
               test_specificity_mean = mean(out$test_specificity),
               test_specificity_std = stats::sd(out$test_specificity),
               stringsAsFactors = FALSE)
  })
  perf <- tibble::as_tibble(perf)

  # make model predictions
  pred <- vapply(seq_len(nrow(feature_data)),
                 FUN.VALUE = numeric(nrow(site_env_data)), function(i) {
    nr <- m[[i]]$parameters$nround
    rowMeans(
      vapply(m[[i]]$models, FUN.VALUE = numeric(nrow(site_env_data)),
             function(x) {
      c(withr::with_package("xgboost", stats::predict(
        x, site_env_data, ntreelimit = x$best_iteration)))
    }))
  })
  colnames(pred) <- site_detection_columns
  pred <- tibble::as_tibble(pred)

  # return results
  list(parameters = lapply(m, `[[`, "parameters"),
       predictions = pred, performance = perf)
}

#' @noRd
tune_model <- function(data, folds, survey_sensitivity, survey_specificity,
  parameters, early_stopping_rounds, n_rounds, n_folds, n_threads, verbose,
  seed) {
  # assert arguments are valid
  assertthat::assert_that(
    isTRUE(n_folds == length(folds$train)),
    isTRUE(n_folds == length(folds$test)),
    assertthat::is.count(n_folds), assertthat::noNA(n_folds),
    assertthat::is.number(survey_sensitivity),
    assertthat::is.number(survey_specificity),
    assertthat::is.count(early_stopping_rounds),
    assertthat::noNA(early_stopping_rounds),
    assertthat::is.count(n_threads),
    assertthat::is.count(seed),
    assertthat::noNA(n_threads),
    is.list(parameters),
    assertthat::is.flag(verbose), assertthat::noNA(verbose))

  # create full parameters
  ## generate all combinations
  full_parameters <- do.call(expand.grid, parameters)
  attr(full_parameters, "out.attrs") <- NULL
  full_parameters <- lapply(full_parameters, function(x) {
    if (is.factor(x)) {
      return(as.character(x))
    }
    x
  })
  full_parameters <- data.frame(full_parameters, stringsAsFactors = FALSE)

  ## add objective if missing
  if (is.null(full_parameters$objective)) {
    full_parameters$objective <- "binary:logistic"
    warning(paste("no objective specified for model fitting,",
                  "assuming binary:logistic"))
  } else {
    assertthat::assert_that(length(unique(full_parameters$objective)) == 1,
      msg = "only one objective can be specified")
  }

  # calculate scale_pos_weight
  spw <- vapply(seq_len(n_folds), FUN.VALUE = numeric(1), function(k) {
    sum(data[[k]]$train$y < 0.5) / sum(data[[k]]$train$y > 0.5)
  })
  assertthat::assert_that(all(is.finite(spw)))

  # find best tuning parameters using k-fold cross validation
  ## fit models using all parameters combinations
  is_parallel <- (n_threads > 1) && (nrow(full_parameters) > 1)
  if (is_parallel) {
    cl <- parallel::makeCluster(n_threads, "FORK")
    doParallel::registerDoParallel(cl)
  }
  cv <- plyr::ldply(seq_len(nrow(full_parameters)), .parallel = is_parallel,
                    function(i)  {
    ## extract tuning parameters
    p <- as.list(full_parameters[i, , drop = FALSE])
    ## run cross-validation
    cv <- lapply(seq_len(n_folds), function(k) {
      ### prepare data for xgboost (note we use the fit data not the train data)
      dtrain <- xgboost::xgb.DMatrix(
        data[[k]]$fit$x, label = data[[k]]$fit$y,
        weight = data[[k]]$fit$w)
      dtest <- xgboost::xgb.DMatrix(
        data[[k]]$test$x, label = data[[k]]$test$y,
        weight = data[[k]]$test$w)
      ### prepare evaluation function
      curr_feval_tss <- make_feval_tss(survey_sensitivity, survey_specificity)
      ### prepare arguments for xgboost call
      args <- list(data = dtrain, verbose = FALSE, scale_pos_weight = spw[k],
                   watchlist = list(test = dtest), eval_metric = curr_feval_tss,
                   maximize = TRUE, nrounds = n_rounds, nthread = 1,
                   early_stopping_rounds = early_stopping_rounds)
      args <- append(args, p)
      ### fit model
      withr::with_seed(seed, {
        model <- do.call(what = xgboost::xgb.train, args)
      })
      ### evaluate model
      yhat_test <- c(withr::with_package("xgboost",
          stats::predict(
            model, data[[k]]$test$x, ntreelimit = model$best_iteration)))
      perf <- rcpp_model_performance(
        data[[k]]$test$y, yhat_test, data[[k]]$test$w,
        survey_sensitivity, survey_specificity)[[1]]
      ## check that model evaluations are consistent
      assertthat::assert_that(abs(perf - model$best_score) < 1e-5)
      ### return result
      list(eval = perf, model = model)
    })
    ## store the model performance
    tibble::tibble(
      eval = mean(vapply(cv, `[[`, numeric(1), "eval")),
      models = list(lapply(cv, `[[`, "model")))
  })
  if (is_parallel) {
    cl <- parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
  }
  ## determine best parameters for i'th species
  cv <- tibble::as_tibble(cv)
  j <- which.max(cv$eval)
  # return best parameters and models
  best_params <- as.list(full_parameters[j, , drop = FALSE])
  best_params$scale_pos_weight <- list(spw)
  best_params$objective <- as.character(best_params$objective)
  list(parameters = best_params, models = cv$models[[j]])
}

#' @noRd
make_feval_tss <- function(sens, spec) {
  function(preds, dtest) {
    labels <- xgboost::getinfo(dtest, "label")
    wts <- xgboost::getinfo(dtest, "weight")
    assertthat::assert_that(
      any(labels >= 0.5), any(labels < 0.5),
      msg = "test labels need at least one presence and one absence")
    assertthat::assert_that(all(preds >= 0), all(preds <= 1),
      msg = "xgboost predictions are not between zero and one")
    value <- rcpp_model_performance(labels, preds, wts, sens, spec)
    list(metric = "tss", value = value[[1]])
  }
}
