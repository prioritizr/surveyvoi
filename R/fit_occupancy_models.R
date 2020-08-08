#' @include internal.R
NULL

#' Fit models to estimate probability of occupancy
#'
#' Estimate probability of occupancy for a set of features in a set of
#' planning units.
#'
#' @param site_data \code{\link[sf]{sf}} object with site data.
#'
#' @param site_occupancy_columns \code{character} names of \code{numeric}
#'   columns in the
#'   argument to \code{site_data} that contain presence/absence data.
#'   Each column should correspond to a different feature, and contain
#'   binary presence/absence data (zeros or ones) indicating if the
#'   feature was detected in a previous survey or not. If a site has not
#'   been surveyed before, then missing (\code{NA}) values should be used.
#'
#' @param site_env_vars_columns \code{character} names of columns in the
#'   argument to \code{site_data} that contain environmental information
#'   for fitting updated occupancy models based on possible survey outcomes.
#'   Each column should correspond to a different environmental variable,
#'   and contain \code{numeric}, \code{factor}, or \code{character} data.
#'   No missing (\code{NA}) values are permitted in these columns.
#'
#' @param xgb_tuning_parameters \code{list} object containing the candidate
#'  parameter values for fitting models. Valid parameters include:
#'  \code{"max_depth"}, \code{"eta"}, \code{"lambda"},
#'  \code{"min_child_weight"}, \code{"subsample"}, \code{"colsample_by_tree"},
#'  \code{"objective"}. See documentation for the \code{params} argument in
#'  \code{\link[xgboost]{xgb.train}} for more information.
#'
#' @param xgb_early_stopping_rounds \code{numeric} model rounds for parameter
#'   tuning. See \code{\link[xgboost]{xgboost}} for more information.
#'   Defaults to 100 for each feature.
#'
#' @param xgb_n_rounds \code{numeric} model rounds for model fitting
#'   See \code{\link[xgboost]{xgboost}} for more information.
#'   Defaults to 1000 for each feature.
#'
#' @param xgb_n_folds \code{numeric} number of folds to split the training
#'   data into when fitting models for each feature.
#'   Defaults to 5 for each feature.
#'
#' @param site_weight_columns \code{character} name of columns in
#'  \code{site_data} containing weights for model fitting. These columns must
#'  contain \code{numeric} values. No missing (\code{NA}) values are
#'  permitted. Defaults to \code{NULL} such that all data are given
#'  equal weight when fitting models.
#'
#' @param n_threads \code{integer} number of threads to use for parameter
#'   tuning. Defaults to 1.
#'
#' @param seed \code{integer} initial random number generator state for model
#'   fitting. Defaults to 500.
#'
#' @param verbose \code{logical} indicating if information should be
#'   printed during computations. Defaults to \code{FALSE}.
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
#'  k-fold cross-validation (set via argument to \code{xgb_n_folds}). The
#'  training and evaluation folds are constructed
#'  in such a manner as to ensure that each training and evaluation
#'  fold contains at least one presence and one absence observation.
#'
#'  \item A grid search method is used to tune the model parameters. The
#'  candidate values for each parameter (specified via \code{parameters}) are
#'  used to generate a full set of parameter combinations, and these
#'  parameter combinations are subsequently used for tuning the models.
#'  To account for unbalanced datasets, the
#'  \code{scale_pos_weight} \code{\link[xgboost]{xgboost}} parameter
#'  is calculated as the mean value across each of the training folds
#'  (i.e. number of absence divided by number of presences per feature).
#'  For a given parameter combination, models are fit using k-fold cross-
#'  validation (via \code{\link[xgboost]{xgb.cv}}) -- using the previously
#'  mentioned training and evaluation folds -- and the True Skill Statistic
#'  (TSS) calculated using the data held out from each fold is
#'  used to quantify the performance (i.e. \code{"test_tss_mean"} column in
#'  output). These models are also fitted using the
#'  \code{early_stopping_rounds} parameter to reduce time-spent
#'  tuning models. If relevant, they are also fitted using the supplied weights
#'  (per by the argument to \code{site_weights_data}). After exploring the
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
#'  \code{site_weights_data}). These performance values are calculated using
#'  the models' training and evaluation folds.
#'
#' }
#'
#' @return \code{list} object containing:
#' \describe{
#'
#' \item{parameters}{\code{list} of \code{list} objects containing the best
#' tuning parameters for each feature.}
#'
#' \item{predictions}{\code{\link[tibble]{tibble}} object containing
#'  predictions for each feature.}
#'
#' \item{performance}{\code{\link[tibble]{tibble}} object containing the
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
#'  }
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
#'  x <- simulate_site_data(200, 2, 0.5, n_env_vars = 3)
#'
#' # create list of possible tuning parameters for modelling
#' all_parameters <- list(eta = seq(0.1, 0.5, length.out = 3),
#'                        lambda = 10 ^ seq(-1.0, 0.0, length.out = 3),
#'                        objective = "binary:logistic")
#'
#' # fit models
#' # note that we use 10 random search iterations here so that the example
#' # finishes quickly, you would probably want something like 1000+
#' results <- fit_occupancy_models(
#'    x, paste0("f", seq_len(2)), paste0("e", seq_len(3)),
#'    xgb_n_folds = rep(5, 2), xgb_early_stopping_rounds = rep(100, 2),
#'    xgb_tuning_parameters = all_parameters, n_threads = 1)
#'
#' # print best found model parameters
#' print(results$parameters)
#'
#' # print model predictions
#' print(results$predictions)
#'
#' # print model performance
#' print(results$performance)
#'
#' @export
fit_occupancy_models <- function(
  site_data, site_occupancy_columns, site_env_vars_columns,
  xgb_tuning_parameters,
  xgb_early_stopping_rounds = rep(100, length(site_occupancy_columns)),
  xgb_n_rounds = rep(1000, length(site_occupancy_columns)),
  xgb_n_folds = rep(5, length(site_occupancy_columns)),
  site_weight_columns = NULL, n_threads = 1, seed = 500, verbose = FALSE) {
  # assert that arguments are valid
  assertthat::assert_that(
    inherits(site_data, "sf"), nrow(site_data) > 0, ncol(site_data) > 0,
    is.character(site_occupancy_columns),
    length(site_occupancy_columns) > 0,
    assertthat::noNA(site_occupancy_columns),
    all(assertthat::has_name(site_data, site_occupancy_columns)),
    is.character(site_env_vars_columns),
    assertthat::noNA(site_env_vars_columns),
    all(assertthat::has_name(site_data, site_env_vars_columns)),
    is.numeric(xgb_n_folds), assertthat::noNA(xgb_n_folds),
    length(xgb_n_folds) == length(site_occupancy_columns),
    is.numeric(xgb_n_rounds), assertthat::noNA(xgb_n_rounds),
    length(xgb_n_rounds) == length(site_occupancy_columns),
    all(xgb_n_rounds > 0),
    is.numeric(xgb_early_stopping_rounds),
    assertthat::noNA(xgb_early_stopping_rounds),
    length(xgb_early_stopping_rounds) == length(site_occupancy_columns),
    all(xgb_early_stopping_rounds > 0),
    assertthat::is.count(seed),
    assertthat::noNA(seed),
    assertthat::is.count(n_threads), assertthat::noNA(n_threads),
    is.list(xgb_tuning_parameters))
  if (!is.null(site_weight_columns)) {
    assertthat::assert_that(
      is.character(site_weight_columns),
      identical(length(site_weight_columns), length(site_occupancy_columns)),
      all(assertthat::has_name(site_data, site_weight_columns)),
      assertthat::noNA(site_weight_columns))
    assertthat::assert_that(
      all(sapply(site_weight_columns,
                 function(x) is.numeric(site_data[[x]]))),
      msg = "site_data values in site_weight_columns must be numeric")
    assertthat::assert_that(
      all(sapply(site_weight_columns,
                 function(x) all(is.finite(site_data[[x]])))),
      msg = "site_data values in site_weight_columns must not be NA")
  }
  ## validate rij values
  assertthat::assert_that(
    all(sapply(site_occupancy_columns,
               function(x) is.numeric(site_data[[x]]))),
    msg = "site_data values in site_occupancy_columns must be numeric")
  assertthat::assert_that(
    all(sapply(site_occupancy_columns,
               function(x) max(site_data[[x]], na.rm = TRUE) <= 1)),
    msg = "site_data values in site_occupancy_columns must be <= 1")
  assertthat::assert_that(
    all(sapply(site_occupancy_columns,
               function(x) min(site_data[[x]], na.rm = TRUE) >= 0)),
    msg = "site_data values in site_occupancy_columns must be >= 0")
  ## validate environmental values
  assertthat::assert_that(
    all(sapply(site_env_vars_columns,
               function(x) is.numeric(site_data[[x]]))),
    msg = "site_data values in site_env_vars_columns must be numeric")
  assertthat::assert_that(
    all(sapply(site_env_vars_columns,
               function(x) assertthat::noNA(site_data[[x]]))),
    msg = "site_data values in site_env_vars_columns must not be NA")
  ## validate tuning parameters
  validate_xgboost_tuning_parameters(xgb_tuning_parameters)
  # drop geometry
  site_data <- sf::st_drop_geometry(site_data)

  # convert env data to model matrix format
  site_env_data <-
    as.matrix(stats::model.matrix(~ . - 1,
                                  data = site_data[, site_env_vars_columns]))

  # convert weight data to matrix format
  if (!is.null(site_weight_columns)) {
    site_weight_data <- as.matrix(site_data[, site_weight_columns])
  } else {
    site_weight_data <- matrix(1, ncol = length(site_occupancy_columns),
                               nrow = nrow(site_data))
  }

  # prepare data
  d <- lapply(seq_along(site_occupancy_columns), function(i) {
    # extract data
    y <- site_data[[site_occupancy_columns[i]]]
    x <- site_env_data
    w <- site_weight_data[, i]
    # exclude NA data
    train_idx <- which(!is.na(y))
    assertthat::assert_that(length(train_idx) > 0,
      msg = "a species is missing data for all planning units")
    # subset values
    y <- y[train_idx]
    x <- x[train_idx, , drop = FALSE]
    w <- w[train_idx]
    # return data
    list(y = y, x = x, w = w)
  })

  # prepare folds
  f <- lapply(seq_along(site_occupancy_columns), function(i) {
    withr::with_seed(seed, {
      create_folds(d[[i]]$y, n = xgb_n_folds[i])
    })
  })

  # tune and fit models
  m <- lapply(seq_along(site_occupancy_columns), function(i) {
    withr::with_seed(seed, {
      tune_model(data = d[[i]],
                 folds = f[[i]],
                 parameters = xgb_tuning_parameters,
                 early_stopping_rounds = xgb_early_stopping_rounds[i],
                 n_rounds = xgb_n_rounds[i],
                 n_folds = xgb_n_folds[i],
                 n_threads = n_threads,
                 verbose = verbose)
    })
  })

  # assess models
  perf <- plyr::ldply(seq_along(site_occupancy_columns), function(i) {
    out <- plyr::ldply(seq_len(xgb_n_folds[i]), function(k) {
      # extract fold training and test data
      m_k <- m[[i]]$models[[k]]
      nround_k <- m[[i]]$parameters$nrounds
      x_train_k <- d[[i]]$x[f[[i]]$train[[k]], , drop = FALSE]
      x_test_k <- d[[i]]$x[f[[i]]$test[[k]], , drop = FALSE]
      y_train_k <- d[[i]]$y[f[[i]]$train[[k]]]
      y_test_k <- d[[i]]$y[f[[i]]$test[[k]]]
      w_train_k <- d[[i]]$w[f[[i]]$train[[k]]]
      w_test_k <-d[[i]]$w[f[[i]]$test[[k]]]
      # make predictions
      p_train_k <- c(withr::with_package("xgboost",
        stats::predict(m_k, x_train_k, ntreelimit = nround_k)))
      p_test_k <- c(withr::with_package("xgboost",
        stats::predict(m_k, x_test_k, ntreelimit = nround_k)))
      ## calculate performance
      data.frame(
        train_tss = weighted_tss(y_train_k, p_train_k, w_train_k),
        train_sensitivity = sensitivity(y_train_k, p_train_k, w_train_k),
        train_specificity = specificity(y_train_k, p_train_k, w_train_k),
        test_tss = weighted_tss(y_test_k, p_test_k, w_test_k),
        test_sensitivity = sensitivity(y_test_k, p_test_k, w_test_k),
        test_specificity = specificity(y_test_k, p_test_k, w_test_k))
    })
    data.frame(feature = site_occupancy_columns[i],
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
  pred <- vapply(seq_along(site_occupancy_columns),
                 FUN.VALUE = numeric(nrow(site_env_data)), function(i) {
    nr <- m[[i]]$parameters$nround
    rowMeans(
      vapply(m[[i]]$models, FUN.VALUE = numeric(nrow(site_env_data)),
             function(x) {
      c(withr::with_package("xgboost", stats::predict(x, site_env_data,
                                                      ntreelimit = nr)))
    }))
  })
  colnames(pred) <- site_occupancy_columns
  pred <- tibble::as_tibble(pred)

  # return results
  list(parameters = lapply(m, `[[`, "parameters"),
       predictions = pred, performance = perf)
}

#' @noRd
tune_model <- function(data, folds, parameters, early_stopping_rounds,
  n_rounds, n_folds, n_threads, verbose) {
  # assert arguments are valid
  assertthat::assert_that(
    isTRUE(n_folds == length(folds$train)),
    isTRUE(n_folds == length(folds$test)),
    assertthat::is.count(n_folds), assertthat::noNA(n_folds),
    assertthat::is.count(early_stopping_rounds),
    assertthat::noNA(early_stopping_rounds),
    assertthat::is.count(n_threads),
    assertthat::noNA(n_threads),
    is.list(parameters),
    assertthat::is.flag(verbose), assertthat::noNA(verbose))

  # create full parameters
  ## generate all combinations
  full_parameters <- do.call(expand.grid, parameters)
  attr(full_parameters, "out.attrs") <- NULL
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
  spw <- mean(vapply(seq_len(n_folds), FUN.VALUE = numeric(1), function(k) {
    y <- data$y[folds$train[[k]]]
    sum(y < 0.5) / sum(y > 0.5)
  }))

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
    p$verbose <- 0
    p$nthread <- 1
    ## run cross-validation
    cv <- xgboost::xgb.cv(
      params = p,
      data = xgboost::xgb.DMatrix(data$x, label = data$y, weight = data$w),
      folds = folds$test, train_folds = folds$train,
      nrounds = n_rounds, early_stopping_rounds = early_stopping_rounds,
      feval = feval_tss, maximize = TRUE,
      scale_pos_weight = spw, prediction = FALSE, showsd = TRUE,
      verbose = verbose,
      callback = list(xgboost::cb.cv.predict(save_models = TRUE)))
    ## store the model performance
    tibble::tibble(
      eval = cv$evaluation_log$test_tss_mean[cv$best_ntreelimit],
      nrounds = cv$best_ntreelimit,
      models = list(cv$models))
  })
  if (is_parallel) {
    cl <- parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
  }
  ## determine best parameters for i'th species
  cv <- tibble::as_tibble(cv)
  k <- which.max(cv$eval)
  # return best parameters and models
  best_params <- as.list(full_parameters[k, , drop = FALSE])
  best_params$nrounds <- cv$nrounds[k]
  best_params$scale_pos_weight <- spw
  best_params$objective <- as.character(best_params$objective)
  list(parameters = best_params, models = cv$models[[k]])
}

#' @noRd
sensitivity <- function(actual, predicted, weights = rep(1, length(actual))) {
  assertthat::assert_that(
    is.numeric(actual), assertthat::noNA(actual),
    is.numeric(predicted), assertthat::noNA(predicted),
    is.numeric(weights), assertthat::noNA(weights),
    identical(length(actual), length(predicted)),
    identical(length(actual), length(weights)))
  # if there are no positives, then return if the predictions were all correct
  if (sum(actual > 0.5) == 0)
    return(as.numeric(all(round(actual) == round(predicted))))
  # calculate weighted sensitivity
  sum(weights * ((predicted >= 0.5) & (actual >= 0.5))) /
  sum(weights * (actual >= 0.5))
}

#' @noRd
specificity <- function(actual, predicted, weights = rep(1, length(actual))) {
  assertthat::assert_that(
    is.numeric(actual), assertthat::noNA(actual),
    is.numeric(predicted), assertthat::noNA(predicted),
    is.numeric(weights), assertthat::noNA(weights),
    identical(length(actual), length(predicted)),
    identical(length(actual), length(weights)))
  # if there are no negatives, then return if the predictions were all correct
  if (sum(actual < 0.5) == 0)
    return(as.numeric(all(round(actual) == round(predicted))))
  # calculate weighted specificity
  sum(weights * ((predicted < 0.5) & (actual < 0.5))) /
  sum(weights * (actual < 0.5))
}

#' @noRd
weighted_tss <- function(actual, predicted, weights = rep(1, length(actual))) {
  assertthat::assert_that(
    is.numeric(actual), assertthat::noNA(actual),
    is.numeric(predicted), assertthat::noNA(predicted),
    is.numeric(weights), assertthat::noNA(weights),
    identical(length(actual), length(predicted)),
    identical(length(actual), length(weights)))
  specificity(actual, predicted, weights) +
  sensitivity(actual, predicted, weights) - 1
}

#' @noRd
feval_tss <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  wts <- xgboost::getinfo(dtrain, "weight")
  value <- weighted_tss(labels, plogis(preds), wts)
  list(metric = "tss", value = value)
}
