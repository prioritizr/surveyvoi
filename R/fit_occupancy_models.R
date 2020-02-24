#' @include internal.R tune_occupancy_models.R
NULL

#' Fit models to estimate probability of occupancy
#'
#' Estimate probability of occupancy for a set of features in a set of
#' planning units.
#'
#' @inheritParams tune_occupancy_models
#'
#' @param parameters \code{list} of \code{list} objects containing the
#'  parameters for fitting models for each
#'  feature. Ideally, these parameters would be determined using the
#'  \code{\link{tune_occupancy_models}} function. Note that arguments must
#'  have a \code{nrounds} element (see example below).
#'  Valid parameters include: \code{"max_depth"}, \code{"eta"},
#'  \code{"nrounds"}, \code{"lambda"}, \code{"subsample"},
#'  \code{"colsample_bytree"}, and \code{"objective"}.
#'  See documentation for the \code{params} argument in
#'  \code{\link[xgboost]{xgb.train}} for more information.
#'
#' @details The models are fitted using the specified tuning parameters and
#'  automatically calibrated \code{scale_pos_weight}
#'  \code{\link[xgboost]{xgboost}} parameters
#'  (i.e. number of absence divided by number of presences per feature)
#'  to account for unbalanced datasets. After the models are fitted, they
#'  are used to make probability of occupancy predictions (via
#'  \code{\link[xgboost]{xgboost}}). The average sensitivity, specificity,
#'  and area under the curve (AUC) statistics are calculated using the data held
#'  are then reported to assess model performance. Note that each feature
#'  is run with the same seed to ensure reproducibility.
#'
#' @return \code{list} object containing the following elements:
#'  \describe{
#'  \item{predictions}{\code{\link[tibble]{tibble}} with estimated
#'    probabilities of each feature occupying each site. Rows
#'    correspond to sites and columns correspond to features.}
#'  \item{performance}{\code{\link[tibble]{tibble}} with model information.
#'    Each row corresponds to a different feature and column contain
#'    different information. Specifically, the columns contain:
#'    (\code{name}) \code{character} name of the features,
#'    (\code{sensitivity}) \code{numeric} model sensitivities,
#'    (\code{specificity}) \code{numeric} model specificities, and
#'    (\code{auc}) \code{numeric} model area under the curve (AUC) statistics.}
#'  }
#'
#' @examples
#' # set seed
#' set.seed(100)
#'
#' # simulate data for 200 sites, 2 features, and 3 environmental variables
#'  x <- simulate_site_data(200, 2, 0.5, n_env_vars = 3)
#'
#' # set parameters manually
#' xgb_params <- list(list(nrounds = 10, objective = "binary:logistic"),
#'                    list(nrounds = 5, objective = "binary:logistic"))
#'
#' # fit models
#' m <- fit_occupancy_models(x, c("f1", "f2"), c("e1", "e2", "e3"),
#'                           xgb_params)
#'
#' # preview models predictions
#' print(m$predictions)
#'
#' # preview models performance
#' print(m$performance)
#'
#' @export
fit_occupancy_models <- function(
  site_data, site_occupancy_columns, site_env_vars_columns, parameters,
  n_folds = rep(ifelse(nrow(site_data) > 1000, 10, 5),
                length(site_occupancy_columns)),
  site_weight_columns = NULL, n_threads = 1, verbose = FALSE) {
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
    is.list(parameters),
    identical(length(parameters), length(site_occupancy_columns)),
    assertthat::is.number(n_threads), assertthat::noNA(n_threads))
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

  # store seed so that each species is run with the same seed
  seed <- .Random.seed

  # fit models separately for each species
  m <- lapply(seq_along(site_occupancy_columns), function(f) {
    set.seed(seed)
    fit_model(y = site_data[[site_occupancy_columns[f]]],
              x = site_env_data,
              w = site_weight_data[, f],
              n_folds = n_folds[f],
              parameters = parameters[[f]],
              verbose = verbose)
  })
  # return results
  pij <- vapply(m, `[[`, numeric(nrow(site_data)), "predictions")
  colnames(pij) <- site_occupancy_columns
  perf <- tibble::as_tibble(plyr::ldply(m, `[[`, "performance"))
  perf$name <- site_occupancy_columns
  list(predictions = tibble::as_tibble(pij), performance = perf)
}

#' @noRd
fit_model <- function(y, x, w, n_folds = 10, parameters = list(), n_threads = 1,
                      verbose = FALSE) {
  # assert arguments are valid
  assertthat::assert_that(
    is.numeric(y),
    is.matrix(x), ncol(x) > 0, nrow(x) == length(y),
    is.numeric(y), length(w) == length(y),
    assertthat::is.count(n_folds), assertthat::noNA(n_folds),
    is.list(parameters),
    assertthat::is.count(n_threads), assertthat::noNA(n_threads))
  # init
  train_idx <- which(!is.na(y))
  assertthat::assert_that(length(train_idx) > 0,
    msg = "a species is missing data for all planning units")
  ## add objective if missing
  if (is.null(parameters$objective)) {
    parameters$objective <- "binary:logistic"
    warning(paste("no objective specified for model fitting,",
                  "assuming binary:logistic"))
  }
  # subset data for model fitting and evaluation
  y_sub <- y[train_idx]
  x_sub <- x[train_idx, , drop = FALSE]
  w_sub <- w[train_idx]
  # create folds
  xgb_folds <- create_folds(y_sub, n = n_folds)
  # generate model predictions k-fold cross-validation
  xgb_models <- lapply(seq_len(n_folds), function(i) {
    # create data for model fitting and tuning
    xgb_data <- xgboost::xgb.DMatrix(
      x_sub[xgb_folds$train[[i]], , drop = FALSE],
      label = y_sub[xgb_folds$train[[i]]],
      weight = w_sub[xgb_folds$train[[i]]])
    # calculate scale pos weight
    spw <- sum(xgboost::getinfo(xgb_data, "label") < 0.5) /
           sum(xgboost::getinfo(xgb_data, "label") > 0.5)
    # fit model
    do.call(xgboost::xgboost,
            append(parameters,
                   list(data = xgb_data,
                        nthread = n_threads,
                        metrics = "auc",
                        eval_metric = "auc",
                        scale_pos_weight = spw,
                        verbose = as.numeric(verbose))))
  })
  # calculate model performance
  fold_perf <- plyr::ldply(seq_along(xgb_folds), function(i) {
    ## extract data
    m <- xgb_models[[i]]
    ## create predictions
    p_train <-
      withr::with_package("xgboost",
        stats::predict(m, x_sub[xgb_folds$train[[i]], , drop = FALSE]))
    p_test <-
      withr::with_package("xgboost",
        stats::predict(m, x_sub[xgb_folds$test[[i]], , drop = FALSE]))
    ## subset observations
    y_train <- y_sub[xgb_folds$train[[i]]]
    y_test <- y_sub[xgb_folds$test[[i]]]
    ## calculate performance
    data.frame(
      train_sensitivity = sensitivity(y_train, p_train),
      train_specificity = specificity(y_train, p_train),
      train_auc = Metrics::auc(y_train, p_train),
      test_sensitivity = sensitivity(y_test, p_test),
      test_specificity = specificity(y_test, p_test),
      test_auc = Metrics::auc(y_test, p_test))
  })
  # calculate average predictions
  fold_pred <- lapply(xgb_models, function(m) {
    matrix(withr::with_package("xgboost", stats::predict(m, x)), nrow = 1)
  })
  # calculate average predictions across k-folds
  avg_pred <- colMeans(do.call(rbind, fold_pred))
  # calculate mean performance across k-folds
  avg_perf <- matrix(colMeans(fold_perf), nrow = 1)
  colnames(avg_perf) <- paste0(names(fold_perf), "_mean")
  # calculate std dev performance across k-folds
  std_perf <- matrix(apply(fold_perf, 2, sd), nrow = 1)
  colnames(std_perf) <- paste0(names(fold_perf), "_std")
  perf <- tibble::as_tibble(cbind(avg_perf, std_perf))
  ord1 <- sort(grep("train", names(perf), value = TRUE))
  ord2 <- sort(grep("test", names(perf), value = TRUE))
  perf <- perf[, c(ord1, ord2)]
  # return results
  list(predictions = avg_pred,
       performance = perf)
}

#' @noRd
sensitivity <- function(actual, predicted) {
  assertthat::assert_that(
    is.numeric(actual), assertthat::noNA(actual),
    is.numeric(predicted), assertthat::noNA(predicted))
  # if there are no positives, then return if the predictions were all correct
  if (sum(actual > 0.5) == 0)
    return(as.numeric(all(round(actual) == round(predicted))))
  sum((predicted > 0.5) & (actual > 0.5)) / sum(actual > 0.5)
}

#' @noRd
specificity <- function(actual, predicted) {
  assertthat::assert_that(
    is.numeric(actual), assertthat::noNA(actual),
    is.numeric(predicted), assertthat::noNA(predicted))
  # if there are no negatives, then return if the predictions were all correct
  if (sum(actual < 0.5) == 0)
    return(as.numeric(all(round(actual) == round(predicted))))
  sum((predicted < 0.5) & (actual < 0.5)) / sum(actual < 0.5)
}
