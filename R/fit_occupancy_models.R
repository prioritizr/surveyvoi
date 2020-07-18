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
#' @param parameters \code{list} object containing the candidate
#'  parameter values for fitting models. Valid parameters include:
#'  \code{"max_depth"}, \code{"eta"}, \code{"lambda"},
#'  \code{"min_child_weight"}, \code{"subsample"}, \code{"colsample_by_tree"},
#'  \code{"objective"}. See documentation for the \code{params} argument in
#'  \code{\link[xgboost]{xgb.train}} for more information.
#'
#' @param early_stopping_rounds \code{numeric} model rounds for parameter
#'   tuning. See \code{\link[xgboost]{xgboost}} for more information.
#'   Defaults to 100.
#'
#' @param tree_method \code{character} method used for constructing trees.
#'   Available options are: \code{"auto"}, \code{"exact"}, \code{"approx"},
#'   and \code{"hist"}. Defaults to \code{"auto"}.
#'
#' @param n_rounds \code{numeric} model rounds for model fitting
#'   See \code{\link[xgboost]{xgboost}} for more information.
#'   Defaults to 1000.
#'
#' @param n_folds \code{numeric} number of folds to split the training
#'   data into when fitting models for each feature.
#'   Defaults to 10 if there are more than 1000
#'   rows in the argument \code{site_data}, otherwise defaults to 5.
#'   Note that if the maximum number of minority cases (e.g. 1 presence and
#'   1000 absences) is fewer than the number of folds then one minority case
#'   is sampled per fold.
#'
#' @param n_random_search_iterations \code{numeric} number of random search
#'    iterations to use when tuning model parameters. Defaults to 10000.
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
#'  k-fold cross-validation (set via argument to \code{n_folds}). The training
#'  and evaluation folds are constructed
#'  in such a manner as to ensure that each training and evaluation
#'  fold contains at least one presence and one absence observation.
#'
#'  \item A random search method is used to tune the model parameters. The
#'  candidate values for each parameter (specified via \code{parameters}) are
#'  used to generate a full set of parameter combinations, and then
#'  a subset of these parameter combinations is randomly selected for
#'  the tuning procedure (specified via \code{n_random_search_iterations}).
#'  To account for unbalanced datasets, the
#'  \code{scale_pos_weight} \code{\link[xgboost]{xgboost}} parameter
#'  is calculated as the mean value across each of the training folds
#'  (i.e. number of absence divided by number of presences per feature).
#'  For a given parameter combination, models are fit using k-fold cross-
#'  validation (via \code{\link[xgboost]{xgb.cv}}) -- using the previously
#'  mentioned training and evaluation folds -- and the average area under the
#'  curve (AUC) statistic calculated using the data held out from each fold is
#'  used to quantify the performance. These models are also fitted using the
#'  \code{early_stopping_rounds} parameter to reduce time-spent
#'  tuning models. If relevant, they are also fitted using the supplied weights
#'  (per by the argument to \code{site_weights_data}).After exploring the
#'  subset of parameter combinations, the best parameter combination is
#'  identified, and the associated parameter values and models are stored for
#'  later use.
#'
#'  \item The cross-validation models associated with the best parameter
#'   combination are used to generate predict the average probability that the
#'   feature occupies each site. These predictions include sites that have
#'   been surveyed before, and also sites that have not been surveyed before.
#'
#'  \item The performance of the cross-validation models is evaluated.
#'  Specifically, the AUC, sensitivity, and specificity statistics are
#'  calculated
#'  (if relevant, weighted by the argument to \code{site_weights_data}).
#'  These performance values are calculated using the models' evaluation folds.
#'
#' }
#'
#' @return \code{list} object containing:
#' \describe{
#' \item{parameters}{\code{list} of \code{list} objects containing the best
#' tuning parameters for each feature.}
#'
#' \item{predictions}{\code{\link[tibble]{tibble}} object containing
#'  predictions for each feature.}
#'
#' \item{performance}{\code{\link[tibble]{tibble}} object containing the
#'  performance of the best cross-validation models for each feature.}
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
#' all_parameters <- list(max_depth = seq(1, 10, 1),
#'                        eta = seq(0.1, 0.5, 0.1),
#'                        lambda = 10 ^ seq(-1.0, 0.0, 0.25),
#'                        subsample = seq(0.5, 1.0, 0.1),
#'                        colsample_bytree = seq(0.4, 1.0, 0.1),
#'                        objective = "binary:logistic")
#'
#' # fit models
#' # note that we use 10 random search iterations here so that the example
#' # finishes quickly, you would probably want something like 1000+
#' results <- fit_occupancy_models(
#'    x, paste0("f", seq_len(2)), paste0("e", seq_len(3)),
#'    n_folds = rep(5, 2), n_random_search_iterations = 10,
#'    early_stopping_rounds = 100, parameters = all_parameters, n_threads = 1)
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
  site_data, site_occupancy_columns, site_env_vars_columns, parameters,
  early_stopping_rounds = 100, tree_method = "auto", n_rounds = 1000,
  n_folds = rep(5, length(site_occupancy_columns)),
  n_random_search_iterations = 10000,
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
    is.numeric(n_folds), assertthat::noNA(n_folds),
    length(n_folds) == length(site_occupancy_columns),
    assertthat::is.count(n_random_search_iterations),
    assertthat::noNA(n_random_search_iterations),
    assertthat::is.count(n_rounds),
    assertthat::noNA(n_rounds),
    assertthat::is.count(early_stopping_rounds),
    assertthat::noNA(early_stopping_rounds),
    assertthat::is.count(seed),
    assertthat::noNA(seed),
    assertthat::is.string(tree_method),
    assertthat::noNA(tree_method),
    tree_method %in% c("auto", "hist", "exact", "approx"),
    assertthat::is.count(n_threads), assertthat::noNA(n_threads),
    is.list(parameters),
    isTRUE(n_random_search_iterations <= prod(lengths(parameters))))
    param_names <- c("max_depth", "eta", "lambda", "subsample",
                     "colsample_bytree", "objective")
    assertthat::assert_that(
      all(names(parameters) %in% param_names),
      msg = paste("argument to parameters has unrecognised elements:",
                 paste(setdiff(names(parameters), param_names),
                       collapse = ", ")))
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
      create_folds(d[[i]]$y, n = n_folds[i])})
  })

  # tune and fit models
  m <- lapply(seq_along(site_occupancy_columns), function(i) {
    withr::with_seed(seed, {
      tune_model(data = d[[i]],
                 folds = f[[i]],
                 parameters = parameters,
                 early_stopping_rounds = early_stopping_rounds,
                 tree_method = tree_method,
                 n_rounds = n_rounds,
                 n_folds = n_folds[i],
                 n_random_search_iterations = n_random_search_iterations,
                 n_threads = n_threads,
                 verbose = verbose)
    })
  })

  # assess models
  perf <- plyr::ldply(seq_along(site_occupancy_columns), function(i) {
    out <- plyr::ldply(seq_len(n_folds[i]), function(k) {
      # extract data
      m_k <- m[[i]]$models[[k]]
      nround_k <- m[[i]]$parameters$nround
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
        train_auc = weighted_auc(y_train_k, p_train_k, w_train_k),
        train_sensitivity = sensitivity(y_train_k, p_train_k, w_train_k),
        train_specificity = specificity(y_train_k, p_train_k, w_train_k),
        test_auc = weighted_auc(y_test_k, p_test_k, w_test_k),
        test_sensitivity = sensitivity(y_test_k, p_test_k, w_test_k),
        test_specificity = specificity(y_test_k, p_test_k, w_test_k))
    })
    data.frame(feature = site_occupancy_columns[i],
               train_auc_mean = mean(out$train_auc),
               train_auc_std = sd(out$train_auc),
               train_sensitivity_mean = mean(out$train_sensitivity),
               train_sensitivity_std = sd(out$train_sensitivity),
               train_specificity_mean = mean(out$train_specificity),
               train_specificity_std = sd(out$train_specificity),
               test_auc_mean = mean(out$test_auc),
               test_auc_std = sd(out$test_auc),
               test_sensitivity_mean = mean(out$test_sensitivity),
               test_sensitivity_std = sd(out$test_sensitivity),
               test_specificity_mean = mean(out$test_specificity),
               test_specificity_std = sd(out$test_specificity),
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
  tree_method, n_rounds, n_folds, n_random_search_iterations, n_threads,
  verbose) {
  # assert arguments are valid
  assertthat::assert_that(
    isTRUE(n_folds == length(folds$train)),
    isTRUE(n_folds == length(folds$test)),
    assertthat::is.count(n_folds), assertthat::noNA(n_folds),
    assertthat::is.count(n_random_search_iterations),
    assertthat::noNA(n_random_search_iterations),
    assertthat::is.count(early_stopping_rounds),
    assertthat::noNA(early_stopping_rounds),
    assertthat::is.string(tree_method),
    assertthat::noNA(tree_method),
    assertthat::is.count(n_threads),
    assertthat::noNA(n_threads),
    is.list(parameters),
    assertthat::is.flag(verbose), assertthat::noNA(verbose))

  # create full parameters
  ## generate all combinations
  full_parameters <- do.call(expand.grid, parameters)
  ## add objective if missing
  if (is.null(full_parameters$objective)) {
    full_parameters$objective <- "binary:logistic"
    warning(paste("no objective specified for model fitting,",
                  "assuming binary:logistic"))
  } else {
    assertthat::assert_that(length(unique(full_parameters$objective)) == 1,
      msg = "only one objective can be specified")
  }
  ## subset tuning parameters to number of random search iterations
  full_parameters <-
    full_parameters[sample.int(nrow(full_parameters),
                               n_random_search_iterations), , drop = FALSE]

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
    p <- full_parameters[i, , drop = FALSE]
    ## run cross-validation
    cv <- xgboost::xgb.cv(
      params = list(objective = p$objective, verbose = 0,
                    max_depth = p$max_depth, eta = p$eta, nthread = 1,
                    lambda = p$lambda, subsample = p$subsample,
                    colsample_bytree = p$colsample_bytree,
                    tree_method = tree_method),
      data = xgboost::xgb.DMatrix(data$x, label = data$y, weight = data$w),
      folds = folds$test, train_folds = folds$train,
      nrounds = n_rounds, early_stopping_rounds = early_stopping_rounds,
      eval_metric = "auc", metrics = "auc", maximize = TRUE,
      scale_pos_weight = spw, prediction = FALSE, showsd = TRUE,
      verbose = verbose,
      callback = list(xgboost::cb.cv.predict(save_models = TRUE)))
    ## store the model performance
    tibble::tibble(
      eval = cv$evaluation_log$test_auc_mean[cv$best_ntreelimit],
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
  list(
    parameters = list(
      max_depth = full_parameters$max_depth[k],
      eta = full_parameters$eta[k],
      nrounds = cv$nrounds[k],
      scale_pos_weight = spw,
      lambda = full_parameters$lambda[k],
      tree_method = tree_method,
      subsample = full_parameters$subsample[k],
      colsample_bytree = full_parameters$colsample_bytree[k],
      objective = as.character(full_parameters$objective[k])),
    models = cv$models[[k]])
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
weighted_auc <- function(actual, predicted, weights = rep(1, length(actual))) {
  assertthat::assert_that(
    is.numeric(actual), assertthat::noNA(actual),
    is.numeric(predicted), assertthat::noNA(predicted),
    is.numeric(weights), assertthat::noNA(weights),
    identical(length(actual), length(predicted)),
    identical(length(actual), length(weights)))
  WeightedROC::WeightedAUC(WeightedROC::WeightedROC(predicted, actual, weights))
}
