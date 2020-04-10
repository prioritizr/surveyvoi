#' @include internal.R
NULL

#' Tune occupancy models
#'
#' Identify suitable tuning parameters for fitting models to predict
#' site occupancy.
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
#'  parameter values for fitting models. Note this argument must
#'  have a \code{nrounds} element (see example below).
#'  Valid parameters include: \code{"max_depth"}, \code{"eta"},
#'  \code{"lambda"}, \code{"subsample"},
#'  \code{"colsample_bytree"}, and \code{"objective"}.
#'  See documentation for the \code{params} argument in
#'  \code{\link[xgboost]{xgb.train}} for more information.
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
#' @param nrounds \code{numeric} model rounds for model fitting.
#'   See \code{\link[xgboost]{xgboost}} for more information.
#'   Defaults to 10.
#'
#' @param early_stopping_rounds \code{numeric} model rounds for parameter
#'   tuning. See \code{\link[xgboost]{xgboost}} for more information.
#'   Defaults to 10.
#'
#' @param site_weight_columns \code{character} name of columns in
#'  \code{site_data} containing weights for model fitting. These columns must
#'  contain \code{numeric} values. No missing (\code{NA}) values are
#'  permitted. Defaults to \code{NULL} such that all data are given
#'  equal weight when fitting models.
#'
#' @param n_threads \code{integer} number of threads to use for fitting models.
#'   Defaults to 1.

#' @param seed \code{integer} seed for random number generation.
#'   Defaults to 500.
#'
#' @param verbose \code{logical} indicating if information should be
#'   printed during computations. Defaults to \code{FALSE}.

#' @details A random search method is used to tune the model
#'  parameters. For a given set of tuning parameters, models are fit
#'  using k-fold cross-validation (via \code{\link[xgboost]{xgb.cv}}) and the
#'  binary classification error rate calculated using the data held
#'  out from each fold is used to quantify the performance. In the event
#'  that the number of folds for a given feature exceeds the number of sites
#'  where it was recorded as present or absent, then a randomly selected
#'  site is added to each fold to ensure that each fold has at least one
#'  present site and one absent site. The models are
#'  fit using the \code{early_stopping_rounds} parameter to reduce time-spent
#'  tuning models. They are also fit using an automatically calibrated
#'  \code{scale_pos_weight} \code{\link[xgboost]{xgboost}} parameter
#'  (i.e. number of absence divided by number of presences per feature)
#'  to account for unbalanced datasets. Note that each feature
#'  is run with the same seed to ensure reproducibility.
#'
#' @return \code{list} of \code{list} objects containing the best tuning
#'   parameters found for each feature.
#'
#' @examples
#' # set seed
#' set.seed(100)
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
#' # identify suitable tuning parameters for each feature,
#' # note that we use 10 iterations here so that the example finishes quickly,
#' # you would probably want something like 1000+
#' parameters <- tune_occupancy_models(
#'    x, paste0("f", seq_len(2)), paste0("e", seq_len(3)),
#'    n_folds = rep(5, 2), n_random_search_iterations = 10,
#'    early_stopping_rounds = 5, nrounds = 1000, parameters = all_parameters,
#'    n_threads = 1)
#'
#' # print best found parameters for each feature
#' print(parameters)
#'
#' @export
tune_occupancy_models <- function(
  site_data, site_occupancy_columns, site_env_vars_columns, parameters,
  n_folds = rep(ifelse(nrow(site_data) > 1000, 10, 5), site_occupancy_columns),
  n_random_search_iterations = 10000, early_stopping_rounds = 10,
  n_rounds = 1000, site_weight_columns = NULL, n_threads = 1, seed = 500,
  verbose = FALSE) {
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
    assertthat::is.number(n_random_search_iterations),
    assertthat::noNA(n_random_search_iterations),
    assertthat::is.number(early_stopping_rounds),
    assertthat::noNA(early_stopping_rounds),
    assertthat::is.number(n_rounds), assertthat::noNA(n_rounds),
    assertthat::is.number(seed), assertthat::noNA(seed),
    assertthat::is.number(n_threads), assertthat::noNA(n_threads),
    is.list(parameters))
  param_names <- c("max_depth", "eta", "lambda", "subsample",
                   "colsample_bytree", "objective")
  assertthat::assert_that(
    all(names(parameters) %in% param_names),
    msg = paste0("unrecognised parameters(s) in argument to parameters:",
      paste(setdiff(names(parameters), param_names), collapse = ", ")))
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

  # fit models separately for each species
  m <- lapply(seq_along(site_occupancy_columns), function(f) {
    withr::with_seed(seed, {
      tune_model(y = site_data[[site_occupancy_columns[f]]],
                 x = site_env_data, w = site_weight_data[, f],
                 n_folds = n_folds[f],
                 n_random_search_iterations = n_random_search_iterations,
                 early_stopping_rounds = early_stopping_rounds,
                 n_rounds = n_rounds, parameters = parameters,
                 n_threads = n_threads, seed = seed, verbose = verbose)
    })
  })
  # return results
  m
}

#' @noRd
tune_model <- function(y, x, w, n_folds = 10, n_random_search_iterations = 5000,
                      early_stopping_rounds = 5, n_rounds = 1000,
                      parameters, n_threads = 1, seed = 500, verbose = FALSE) {
  # assert arguments are valid
  assertthat::assert_that(
    is.numeric(y),
    is.matrix(x), ncol(x) > 0, nrow(x) == length(y),
    is.numeric(w), length(w) == length(y),
    assertthat::is.count(n_folds), assertthat::noNA(n_folds),
    assertthat::is.count(n_random_search_iterations),
    assertthat::noNA(n_random_search_iterations),
    assertthat::is.count(early_stopping_rounds),
    assertthat::noNA(early_stopping_rounds),
    assertthat::is.count(n_threads),
    assertthat::noNA(n_threads),
    is.list(parameters),
    assertthat::is.flag(verbose), assertthat::noNA(verbose))
  # init
  assertthat::assert_that(!all(is.na(y)),
    msg = "a species is missing data for all planning units")
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
  full_parameters <- full_parameters[sample.int(nrow(full_parameters),
                                                n_random_search_iterations), ,
                                                drop = FALSE]
  # exclude NA data
  train_idx <- which(!is.na(y))
  y <- y[train_idx]
  x <- x[train_idx, , drop = FALSE]
  w <- w[train_idx]
  # create folds
  withr::with_seed(parameters$seed, {
    xgb_folds <- create_folds(y, n = n_folds, seed = seed)
  })
  # find best tuning parameters using k-fold cross validation
  ## fit models using all parameters combinations
  is_parallel <- n_threads > 1
  if (is_parallel) {
    cl <- parallel::makeCluster(n_threads, "FORK")
    doParallel::registerDoParallel(cl)
  }
  cv <- plyr::ldply(seq_len(nrow(full_parameters)), .parallel = is_parallel,
                    function(i)  {
    xgb_data <- xgboost::xgb.DMatrix(x, label = y, weight = w)
    xgboost_cv(
      data = xgb_data, parameters = full_parameters[i, , drop = FALSE],
      folds = xgb_folds, n_rounds = n_rounds,
      early_stopping_rounds = early_stopping_rounds, seed = seed,
      verbose = verbose)
  })
  if (is_parallel) {
    cl <- parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
  }
  ## determine best parameters for i'th species
  cv <- tibble::as_tibble(cv)
  k <- which.max(cv$test_eval_mean)
  o1 <<- cv[k, ]
  # return best parameters
  list(objective = as.character(full_parameters$objective[[k]]),
       seed = seed,
       max_depth = full_parameters$max_depth[[k]],
       eta = full_parameters$eta[[k]],
       lambda = full_parameters$lambda[[k]],
       subsample = full_parameters$subsample[[k]],
       colsample_bytree = full_parameters$colsample_bytree[[k]],
       nrounds = cv$nrounds[[k]],
       scale_pos_weight = cv$scale_pos_weight[[k]])
}

xgboost_cv <- function(data, parameters, folds, n_rounds,
                       early_stopping_rounds, seed, verbose) {
  cv <- plyr::ldply(seq_along(folds[[1]]), function(i) {
    ## extract data
    train_data <- xgboost::slice(data, folds$train[[i]])
    test_data <- xgboost::slice(data, folds$test[[i]])
    ## calculate scale pos weight
    spw <- sum(xgboost::getinfo(train_data, "label") < 0.5) /
           sum(xgboost::getinfo(train_data, "label") > 0.5)
    ## fit model
    withr::with_seed(seed, {
      m <- xgboost::xgb.train(
        params = list(
          objective = as.character(parameters$objective),
          max_depth = parameters$max_depth, eta = parameters$eta, nthread = 1,
          lambda = parameters$lambda, subsample = parameters$subsample,
          colsample_bytree = parameters$colsample_bytree),
        data = train_data,
        watchlist = list(train = train_data, eval = test_data),
        eval_metric = "auc", metric = "auc",
        nrounds = n_rounds, early_stopping_rounds = early_stopping_rounds,
        scale_pos_weight = spw, prediction = FALSE, showsd = TRUE,
        verbose = verbose)
    })
    ## make predictions
    p_train <- withr::with_package("xgboost",
      stats::predict(m, train_data, ntreelimit = m$best_ntreelimit))
    p_test <- withr::with_package("xgboost",
      stats::predict(m, test_data, ntreelimit = m$best_ntreelimit))
    ## return result
    tibble::tibble(
      nrounds = m$best_iteration,
      scale_pos_weight = spw,
      eval_train = weighted_auc(
        actual = xgboost::getinfo(train_data, "label"),
        predicted = p_train,
        weights = xgboost::getinfo(train_data, "weight")),
      eval_test = weighted_auc(
        actual = xgboost::getinfo(test_data, "label"),
        predicted = p_test,
        weights = xgboost::getinfo(test_data, "weight")))
  })
  tibble::tibble(train_eval_mean = mean(cv$eval_train),
                 train_eval_sd = sd(cv$eval_train),
                 test_eval_mean = mean(cv$eval_test),
                 test_eval_sd = sd(cv$eval_test),
                 nrounds = list(cv$nrounds),
                 scale_pos_weight = list(cv$scale_pos_weight))
}
