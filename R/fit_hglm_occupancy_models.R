#' @include internal.R
NULL

#' Fit hierarchical generalized linear models to predict occupancy
#'
#' Estimate probability of occupancy for a set of features in a set of
#' planning units. Models are fitted as hierarchical generalized linear models
#' that account for for imperfect detection (following Royle & Link 2006)
#' using JAGS (via \code{\link[runjags]{run.jags}}). To limit over-fitting,
#' spike-and-slab coefficients are used
#' to perform stochastic search variable selection (Kuo & Mallick 1998).
#'
#' @inheritParams fit_xgb_occupancy_models
#'
#' @param jags_n_iter \code{integer} total number of iterations per chain
#'   for MCMC analyses.
#'
#' @param jags_n_thin \code{integer} number of thinning iterations for MCMC
#'   analyses.
#'
#' @param jags_n_burnin \code{integer} number of warm up iterations for MCMC
#'   analyses.
#'
#' @param jags_n_chains \code{integer} total number of chains for MCMC analyses.
#'
#' @param jags_n_adapt \code{integer} number of adapting iterations for MCMC
#'   analyses.
#'
#' @details
#'  This function (i) prepares the data for model fitting,
#'  (ii) fits the models, and (iii) assesses the performance of the
#'   models. These analyses are performed separately for each
#'  feature. For a given feature:
#'
#'  \enumerate{
#'
#'  \item The data are prepared for model fitting by partitioning the data using
#'  k-fold cross-validation (set via argument to \code{n_folds}). The
#'  training and evaluation folds are constructed
#'  in such a manner as to ensure that each training and evaluation
#'  fold contains at least one presence and one absence observation.
#'
#'  \item A model for fit separately for each fold (see
#'  \code{inst/jags/model.jags} for model code). To assess convergence,
#'  the multi-variate potential scale reduction factor
#'  (MPSRF) statistic is calculated for each model.
#'
#'  \item The performance of the cross-validation models is evaluated.
#'  Specifically, the TSS, sensitivity, and specificity statistics are
#'  calculated (if relevant, weighted by the argument to
#'  \code{site_weights_data}). These performance values are calculated using
#'  the models' training and evaluation folds. To assess convergence,
#'  the maximum MPSRF statistic for the models fit for each feature
#'  is calculated.
#'
#' }
#'
#' @references
#' Kuo L & Mallick B (1998) Variable selection for regression models.
#' \emph{Sankhyā: The Indian Journal of Statistics, Series B}, 60: 65--81.
#'
#' Royle JA & Link WA (2006) Generalized site occupancy models allowing for
#' false positive and false negative errors. \emph{Ecology}, 87: 835--841.
#'
#' @return \code{list} object containing:
#' \describe{
#'
#' \item{models}{\code{list} of \code{list} objects containing the models.}
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
#'  \item{max_mpsrf}{maximum multi-variate potential scale reduction factor
#'    (MPSRF) value for the models. A MPSRF value less than 1.05 means that all
#'    coefficients in a given model have converged, and so a value less than
#'    1.05 in this column means that all the models fit for a given feature
#'    have successfully converged.}
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
#' # print JAGS model code
#' cat(readLines(system.file("inst/jags/model.jags", package = "surveyvoi")),
#'     sep = "\n")
#'
#' \dontrun{
#' # fit models
#' # note that we use a limited number of MCMC iterations so that the example
#' # finishes quickly, you would probably want something like 10,000+
#' results <- fit_hglm_occupancy_models(
#'    site_data, feature_data,
#'    c("f1", "f2"), c("n1", "n2"), c("e1", "e2", "e3"),
#'    "survey_sensitivity", "survey_specificity",
#'    n_folds = rep(5, 2),
#'    jags_n_iter = rep(500, 2), jags_n_burnin = rep(250, 2),
#'    jags_n_thin = rep(1, 2), jags_n_adapt = rep(50, 2),
#'    n_threads = 1)
#'
#' # print best found model parameters
#' print(results$parameters)
#'
#' # print model predictions
#' print(results$predictions)
#'
#' # print model performance
#' print(results$performance)
#' }
#' @export
fit_hglm_occupancy_models <- function(
  site_data, feature_data,
  site_detection_columns, site_n_surveys_columns,
  site_env_vars_columns,
  feature_survey_sensitivity_column, feature_survey_specificity_column,
  jags_n_iter = rep(150000, length(site_detection_columns)),
  jags_n_burnin = rep(100000, length(site_detection_columns)),
  jags_n_thin = rep(100, length(site_detection_columns)),
  jags_n_adapt = rep(1000, length(site_detection_columns)),
  jags_n_chains = rep(4, length(site_detection_columns)),
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
    ## n_folds
    is.numeric(n_folds),
    assertthat::noNA(n_folds),
    all(n_folds > 0),
    length(n_folds) == nrow(feature_data),
    ## jags_n_iter
    is.numeric(jags_n_iter),
    assertthat::noNA(jags_n_iter),
    length(jags_n_iter) == nrow(feature_data),
    all(jags_n_iter > 0),
    ## jags_n_thin
    is.numeric(jags_n_thin),
    assertthat::noNA(jags_n_thin),
    length(jags_n_thin) == nrow(feature_data),
    all(jags_n_thin > 0),
    ## jags_n_burnin
    is.numeric(jags_n_burnin),
    assertthat::noNA(jags_n_burnin),
    length(jags_n_burnin) == nrow(feature_data),
    all(jags_n_burnin > 0),
    all(jags_n_burnin < jags_n_iter),
    ## jags_n_chains
    is.numeric(jags_n_chains),
    assertthat::noNA(jags_n_chains),
    length(jags_n_chains) == nrow(feature_data),
    all(jags_n_chains > 0),
    ## jags_n_adapt
    is.numeric(jags_n_adapt),
    assertthat::noNA(jags_n_adapt),
    length(jags_n_adapt) == nrow(feature_data),
    all(jags_n_adapt > 0),
    ## seed
    assertthat::is.count(seed),
    assertthat::noNA(seed),
    ## n_threads
    assertthat::is.count(n_threads), assertthat::noNA(n_threads))
  ## validate survey data
  validate_site_detection_data(site_data, site_detection_columns)
  validate_site_n_surveys_data(site_data, site_n_surveys_columns)

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
      # prepare fitting fold data (actually used for model fitting)
      train_data <- site_data[f[[i]]$train[[k]], , drop = FALSE]
      y_fit <- matrix(-1, nrow = nrow(train_data),
                        ncol = max(train_data$n_surveys))
      for (j in seq_len(nrow(train_data))) {
        y_fit[j, seq_len(train_data$n_surveys[j])] <-
          c(rep(1, train_data$n_det[j]), rep(0, train_data$n_nondet[j]))
      }
      x_fit <- site_env_data[train_data$idx, , drop = FALSE]
      # prepare train fold data (only used for model evaluation)
      y_train <- c(rep(1, nrow(train_data)), rep(0, nrow(train_data)))
      x_train <- site_env_data[c(train_data$idx, train_data$idx), ,
                               drop = FALSE]
      w_train <- c(train_data$n_det / train_data$n_surveys,
                  train_data$n_nondet / train_data$n_surveys)
      # prepare test fold data
      test_data <- site_data[f[[i]]$test[[k]], , drop = FALSE]
      y_test <- c(rep(1, nrow(test_data)), rep(0, nrow(test_data)))
      x_test <- site_env_data[c(test_data$idx, test_data$idx), ,
                               drop = FALSE]
      w_test <- c(test_data$n_det / test_data$n_surveys,
                  test_data$n_nondet / test_data$n_surveys)
      # return data
      list(fit = list(y = y_fit, x = x_fit),
           train = list(y = y_train, x = x_train, w = w_train),
           test = list(y = y_test, x = x_test, w = w_test))
    })
  })

  # fit models
  ## create iterator
  model_cmbs <- plyr::ldply(seq_len(nrow(feature_data)), function(i) {
    expand.grid(spp = i, fold = seq_len(n_folds[i]),
                chain = seq_len(jags_n_chains[i]))
  })
  model_cmbs$idx <- seq_len(nrow(model_cmbs))
  ## prepare cluster
  if (n_threads > 1) {
    cl <- parallel::makeCluster(n_threads, "FORK")
    doParallel::registerDoParallel(cl)
  }
  ## main processing
  m_raw <- plyr::llply(
    seq_len(nrow(model_cmbs)),
    .parallel = n_threads > 1,
    .progress = ifelse(n_threads > 1, "none", "text"),
    function(cm) {
    i <- model_cmbs$spp[[cm]]
    f <- model_cmbs$fold[[cm]]
    ch <- model_cmbs$chain[[cm]]
    fit_hglm_model(
      data = d[[i]][[f]]$fit,
      survey_sensitivity =
        feature_data[[feature_survey_sensitivity_column]][i],
      survey_specificity =
        feature_data[[feature_survey_specificity_column]][i],
      n_iter = jags_n_iter[i],
      n_burnin = jags_n_burnin[i],
      n_thin = jags_n_thin[i],
      n_chains = 1,
      n_adapt = jags_n_adapt[i],
      verbose = verbose,
      seed = seed + i + f + ch)
  })
  ## kill cluster
  if (n_threads > 1) {
    doParallel::stopImplicitCluster()
    cl <- parallel::stopCluster(cl)
  }

  # combine models into single model object (quietly)
  m <- lapply(seq_len(nrow(feature_data)), function(i) {
    lapply(seq_len(n_folds[i]), function(k) {
      i <- model_cmbs$idx[model_cmbs$spp == i & model_cmbs$fold == k]
      runjags::combine.jags(m_raw[i])
   })
  })

  # assess models
  perf <- plyr::ldply(seq_len(nrow(feature_data)), function(i) {
    out <- plyr::ldply(seq_len(n_folds[i]), function(k) {
      ## extract fold training and test data
      x_train_k <- d[[i]][[k]]$train$x
      x_test_k <- d[[i]][[k]]$test$x
      y_train_k <- d[[i]][[k]]$train$y
      y_test_k <- d[[i]][[k]]$test$y
      w_train_k <- d[[i]][[k]]$train$w
      w_test_k <- d[[i]][[k]]$test$w
      survey_sens <- feature_data[[feature_survey_sensitivity_column]][[i]]
      survey_spec <- feature_data[[feature_survey_specificity_column]][[i]]
      ## make predictions
      p_train_k <- hglm_predict(m[[i]][[k]], x_train_k)
      p_test_k <- hglm_predict(m[[i]][[k]], x_test_k)
      ## calculate performance
      perf_train_k <- rcpp_model_performance(
        y_train_k, p_train_k, w_train_k, survey_sens, survey_spec)
      perf_test_k <- rcpp_model_performance(
        y_test_k, p_test_k, w_test_k, survey_sens, survey_spec)
      ## calculate performance
      data.frame(
        mpsrf = m[[i]][[k]]$psrf$mpsrf,
        train_tss = perf_train_k[[1]],
        train_sensitivity = perf_train_k[[2]],
        train_specificity = perf_train_k[[3]],
        test_tss = perf_test_k[[1]],
        test_sensitivity = perf_test_k[[2]],
        test_specificity = perf_test_k[[3]])
    })
    data.frame(feature = site_detection_columns[i],
               max_mpsrf = max(out$mpsrf, na.rm = TRUE),
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
  pred <- vapply(m, FUN.VALUE = numeric(nrow(site_env_data)), function(x) {
    rowMeans(vapply(x, FUN.VALUE = numeric(nrow(site_env_data)), function(y) {
      c(hglm_predict(y, site_env_data))
    }))
  })
  colnames(pred) <- site_detection_columns
  pred <- tibble::as_tibble(pred)

  # return results
  list(models = m,
       predictions = pred,
       performance = perf)
}

#' @noRd
fit_hglm_model <- function(data, survey_sensitivity, survey_specificity,
  n_iter, n_thin, n_burnin, n_chains, n_adapt, verbose, seed) {
  # assert arguments are valid
  assertthat::assert_that(
    is.list(data),
    assertthat::is.number(survey_sensitivity),
    assertthat::is.number(survey_specificity),
    assertthat::is.count(n_iter),
    assertthat::noNA(n_iter),
    assertthat::is.count(n_burnin),
    assertthat::noNA(n_burnin),
    assertthat::is.count(n_adapt),
    assertthat::noNA(n_adapt),
    assertthat::is.count(n_chains),
    assertthat::noNA(n_chains),
    assertthat::is.count(seed),
    assertthat::is.flag(verbose), assertthat::noNA(verbose))
  # prepare JAGS data
  jags_data <- list(
    n_vars = ncol(data$x),
    train_n_sites = nrow(data$x),
    train_n_obs = vapply(
      seq_len(nrow(data$y)), FUN.VALUE = integer(1), function(i) {
        sum(data$y[i, ] > -0.5)
      }),
    sens = survey_sensitivity,
    spec = survey_specificity,
    train_model_matrix = data$x,
    train_obs = data$y)
  # run JAGS
  withr::with_seed(seed, {
    runjags::run.jags(
      model = system.file("jags", "model.jags", package = "surveyvoi"),
      data = jags_data, monitor = c("coefs"),
      n.chains = 1, sample = n_iter - n_burnin, burnin = n_burnin,
      thin = n_thin, adapt = n_adapt)
  })
}

#' @noRd
hglm_predict <- function(model, x) {
  c(stats::plogis(x %*% unname(model$summaries[, "Median"])))
}
