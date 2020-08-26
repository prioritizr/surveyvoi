#' Approximate expected value of the decision given survey information
#'
#' Calculate the \emph{expected value of the conservation management decision
#' given survey information}. This metric describes the value of the management
#' decision that is expected when the decision maker conducts a surveys a
#' set of sites to inform the decision.
#'
#' @inheritParams evdsi
#'
#' @param n_approx_replicates \code{integer} number of replicates to use for
#'   approximating the expected value calculations. Defaults to 100.
#'
#' @param n_approx_outcomes_per_replicate \code{integer} number of outcomes to
#'   use per replicate for approximation calculations. Defaults to 10000.
#'
#' @param method_approx_outcomes \code{character} name of method that is
#'   used to sample outcomes for approximating the expected value
#'   calculations. Available options are:
#'   \code{"uniform_with_replacement"}, \code{"uniform_without_replacement"},
#'   \code{"weighted_with_replacement"}, \code{"weighted_without_replacement"}.
#'   Uniform sampling methods have an equal chance of returning each
#'   state, and weighted sampling methods are more likely to return
#'   outcomes with a higher prior probability of occurring.
#'   Defaults to \code{"weighted_without_replacement"}.
#'
#' @param seed \code{integer} state of the random number generator for
#'  partitioning data into folds cross-validation and fitting \pkg{xgboost}
#'  models. It is also used for generating outcomes.
#'  This parameter must remain the same to compare results from different
#'  functions using the approximation methods. Defaults to 500.
#;
#' @details This function uses approximation methods to estimate the
#'   expected value calculations. As such, you will need to ensure that
#'   the same seed is used when comparing results to other functions that
#'   use approximation methods. Additionally, the accuracy of these
#'   calculations depend on the arguments to
#'   \code{n_approx_replicates} and \code{n_approx_outcomes_per_replicate}, and
#'   so you may need to increase these parameters for large problems.
#'
#' @return \code{numeric} vector containing the expected values for each
#' replicate.
#'
#' @inherit evdci seealso
#'
#' @examples
#' # set seeds for reproducibility
#' library(RandomFields)
#' set.seed(123)
#' RFoptions(seed = 123)
#'
#' # simulate data
#' site_data <- simulate_site_data(n_sites = 30, n_features = 2, prop = 0.1)
#' feature_data <- simulate_feature_data(n_features = 2, prop = 1)
#' feature_data$target <- c(15, 15)
#'
#' # preview simulated data
#' print(site_data)
#' print(feature_data)
#'
#' # set total budget for managing sites for conservation
#'  # (i.e. 50% of the cost of managing all sites)
#' total_budget <- sum(site_data$management_cost) * 0.5
#'
#' # create a survey scheme that samples the first two sites that
#' # are missing data
#' site_data$survey_site <- FALSE
#' site_data$survey_site[which(is.na(site_data$f1))[1:2]] <- TRUE
#'
#' # define xgboost tuning parameters
#' xgb_parameters <- list(eta = seq(0.1, 0.5, length.out = 3),
#'                        lambda = 10 ^ seq(-1.0, 0.0, length.out = 3),
#'                        objective = "binary:logistic")
#'
#' # calculate expected value of management decision given the survey
#' # information using exact method
#' approx_ev_survey <- approx_evdsi(
#'   site_data, feature_data, c("f1", "f2"), c("p1", "p2"),
#'   c("e1", "e2", "e3"), "management_cost", "survey_site",
#'   "survey_cost", "survey", "survey_sensitivity", "survey_specificity",
#'   "model_sensitivity", "model_specificity",
#'   "target", total_budget, xgb_parameters)
#'
#' # print exact value
#' print(approx_ev_survey)
#'
#' @export
approx_evdsi <- function(
  site_data, feature_data,
  site_detection_columns, site_n_surveys_columns, site_probability_columns,
  site_env_vars_columns,
  site_management_cost_column,
  site_survey_scheme_column,
  site_survey_cost_column,
  feature_survey_column,
  feature_survey_sensitivity_column,
  feature_survey_specificity_column,
  feature_model_sensitivity_column,
  feature_model_specificity_column,
  feature_target_column,
  total_budget,
  xgb_tuning_parameters,
  xgb_early_stopping_rounds = rep(10, length(site_detection_columns)),
  xgb_n_rounds = rep(100, length(site_detection_columns)),
  xgb_n_folds = rep(5, length(site_detection_columns)),
  site_management_locked_in_column = NULL,
  site_management_locked_out_column = NULL,
  prior_matrix = NULL,
  n_approx_replicates = 100,
  n_approx_outcomes_per_replicate = 10000,
  method_approx_outcomes = "weighted_without_replacement",
  seed = 500) {
  # assert arguments are valid
  assertthat::assert_that(
    ## site_data
    inherits(site_data, "sf"), ncol(site_data) > 0,
    nrow(site_data) > 0,
    ## feature_data
    inherits(feature_data, "data.frame"), ncol(feature_data) > 0,
    nrow(feature_data) > 0,
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
    ## site_probability_columns
    is.character(site_probability_columns),
    identical(nrow(feature_data), length(site_probability_columns)),
    assertthat::noNA(site_probability_columns),
    all(assertthat::has_name(site_data, site_probability_columns)),
    ## site_env_vars_columns
    is.character(site_env_vars_columns),
    assertthat::noNA(site_env_vars_columns),
    all(assertthat::has_name(site_data, site_env_vars_columns)),
    ## site_management_cost_column
    assertthat::is.string(site_management_cost_column),
    all(assertthat::has_name(site_data, site_management_cost_column)),
    is.numeric(site_data[[site_management_cost_column]]),
    assertthat::noNA(site_data[[site_management_cost_column]]),
    ## site_survey_scheme_column
    assertthat::is.string(site_survey_scheme_column),
    all(assertthat::has_name(site_data, site_survey_scheme_column)),
    is.logical(site_data[[site_survey_scheme_column]]),
    assertthat::noNA(site_data[[site_survey_scheme_column]]),
    sum(site_data[[site_survey_scheme_column]]) >= 1,
    ## site_survey_cost_column
    assertthat::is.string(site_survey_cost_column),
    all(assertthat::has_name(site_data, site_survey_cost_column)),
    is.numeric(site_data[[site_survey_cost_column]]),
    assertthat::noNA(site_data[[site_survey_cost_column]]),
    ## feature_survey_column
    assertthat::is.string(feature_survey_column),
    all(assertthat::has_name(feature_data, feature_survey_column)),
    is.logical(feature_data[[feature_survey_column]]),
    assertthat::noNA(feature_data[[feature_survey_column]]),
    sum(feature_data[[feature_survey_column]]) >= 1,
    ## feature_survey_sensitivity_column
    assertthat::is.string(feature_survey_sensitivity_column),
    all(assertthat::has_name(feature_data, feature_survey_sensitivity_column)),
    is.numeric(feature_data[[feature_survey_sensitivity_column]]),
    assertthat::noNA(
      feature_data[[feature_survey_sensitivity_column]]),
    all(feature_data[[feature_survey_sensitivity_column]] >= 0),
    all(feature_data[[feature_survey_sensitivity_column]] <= 1),
    ## feature_survey_specificity_column
    assertthat::is.string(feature_survey_specificity_column),
    all(assertthat::has_name(feature_data, feature_survey_specificity_column)),
    is.numeric(feature_data[[feature_survey_specificity_column]]),
    assertthat::noNA(feature_data[[feature_survey_specificity_column]]),
    all(feature_data[[feature_survey_specificity_column]] >= 0),
    all(feature_data[[feature_survey_specificity_column]] <= 1),
    ## feature_model_sensitivity_column
    assertthat::is.string(feature_model_sensitivity_column),
    all(assertthat::has_name(feature_data, feature_model_sensitivity_column)),
    is.numeric(feature_data[[feature_model_sensitivity_column]]),
    assertthat::noNA(feature_data[[feature_model_sensitivity_column]]),
    all(feature_data[[feature_model_sensitivity_column]] >= 0),
    all(feature_data[[feature_model_sensitivity_column]] <= 1),
    ## feature_model_specificity_column
    assertthat::is.string(feature_model_specificity_column),
    all(assertthat::has_name(feature_data, feature_model_specificity_column)),
    is.numeric(feature_data[[feature_model_specificity_column]]),
    assertthat::noNA(feature_data[[feature_model_specificity_column]]),
    all(feature_data[[feature_model_specificity_column]] >= 0),
    all(feature_data[[feature_model_specificity_column]] <= 1),
    ## feature_target_column
    assertthat::is.string(feature_target_column),
    all(assertthat::has_name(feature_data, feature_target_column)),
    is.numeric(feature_data[[feature_target_column]]),
    assertthat::noNA(feature_data[[feature_target_column]]),
    all(feature_data[[feature_target_column]] >= 0),
    ## total_budget
    assertthat::is.number(total_budget), assertthat::noNA(total_budget),
    isTRUE(total_budget > 0),
    ## xgb_tuning_parameters
    is.list(xgb_tuning_parameters),
    ## xgb_n_folds
    is.numeric(xgb_n_folds), assertthat::noNA(xgb_n_folds),
    length(xgb_n_folds) == length(site_detection_columns),
    is.numeric(xgb_n_rounds), assertthat::noNA(xgb_n_rounds),
    length(xgb_n_rounds) == length(site_detection_columns),
    all(xgb_n_rounds > 0),
    ## xgb_early_stopping_rounds
    is.numeric(xgb_early_stopping_rounds),
    assertthat::noNA(xgb_early_stopping_rounds),
    length(xgb_early_stopping_rounds) == length(site_detection_columns),
    all(xgb_early_stopping_rounds > 0),
    ## prior_matrix
    inherits(prior_matrix, c("matrix", "NULL")),
    ## n_approx_replicates
    assertthat::is.count(n_approx_replicates),
    assertthat::noNA(n_approx_replicates),
    ## n_approx_outcomes_per_replicate
    assertthat::is.count(n_approx_outcomes_per_replicate),
    assertthat::noNA(n_approx_outcomes_per_replicate),
    ## method_approx_outcomes
    assertthat::is.string(method_approx_outcomes),
    assertthat::noNA(method_approx_outcomes),
    isTRUE(method_approx_outcomes %in%
      c("uniform_with_replacement", "uniform_without_replacement",
        "weighted_with_replacement", "weighted_without_replacement")),
    ## seed
    assertthat::is.number(seed))
  ## site_management_locked_in_column
  if (!is.null(site_management_locked_in_column)) {
    assertthat::assert_that(
      assertthat::is.string(site_management_locked_in_column),
      all(assertthat::has_name(site_data, site_management_locked_in_column)),
      is.logical(site_data[[site_management_locked_in_column]]),
      assertthat::noNA(site_data[[site_management_locked_in_column]]))
    assertthat::assert_that(
      sum(site_data[[site_management_locked_in_column]] *
          site_data[[site_management_cost_column]]) <=
      total_budget,
      msg = "cost of managing locked in sites exceeds total budget")
  }
  ## site_management_locked_out_column
  if (!is.null(site_management_locked_out_column)) {
    assertthat::assert_that(
      assertthat::is.string(site_management_locked_out_column),
      all(assertthat::has_name(site_data, site_management_locked_out_column)),
      is.logical(site_data[[site_management_locked_out_column]]),
      assertthat::noNA(site_data[[site_management_locked_out_column]]))
    if (all(site_data[[site_management_locked_out_column]]))
      warning("all sites locked out")
  }
  ## validate locked arguments if some locked in and some locked out
  if (!is.null(site_management_locked_in_column) &&
      !is.null(site_management_locked_out_column)) {
    assertthat::assert_that(
      all(site_data[[site_management_locked_in_column]] +
          site_data[[site_management_locked_out_column]] <= 1),
      msg = "at least one planning unit is locked in and locked out")
  }
  ## n_approx_outcomes_per_replicate
  if ((nrow(site_data) * nrow(feature_data)) < 50)
    assertthat::assert_that(
    isTRUE(n_approx_outcomes_per_replicate <=
           n_states(nrow(site_data), nrow(feature_data))))
  ## validate targets
  validate_target_data(feature_data, feature_target_column)
  ## validate survey data
  validate_site_detection_data(site_data, site_detection_columns)
  validate_site_n_surveys_data(site_data, site_n_surveys_columns)
  ## validate model probability values
  validate_site_probability_data(site_data, site_probability_columns)
  ## validate xgboost tuning parameters
  validate_xgboost_tuning_parameters(xgb_tuning_parameters)
  ## verify targets
  assertthat::assert_that(
    all(feature_data[[feature_target_column]] <= nrow(site_data)))
  if (!is.null(site_management_locked_out_column)) {
    assertthat::assert_that(
      all(feature_data[[feature_target_column]] <=
          sum(!site_data[[site_management_locked_out_column]])))
  }

  # prepare data for analysis
  ## drop spatial information
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)
  ## calculate prior matrix
  if (is.null(prior_matrix)) {
    pij <- prior_probability_matrix(
      site_data, feature_data, site_detection_columns,
      site_n_surveys_columns, site_probability_columns,
      feature_survey_sensitivity_column, feature_survey_specificity_column,
      feature_model_sensitivity_column, feature_model_specificity_column)
  } else {
    validate_prior_data(prior_matrix, nrow(site_data), nrow(feature_data))
    pij <- prior_matrix
  }
  ## prepare locked in data
  if (!is.null(site_management_locked_in_column)) {
    site_management_locked_in <- site_data[[site_management_locked_in_column]]
  } else {
    site_management_locked_in <- rep(FALSE, nrow(site_data))
  }
  ## prepare locked out data
  if (!is.null(site_management_locked_out_column)) {
    site_management_locked_out <- site_data[[site_management_locked_out_column]]
  } else {
    site_management_locked_out <- rep(FALSE, nrow(site_data))
  }
  ## validate that targets are feasible given budget and locked out units
  sorted_costs <- sort(
    site_data[[site_management_cost_column]][!site_management_locked_out])
  sorted_costs <- sorted_costs[
    seq_len(max(feature_data[[feature_target_column]]))]
  assertthat::assert_that(
    sum(sorted_costs) <= total_budget,
    msg = paste("targets cannot be achieved given budget and locked out",
                "planning units"))
  ## prepare parameter combinations for model tuning
  xgb_full_parameters <- do.call(expand.grid, xgb_tuning_parameters)
  attr(xgb_full_parameters, "out.attrs") <- NULL
  xgb_full_parameters$nthread <- "1" # force single thread for reproducibility
  xgb_full_parameters$verbose <- "0" # force quiet
  xgb_full_parameters$seed <- as.character(seed) # set seed
  if (is.null(xgb_full_parameters$objective)) {
    xgb_full_parameters$objective <- "binary:logistic"
    warning(paste("no objective specified for model fitting,",
                  "assuming binary:logistic"))
  }
  ## extract site data
  dij <- t(as.matrix(site_data[, site_detection_columns, drop = FALSE]))
  nij <- t(as.matrix(site_data[, site_n_surveys_columns, drop = FALSE]))
  ejx <- as.matrix(site_data[, site_env_vars_columns])
  ## identify sites that need model predictions for each feature
  pu_model_prediction <- lapply(seq_len(nrow(feature_data)), function(i) {
    which(!site_data[[site_survey_scheme_column]] & (nij[i, ] < 0.5))
  })
  ## folds for training and testing models
  xgb_folds <- lapply(seq_len(nrow(feature_data)), function(i) {
    has_data_idx <- which(
      (nij[i, ] > 0.5) | site_data[[site_survey_scheme_column]])
    create_site_folds(
      dij[i, has_data_idx], nij[i, has_data_idx],
      index = has_data_idx, xgb_n_folds[i], seed = seed)
  })

  # main calculation
  withr::with_seed(seed, {
    out <- rcpp_approx_expected_value_of_decision_given_survey_scheme(
      dij = dij, nij = nij, pij = pij,
      survey_features = feature_data[[feature_survey_column]],
      survey_sensitivity = feature_data[[feature_survey_sensitivity_column]],
      survey_specificity = feature_data[[feature_survey_specificity_column]],
      pu_survey_solution = site_data[[site_survey_scheme_column]],
      pu_model_prediction = pu_model_prediction,
      pu_survey_costs = site_data[[site_survey_cost_column]],
      pu_purchase_costs = site_data[[site_management_cost_column]],
      pu_purchase_locked_in = site_management_locked_in,
      pu_purchase_locked_out = site_management_locked_out,
      pu_env_data = ejx,
      xgb_parameter_names = names(xgb_full_parameters),
      xgb_parameter_values = as.matrix(xgb_full_parameters),
      n_xgb_rounds = xgb_n_rounds,
      n_xgb_early_stopping_rounds = xgb_early_stopping_rounds,
      xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
      xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
      obj_fun_target = round(feature_data[[feature_target_column]]),
      total_budget = total_budget,
      n_approx_replicates = n_approx_replicates,
      n_approx_outcomes_per_replicate = n_approx_outcomes_per_replicate,
      method_approx_outcomes = method_approx_outcomes,
      seed = seed)
  })
  # return result
  out
}
