#' Expected value of the decision given survey information
#'
#' Calculate the \emph{expected value of the conservation management decision
#' given survey information}. This metric describes the value of the management
#' decision that is expected when the decision maker conducts a surveys a
#' set of sites to inform the decision.
#'
#' @param site_data \code{\link[sf]{sf}} object with site data.
#'
#' @param feature_data \code{\link[base]{data.frame}} object with feature data.
#'
#' @param site_occupancy_columns \code{character} names of \code{numeric}
#'   columns in the argument to \code{site_data} that contain detections
#'   and non-detections of features at each site.
#'   Each column should correspond to a different feature, and contain
#'   binary detection/non-detection value (zeros or ones) indicating if the
#'   feature was detected in a previous survey or not. If a site has not
#'   been surveyed before, then missing (\code{NA}) values should be used.
#'   Additionally, if a feature was not looked for when surveying a specific
#'   site, then missing (\code{NA}) value should be used.
#'
#' @param site_probability_columns \code{character} names of \code{numeric}
#'   columns in the argument to \code{site_data} that contain modelled
#'   probabilities of occupancy for each feature in each site.
#'   Each column should correspond to a different feature, and contain
#'   probability data (values between zero and one). No missing (\code{NA})
#'   values are permitted in these columns.
#'
#' @param site_survey_scheme_column \code{character} name of \code{logical}
#'  (\code{TRUE} / \code{FALSE}) column in the argument to \code{site_data}
#'  that indicates which sites are selected in the scheme or not.
#'  No missing \code{NA} values are permitted. Additionally, only sites
#'  that are missing data can be selected or surveying (as per the
#'  argument to \code{site_occupancy_columns}).
#'
#' @param feature_survey_column \code{character} name of the column in the
#'   argument to \code{feature_data} that contains \code{logical} (\code{TRUE} /
#'   \code{FALSE}) values indicating if the feature will be surveyed in
#'   the planned surveys or not. Note that considering additional features will
#'   rapidly increase computational burden, and so it is only recommended to
#'   consider features that are of specific conservation interest.
#'   No missing (\code{NA}) values are permitted in this column.
#'
#' @param site_survey_cost_column \code{character} name of column in the
#'   argument to  \code{site_data} that contains costs for surveying each
#'   site. This column should have \code{numeric} values that are equal to
#'   or greater than zero. No missing (\code{NA}) values are permitted in this
#'   column.
#'
#' @param site_env_vars_columns \code{character} names of columns in the
#'   argument to  \code{site_data} that contain environmental information
#'   for fitting updated occupancy models based on possible survey outcomes.
#'   Each column should correspond to a different environmental variable,
#'   and contain \code{numeric}, \code{factor}, or \code{character} data.
#'   No missing (\code{NA}) values are permitted in these columns.
#'
#' @param site_weight_columns \code{character} name of columns in
#'  \code{site_data} containing weights for model fitting for each
#'  feature. These columns must contain \code{numeric} values greater
#'  than or equal to zero. No missing (\code{NA}) values are
#'  permitted. Defaults to \code{NULL} such that all data are given
#'  equal weight when fitting models.
#'
#' @param site_management_cost_column \code{character} name of column in the
#'   argument to \code{site_data} that contains costs for managing each
#'   site for conservation. This column should have \code{numeric} values that
#'   are equal to or greater than zero. No missing (\code{NA}) values are
#'   permitted in this column.
#'
#' @param feature_survey_sensitivity_column \code{character} name of the
#'   column in the argument to \code{feature_data} that contains
#'   probability of future surveys correctly detecting a presence of each
#'   feature in a given site (i.e. the sensitivity of the survey methodology).
#'   This column should have \code{numeric} values that are between zero and
#'   one. No missing (\code{NA}) values are permitted in this column.
#'
#' @param feature_survey_specificity_column \code{character} name of the
#'   column in the argument to \code{feature_data} that contains
#'   probability of future surveys correctly detecting an absence of each
#'   feature in a given site (i.e. the specificity of the survey methodology).
#'   This column should have \code{numeric} values that are between zero and
#'   one. No missing (\code{NA}) values are permitted in this column.
#'
#' @param feature_model_sensitivity_column \code{character} name of the
#'   column in the argument to \code{feature_data} that contains
#'   probability of the initial models correctly predicting a presence of each
#'   feature in a given site (i.e. the sensitivity of the models).
#'   This column should have \code{numeric} values that are between zero and
#'   one. No missing (\code{NA}) values are permitted in this column.
#'   This should ideally be calculated using \code{\link{fit_occupancy_models}}.
#'
#' @param feature_model_specificity_column \code{character} name of the
#'   column in the argument to \code{feature_data} that contains
#'   probability of the initial models correctly predicting an absence of each
#'   feature in a given site (i.e. the specificity of the models).
#'   This column should have \code{numeric} values that are between zero and
#'   one. No missing (\code{NA}) values are permitted in this column.
#'   This should ideally be calculated using \code{\link{fit_occupancy_models}}.
#'
#' @param feature_target_column \code{character} name of the column in the
#'   argument to \code{feature_data} that contains the \eqn{target}
#'   values used to parametrize the conservation benefit of managing of each
#'   feature.
#'   This column should have \code{numeric} values that
#'   are equal to or greater than zero. No missing (\code{NA}) values are
#'   permitted in this column.
#'
#' @param total_budget \code{numeric} maximum expenditure permitted
#'   for conducting surveys and managing sites for conservation.
#'
#' @param site_management_locked_in_column \code{character} name of the column
#'   in the argument to \code{site_data} that contains \code{logical}
#'   (\code{TRUE} / \code{FALSE}) values indicating which sites should
#'   be locked in for (\code{TRUE}) being managed for conservation or
#'   (\code{FALSE}) not. No missing (\code{NA}) values are permitted in this
#'   column. This is useful if some sites have already been earmarked for
#'   conservation, or if some sites are already being managed for conservation.
#'   Defaults to \code{NULL} such that no sites are locked in.
#'
#' @param site_management_locked_out_column \code{character} name of the column
#'   in the argument to \code{site_data} that contains \code{logical}
#'   (\code{TRUE} / \code{FALSE}) values indicating which sites should
#'   be locked out for (\code{TRUE}) being managed for conservation or
#'   (\code{FALSE}) not. No missing (\code{NA}) values are permitted in this
#'   column. This is useful if some sites could potentially be surveyed
#'   to improve model predictions even if they cannot be managed for
#'   conservation. Defaults to \code{NULL} such that no sites are locked out.
#'
#' @param prior_matrix \code{numeric} \code{matrix} containing
#'  the prior probability of each feature occupying each site.
#'  Rows correspond to features, and columns correspond to sites.
#'  Defaults to \code{NULL} such that prior data is calculated automatically
#'  using \code{\link{prior_probability_matrix}}.
#'
#' @param seed \code{integer} state of the random number generator for
#'  partitioning data into folds cross-validation and fitting \pkg{xgboost}
#'  models. This parameter must remain the same to compare results for
#'  different survey schemes. Defaults to 500.
#'
#' @inheritParams fit_occupancy_models
#'
#' @details This function calculates the expected value and does not
#'  use approximation methods. As such, this function can only be applied
#'  to very small problems.
#'
#' @return \code{numeric} value.

#' @seealso \code{\link{prior_probability_matrix}}.
#'
#' @examples
#' # set seeds for reproducibility
#' library(RandomFields)
#' set.seed(123)
#' RFoptions(seed = 123)
#'
#' # simulate data
#' site_data <- simulate_site_data(n_sites = 15, n_features = 2, prop = 0.5)
#' feature_data <- simulate_feature_data(n_features = 2, prop = 0.5)
#' feature_data$target <- c(3, 3)
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
#' # these should ideally be determined using fit_occupancy_models
#' xgb_parameters <-
#'  list(list(objective = "binary:logistic", nrounds = 8,
#'            scale_pos_weight = 1, eta = 0.1))[rep(1, 2)]
#'
#' # calculate expected value of management decision given the survey
#' # information using exact method
#' ev_survey <- evdsi(
#'   site_data, feature_data, c("f1", "f2"), c("p1", "p2"),
#'   c("e1", "e2", "e3"), "management_cost", "survey_site",
#'   "survey_cost", "survey", "survey_sensitivity", "survey_specificity",
#'   "model_sensitivity", "model_specificity",
#'   "target", total_budget, xgb_parameters)
#'
#' # print exact value
#' print(ev_survey)
#'
#' @export
evdsi <- function(
  site_data,
  feature_data,
  site_occupancy_columns,
  site_probability_columns,
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
  xgb_early_stopping_rounds = rep(100, length(site_occupancy_columns)),
  xgb_n_rounds = rep(1000, length(site_occupancy_columns)),
  xgb_n_folds = rep(5, length(site_occupancy_columns)),
  site_management_locked_in_column = NULL,
  site_management_locked_out_column = NULL,
  prior_matrix = NULL,
  site_weight_columns = NULL,
  seed = 500) {
  # assert arguments are valid
  assertthat::assert_that(
    ## site_data
    inherits(site_data, "sf"), ncol(site_data) > 0,
    nrow(site_data) > 0,
    ## feature_data
    inherits(feature_data, "data.frame"), ncol(feature_data) > 0,
    nrow(feature_data) > 0,
    ## site_occupancy_columns
    is.character(site_occupancy_columns),
    identical(nrow(feature_data), length(site_occupancy_columns)),
    assertthat::noNA(site_occupancy_columns),
    all(assertthat::has_name(site_data, site_occupancy_columns)),
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
    length(xgb_n_folds) == length(site_occupancy_columns),
    is.numeric(xgb_n_rounds), assertthat::noNA(xgb_n_rounds),
    length(xgb_n_rounds) == length(site_occupancy_columns),
    all(xgb_n_rounds > 0),
    ## xgb_early_stopping_rounds
    is.numeric(xgb_early_stopping_rounds),
    assertthat::noNA(xgb_early_stopping_rounds),
    length(xgb_early_stopping_rounds) == length(site_occupancy_columns),
    all(xgb_early_stopping_rounds > 0),
    ## prior_matrix
    inherits(prior_matrix, c("matrix", "NULL")),
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
  ## validate targets
  validate_target_data(feature_data, feature_target_column)
  ## validate rij values
  validate_site_occupancy_data(site_data, site_occupancy_columns)
  ## validate pij values
  validate_site_prior_data(site_data, site_probability_columns)
  ## validate wij values
  if (!is.null(site_weight_columns))
    validate_site_weight_data(site_data, site_occupancy_columns,
      site_weight_columns)
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
      site_data, feature_data, site_occupancy_columns, site_probability_columns,
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
  ## extract site occupancy data
  rij <- t(as.matrix(site_data[, site_occupancy_columns, drop = FALSE]))
  ## identify sites that need model predictions for each feature
  pu_model_prediction <- lapply(seq_len(nrow(feature_data)), function(i) {
    which(!site_data[[site_survey_scheme_column]] & is.na(rij[i, ]))
  })
  ## folds for training and testing models
  xgb_folds <- lapply(seq_len(nrow(feature_data)),
    function(i) {
      pu_train_idx <-
        which(site_data[[site_survey_scheme_column]] | !is.na(rij[i, ]))
      withr::with_seed(seed, {
        create_folds(unname(rij[i, pu_train_idx]),
                     n = xgb_n_folds[i],
                     index = pu_train_idx,
                     na.fail = FALSE, seed = seed)
      })
  })
  ## extract environmental data
  ejx <- as.matrix(site_data[, site_env_vars_columns])
  ## prepare rij matrix for Rcpp
  rij[is.na(rij)] <- -1

  # prepare model weights
  ## initialize weights
  if (!is.null(site_weight_columns)) {
    ## extract user-specified weights
    wij <- t(as.matrix(site_data[, site_weight_columns]))
  } else {
    ## set weights based on if data are missing or not
    wij <- t(as.matrix(site_data[, site_occupancy_columns]))
    wij[] <- as.numeric(!is.na(wij))
  }

  # main calculation
  withr::with_seed(seed, {
    out <- rcpp_expected_value_of_decision_given_survey_scheme(
      rij = rij, pij = pij, wij = wij,
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
      total_budget = total_budget)
  })
  # return result
  out
}
