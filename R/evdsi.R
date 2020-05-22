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
#' @param feature_preweight_column \code{character} name of the column in the
#'   argument to \code{feature_data} that contains the \eqn{preweight}
#'   values used to parametrize  the conservation benefit of managing of each
#'   feature.
#'   This column should have \code{numeric} values that
#'   are equal to or greater than zero. No missing (\code{NA}) values are
#'   permitted in this column.
#'
#' @param feature_postweight_column \code{character} name of the column in the
#'   argument to \code{feature_data} that contains the \eqn{postweight}
#'   values used to parametrize  the conservation benefit of managing of each
#'   feature.
#'   This column should have \code{numeric} values that
#'   are equal to or greater than zero. No missing (\code{NA}) values are
#'   permitted in this column.
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
#' @param prior_matrix \code{numeric} \code{matrix} containing
#'  the prior probability of each feature occupying each site.
#'  Rows correspond to features, and columns correspond to sites.
#'  Defaults to \code{NULL} such that prior data is calculated automatically
#'  using \code{\link{prior_probability_matrix}}.
#'
#' @param xgb_parameters \code{list} of \code{list} objects
#'   containing the parameters for fitting models for each
#'   feature. See documentation for the \code{params} argument in
#'   \code{\link[xgboost]{xgb.train}} for available parameters. Ideally,
#'   these parameters would be determined using the
#'   \code{\link{fit_occupancy_models}} function. Note that arguments must
#'   have \code{"nrounds"}, \code{"objective"}, \code{"scale_pos_weight"}
#'   elements (see example below).
#'
#' @param xgb_n_folds \code{integer} vector containing the number of
#'   k-fold cross-validation folds to use for fitting models and
#'   assessing model performance for each feature. Ideally, the number of folds
#'   should be exactly the same as the number used for tuning the
#'   model parameters (i.e. same parameter to the \code{n_folds}
#'   argument in \link{fit_occupancy_models} when generating parameters
#'  for \code{xgb_parameters}).
#'
#' @param optimality_gap \code{numeric} relative optimality gap for generating
#'   conservation prioritizations. A value of zero indicates that
#'   prioritizations must be solved to optimality. A value of 0.1 indicates
#'   prioritizations must be within 10\% of optimality. Defaults to 0.
#'
#' @param seed \code{integer} state of the random number generator for
#'  partitioning data into folds cross-validation and fitting \pkg{xgboost}
#'  models. This parameter must remain the same to compare results for
#'  different survey schemes. Defaults to 500.
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
#' feature_data <- simulate_feature_data(n_features = 2, prop = 1)
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
#'   "preweight", "postweight", "target",
#'   total_budget, xgb_parameters)
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
  feature_preweight_column,
  feature_postweight_column,
  feature_target_column,
  total_budget,
  xgb_parameters,
  site_management_locked_in_column = NULL,
  prior_matrix = NULL,
  optimality_gap = 0,
  site_weight_columns = NULL,
  xgb_n_folds = rep(5, nrow(feature_data)),
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
    ## feature_preweight_column
    assertthat::is.string(feature_preweight_column),
    all(assertthat::has_name(feature_data, feature_preweight_column)),
    is.numeric(feature_data[[feature_preweight_column]]),
    assertthat::noNA(feature_data[[feature_preweight_column]]),
    all(feature_data[[feature_preweight_column]] >= 0),
    ## feature_postweight_column
    assertthat::is.string(feature_postweight_column),
    all(assertthat::has_name(feature_data, feature_postweight_column)),
    is.numeric(feature_data[[feature_postweight_column]]),
    assertthat::noNA(feature_data[[feature_postweight_column]]),
    all(feature_data[[feature_postweight_column]] >= 0),
    ## feature_target_column
    assertthat::is.string(feature_target_column),
    all(assertthat::has_name(feature_data, feature_target_column)),
    is.numeric(feature_data[[feature_target_column]]),
    assertthat::noNA(feature_data[[feature_target_column]]),
    all(feature_data[[feature_target_column]] >= 0),
    ## total_budget
    assertthat::is.number(total_budget), assertthat::noNA(total_budget),
    isTRUE(total_budget > 0),
    ## xgb_parameters
    is.list(xgb_parameters),
    identical(length(xgb_parameters), nrow(feature_data)),
    ## xgb_n_folds
    is.numeric(xgb_n_folds),
    all(xgb_n_folds > 0), identical(length(xgb_n_folds), nrow(feature_data)),
    assertthat::noNA(xgb_n_folds),
    ## prior_matrix
    inherits(prior_matrix, c("matrix", "NULL")),
    ## optimality_gap
    assertthat::is.number(optimality_gap),
    assertthat::noNA(optimality_gap),
    isTRUE(optimality_gap >= 0),
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
  ## validate rij values
  validate_site_occupancy_data(site_data, site_occupancy_columns)
  ## validate pij values
  validate_site_prior_data(site_data, site_probability_columns)
  ## validate wij values
  if (!is.null(site_weight_columns))
    validate_site_weight_data(site_data, site_occupancy_columns,
      site_weight_columns)
  ## validate xgboost parameters
  validate_xgboost_parameters(xgb_parameters)

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
  ## xgb_nrounds
  xgb_nrounds <- vapply(xgb_parameters, `[[`,  FUN.VALUE = numeric(1),
                        "nrounds")
  ## format xgb_parameters
  xgb_parameters <- lapply(xgb_parameters, function(x) {
    out <- x[names(x) != "nrounds"]
    out <- lapply(out, as.character)
    out$nthread <- "1" # force single thread for reproducibility
    out$seed <- as.character(seed)
    out
  })
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
        create_folds(unname(rij[i, pu_train_idx]), xgb_n_folds[i],
                     index = pu_train_idx,
                     na.fail = FALSE, seed = seed)
      })
  })
  ## extract site weight data
  if (!is.null(site_weight_columns)) {
    wij <- t(as.matrix(site_data[, site_weight_columns]))
  } else {
    wij <- matrix(1, ncol = ncol(rij), nrow = nrow(rij))
  }
  ## extract environmental data
  ejx <- as.matrix(site_data[, site_env_vars_columns])
  ## prepare rij matrix for Rcpp
  rij[is.na(rij)] <- -1

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
      pu_env_data = ejx,
      xgb_parameters = xgb_parameters,
      xgb_train_folds = lapply(xgb_folds, `[[`, "train"),
      xgb_test_folds = lapply(xgb_folds, `[[`, "test"),
      n_xgb_nrounds = xgb_nrounds,
      obj_fun_preweight = feature_data[[feature_preweight_column]],
      obj_fun_postweight = feature_data[[feature_postweight_column]],
      obj_fun_target = feature_data[[feature_target_column]],
      total_budget = total_budget,
      optim_gap = optimality_gap)
  })
  # return result
  out
}
