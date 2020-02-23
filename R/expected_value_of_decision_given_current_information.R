#' Expected value of decision given current information
#'
#' Calculate the expected value of the management decision given current
#' information.
#'
#' @param site_data \code{\link[sf]{sf}} object with site data.
#'
#' @param feature_data \code{\link[base]{data.frame}} object with feature data.
#'
#' @param site_occupancy_columns \code{character} names of \code{numeric}
#'   columns in the
#'   argument to \code{site_data} that contain presence/absence data.
#'   Each column should correspond to a different feature, and contain
#'   binary presence/absence data (zeros or ones) indicating if the
#'   feature was detected in a previous survey or not. If a site has not
#'   been surveyed before, then missing (\code{NA}) values should be used.
#'
#' @param site_probability_columns \code{character} names of \code{numeric}
#'   columns in the argument to \code{site_data} that contain modelled
#'   probabilities of occupancy for each feature in each site.
#'   Each column should correspond to a different feature, and contain
#'   probability data (values between zero and one). No missing (\code{NA})
#'   values are permitted in these columns.
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
#' @param feature_alpha_column \code{character} name of the column in the
#'   argument to \code{feature_data} that contains the \eqn{\alpha}
#'   values used to parametrize
#'   the conservation benefit of managing of each feature.
#'   This column should have \code{numeric} values that
#'   are equal to or greater than zero. No missing (\code{NA}) values are
#'   permitted in this column.
#'
#' @param feature_gamma_column \code{character} name of the column in the
#'   argument to \code{feature_data} that contains the \eqn{\gamma}
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
#' @param n_approx_obj_fun_points \code{integer} number of points to use
#'  for approximating the piecewise-linear components of the objective
#'  function. Greater values result in more precise calculations, but
#'  also incur more computational costs. Defaults to 1000.
#'
#' @param n_approx_states \code{integer} number of states to use for
#'   approximating the expected value of a given management action.
#'   Valid arguments include an integer value greater than zero (e.g. 10000), or
#'   \code{NULL}. If \code{NULL} is supplied as an argument, then the expected
#'   value of a given management action is calculated precisely and not using
#'   approximation methods. Defaults to \code{NULL}.

#' @param optimality_gap \code{numeric} relative optimality gap for generating
#'   conservation prioritizations. A value of zero indicates that
#'   prioritizations must be solved to optimality. A value of 0.1 indicates
#'   prioritizations must be within 10\% of optimality. Defaults to 0.
#'
#' @param seed \code{integer} random seed used for calculations. This
#'   parameter only has an influence when approximation methods are used
#'   for calculating the expected value of a given management action.
#'   Defaults to 500.
#'
#' @details
#' TODO.
#'
#' @return \code{numeric} value.
#'
#' @examples
#' # TODO
#'
#' @export
expected_value_of_decision_given_current_information <- function(
  site_data,
  feature_data,
  site_occupancy_columns,
  site_probability_columns,
  site_management_cost_column,
  feature_survey_sensitivity_column,
  feature_survey_specificity_column,
  feature_model_sensitivity_column,
  feature_alpha_column,
  feature_gamma_column,
  total_budget,
  site_management_locked_in_column = NULL,
  prior_matrix = NULL,
  n_approx_obj_fun_points = 1000,
  n_approx_states = NULL,
  optimality_gap = 0,
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
    ## site_management_cost_column
    assertthat::is.string(site_management_cost_column),
    all(assertthat::has_name(site_data, site_management_cost_column)),
    is.numeric(site_data[[site_management_cost_column]]),
    assertthat::noNA(site_data[[site_management_cost_column]]),
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
    ## feature_alpha_column
    assertthat::is.string(feature_alpha_column),
    all(assertthat::has_name(feature_data, feature_alpha_column)),
    is.numeric(feature_data[[feature_alpha_column]]),
    assertthat::noNA(feature_data[[feature_alpha_column]]),
    all(feature_data[[feature_alpha_column]] >= 0),
    ## feature_gamma_column
    assertthat::is.string(feature_gamma_column),
    all(assertthat::has_name(feature_data, feature_gamma_column)),
    is.numeric(feature_data[[feature_gamma_column]]),
    assertthat::noNA(feature_data[[feature_gamma_column]]),
    all(feature_data[[feature_gamma_column]] >= 0),
    ## total_budget
    assertthat::is.number(total_budget), assertthat::noNA(total_budget),
    isTRUE(total_budget >= 0),
    ## prior_matrix
    inherits(prior_matrix, c("matrix", "NULL")),
    ## n_approx_obj_fun_points
    assertthat::is.number(n_approx_obj_fun_points),
    assertthat::noNA(n_approx_obj_fun_points),
    isTRUE(n_approx_obj_fun_points > 0),
    ## n_approx_states
    inherits(n_approx_states, c("numeric", "NULL")),
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
  ## validate n_approx_states
  if (!is.null(n_approx_states)) {
    assertthat::assert_that(
      assertthat::is.count(n_approx_states), assertthat::noNA(n_approx_states),
      isTRUE(n_approx_states <=
             rcpp_n_states(nrow(site_data) * nrow(feature_data))))
  }
  ## validate rij values
  validate_site_occupancy_data(site_data, site_occupancy_columns)
  ## validate pij values
  validate_site_prior_data(site_data, site_probability_columns)

  # drop spatial information
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)

  # calculate prior matrix
  if (is.null(prior_matrix)) {
    pij <- prior_probability_matrix(
      site_data, feature_data, site_occupancy_columns, site_probability_columns,
      feature_survey_sensitivity_column, feature_survey_specificity_column,
      feature_model_sensitivity_column)
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

  # set the seed
  rng_state <- .Random.seed
  set.seed(seed)

  # main calculation
  if (is.null(n_approx_states)) {
    out <- rcpp_expected_value_of_decision_given_current_info(
      pij = pij,
      pu_costs = site_data[[site_management_cost_column]],
      pu_locked_in = site_management_locked_in,
      alpha = feature_data[[feature_alpha_column]],
      gamma = feature_data[[feature_gamma_column]],
      n_approx_obj_fun_points = n_approx_obj_fun_points,
      budget = total_budget,
      gap = optimality_gap)
  } else {
    out <- rcpp_approx_expected_value_of_decision_given_current_info_n_states(
      pij = pij,
      pu_costs = site_data[[site_management_cost_column]],
      pu_locked_in = site_management_locked_in,
      alpha = feature_data[[feature_alpha_column]],
      gamma = feature_data[[feature_gamma_column]],
      n_approx_obj_fun_points = n_approx_obj_fun_points,
      budget = total_budget,
      gap = optimality_gap,
      n_approx_states = n_approx_states)
  }

  # restore the previous state
  set.seed(rng_state)

  # return result
  out
}
