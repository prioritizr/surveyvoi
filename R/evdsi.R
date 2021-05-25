#' Expected value of the decision given survey information
#'
#' Calculate the *expected value of the management decision
#' given survey information*. This metric describes the value of the management
#' decision that is expected when the decision maker surveys a
#' set of sites to help inform the decision.
#'
#' @inheritParams fit_xgb_occupancy_models
#'
#' @param site_probability_columns `character` names of `numeric`
#'   columns in the argument to `site_data` that contain modelled
#'   probabilities of occupancy for each feature in each site.
#'   Each column should correspond to a different feature, and contain
#'   probability data (values between zero and one). No missing (`NA`)
#'   values are permitted in these columns.
#'
#' @param site_survey_scheme_column `character` name of `logical`
#'  (`TRUE` / `FALSE`) column in the argument to `site_data`
#'  that indicates which sites are selected in the scheme or not.
#'  No missing `NA` values are permitted. Additionally, only sites
#'  that are missing data can be selected or surveying (as per the
#'  argument to `site_detection_columns`).
#'
#' @param feature_survey_column `character` name of the column in the
#'   argument to `feature_data` that contains `logical` (`TRUE` /
#'   `FALSE`) values indicating if the feature will be surveyed in
#'   the planned surveys or not. Note that considering additional features will
#'   rapidly increase computational burden, and so it is only recommended to
#'   consider features that are of specific conservation interest.
#'   No missing (`NA`) values are permitted in this column.
#'
#' @param site_survey_cost_column `character` name of column in the
#'   argument to  `site_data` that contains costs for surveying each
#'   site. This column should have `numeric` values that are equal to
#'   or greater than zero. No missing (`NA`) values are permitted in this
#'   column.
#'
#' @param site_management_cost_column `character` name of column in the
#'   argument to `site_data` that contains costs for managing each
#'   site for conservation. This column should have `numeric` values that
#'   are equal to or greater than zero. No missing (`NA`) values are
#'   permitted in this column.
#'
#' @param feature_model_sensitivity_column `character` name of the
#'   column in the argument to `feature_data` that contains
#'   probability of the initial models correctly predicting a presence of each
#'   feature in a given site (i.e. the sensitivity of the models).
#'   This column should have `numeric` values that are between zero and
#'   one. No missing (`NA`) values are permitted in this column.
#'   This should ideally be calculated using
#'   [fit_xgb_occupancy_models()] or
#'   [fit_hglm_occupancy_models()].
#'
#' @param feature_model_specificity_column `character` name of the
#'   column in the argument to `feature_data` that contains
#'   probability of the initial models correctly predicting an absence of each
#'   feature in a given site (i.e. the specificity of the models).
#'   This column should have `numeric` values that are between zero and
#'   one. No missing (`NA`) values are permitted in this column.
#'   This should ideally be calculated using
#'   [fit_xgb_occupancy_models()] or
#'   [fit_hglm_occupancy_models()].
#'
#' @param feature_target_column `character` name of the column in the
#'   argument to `feature_data` that contains the \eqn{target}
#'   values used to parametrize the conservation benefit of managing of each
#'   feature.
#'   This column should have `numeric` values that
#'   are equal to or greater than zero. No missing (`NA`) values are
#'   permitted in this column.
#'
#' @param total_budget `numeric` maximum expenditure permitted
#'   for conducting surveys and managing sites for conservation.
#'
#' @param site_management_locked_in_column `character` name of the column
#'   in the argument to `site_data` that contains `logical`
#'   (`TRUE` / `FALSE`) values indicating which sites should
#'   be locked in for (`TRUE`) being managed for conservation or
#'   (`FALSE`) not. No missing (`NA`) values are permitted in this
#'   column. This is useful if some sites have already been earmarked for
#'   conservation, or if some sites are already being managed for conservation.
#'   Defaults to `NULL` such that no sites are locked in.
#'
#' @param site_management_locked_out_column `character` name of the column
#'   in the argument to `site_data` that contains `logical`
#'   (`TRUE` / `FALSE`) values indicating which sites should
#'   be locked out for (`TRUE`) being managed for conservation or
#'   (`FALSE`) not. No missing (`NA`) values are permitted in this
#'   column. This is useful if some sites could potentially be surveyed
#'   to improve model predictions even if they cannot be managed for
#'   conservation. Defaults to `NULL` such that no sites are locked out.
#'
#' @param prior_matrix `numeric` `matrix` containing
#'  the prior probability of each feature occupying each site.
#'  Rows correspond to features, and columns correspond to sites.
#'  Defaults to `NULL` such that prior data is calculated automatically
#'  using [prior_probability_matrix()].
#'
#' @details This function calculates the expected value and does not
#'  use approximation methods. As such, this function can only be applied
#'  to very small problems.
#'
#' @return `numeric` value.
#'
#' @seealso [prior_probability_matrix()].
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
#' feature_data$target <- c(10, 10)
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
#' site_data$survey_site[which(site_data$n1 < 0.5)[1:2]] <- TRUE
#'
#' # calculate expected value of management decision given the survey
#' # information using exact method
#' ev_survey <- evdsi(
#'   site_data, feature_data,
#'   c("f1", "f2"), c("n1", "n2"), c("p1", "p2"),
#'   "management_cost", "survey_site",
#'   "survey_cost", "survey", "survey_sensitivity", "survey_specificity",
#'   "model_sensitivity", "model_specificity",
#'   "target", total_budget)
#'
#' # print value
#' print(ev_survey)
#' @export
evdsi <- function(
  site_data, feature_data,
  site_detection_columns, site_n_surveys_columns, site_probability_columns,
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
  site_management_locked_in_column = NULL,
  site_management_locked_out_column = NULL,
  prior_matrix = NULL) {
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
    ## prior_matrix
    inherits(prior_matrix, c("matrix", "NULL")))
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
  ## validate survey data
  validate_site_detection_data(site_data, site_detection_columns)
  validate_site_n_surveys_data(site_data, site_n_surveys_columns)
  ## validate model probability values
  validate_site_probability_data(site_data, site_probability_columns)
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
  # main calculation
  out <- rcpp_expected_value_of_decision_given_survey_scheme(
    pij = pij,
    survey_features = feature_data[[feature_survey_column]],
    survey_sensitivity = feature_data[[feature_survey_sensitivity_column]],
    survey_specificity = feature_data[[feature_survey_specificity_column]],
    pu_survey_solution = site_data[[site_survey_scheme_column]],
    pu_survey_costs = site_data[[site_survey_cost_column]],
    pu_purchase_costs = site_data[[site_management_cost_column]],
    pu_purchase_locked_in = site_management_locked_in,
    pu_purchase_locked_out = site_management_locked_out,
    obj_fun_target = round(feature_data[[feature_target_column]]),
    total_budget = total_budget)
  # return result
  out
}
