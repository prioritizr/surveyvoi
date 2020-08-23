#' Expected value of the decision given current information
#'
#' Calculate the \emph{expected value of the conservation management decision
#' given current information}. This metric describes the value of the management
#' decision that is expected when the decision maker is limited to
#' existing biodiversity data (i.e. survey data and environmental niche models).
#'
#' @inheritParams evdsi
#'
#' @inherit evdsi details return seealso
#'
#' @examples
#' # set seeds for reproducibility
#' library(RandomFields)
#' set.seed(123)
#' RFoptions(seed = 123)
#'
#' # simulate data
#' site_data <- simulate_site_data(n_sites = 5, n_features = 2, prop = 0.5)
#' feature_data <- simulate_feature_data(n_features = 2, prop = 1)
#' feature_data$target <- c(1, 1)
#'
#' # preview simulated data
#' print(site_data)
#' print(feature_data)
#'
#' # set total budget for managing sites for conservation
#  # (i.e. 50% of the cost of managing all sites)
#' total_budget <- sum(site_data$management_cost) * 0.5
#'
#' # calculate expected value of management decision given current information
#' # using exact method
#' ev_current <- evdci(
#'   site_data, feature_data, c("f1", "f2"), c("p1", "p2"),
#'   "management_cost", "survey_sensitivity", "survey_specificity",
#'   "model_sensitivity", "model_specificity",
#'   "target", total_budget)
#'
#' # print exact value
#' print(ev_current)
#'
#' @export
evdci <- function(
  site_data, feature_data,
  site_detection_columns, site_n_surveys_columns, site_probability_columns,
  site_management_cost_column,
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
  ## validate survey data
  validate_site_detection_data(site_data, site_detection_columns)
  validate_site_n_surveys_data(site_data, site_n_surveys_columns)
  ## validate targets
  validate_target_data(feature_data, feature_target_column)
  ## validate pij values
  validate_site_prior_data(site_data, site_probability_columns)
  ## verify targets
  assertthat::assert_that(
    all(feature_data[[feature_target_column]] <= nrow(site_data)))
  if (!is.null(site_management_locked_out_column)) {
    assertthat::assert_that(
      all(feature_data[[feature_target_column]] <=
          sum(!site_data[[site_management_locked_out_column]])))
  }

  # drop spatial information
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)

  # calculate prior matrix
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
  out <- rcpp_expected_value_of_decision_given_current_info(
    pij = pij,
    pu_costs = site_data[[site_management_cost_column]],
    pu_locked_in = site_management_locked_in,
    pu_locked_out = site_management_locked_out,
    target = round(feature_data[[feature_target_column]]),
    budget = total_budget)

  # return result
  out
}
