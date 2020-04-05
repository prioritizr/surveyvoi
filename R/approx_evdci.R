#' Approximate expected value of the decision given current information
#'
#' Calculate the \emph{expected value of the conservation management decision
#' given current information}. This metric describes the value of the management
#' decision that is expected when the decision maker is limited to
#' existing biodiversity data (i.e. survey data and environmental niche models).
#'
#' @inheritParams approx_evdsi
#'
#' @param seed \code{integer} random seed used for generating states. Defaults
#'   to 500.
#'
#' @return \code{numeric} \code{vector} containing the estimate for each
#'   replicate.
#'
#' @seealso \code{\link{prior_probability_matrix}},
#' \code{\link{evdci}}.
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
#' # using approximate method with 100 replicates and 50 states per replicate
#' ev_prime_current <- approx_evdci(
#'     site_data, feature_data, c("f1", "f2"), c("p1", "p2"),
#'     "management_cost", "survey_sensitivity", "survey_specificity",
#'     "model_sensitivity", "model_specificity",
#'     "preweight", "postweight", "target",
#'     total_budget, n_approx_replicates = 100,
#'     n_approx_states_per_replicate = 50)
#'
#' # print approximate value
#' print(ev_prime_current)
#'
#' @export
approx_evdci <- function(
  site_data,
  feature_data,
  site_occupancy_columns,
  site_probability_columns,
  site_management_cost_column,
  feature_survey_sensitivity_column,
  feature_survey_specificity_column,
  feature_model_sensitivity_column,
  feature_model_specificity_column,
  feature_preweight_column,
  feature_postweight_column,
  feature_target_column,
  total_budget,
  site_management_locked_in_column = NULL,
  prior_matrix = NULL,
  n_approx_obj_fun_points = 1000,
  optimality_gap = 0,
  n_approx_replicates = 100,
  n_approx_states_per_replicate = 1000,
  method_approx_states = "weighted_without_replacement",
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
    ## prior_matrix
    inherits(prior_matrix, c("matrix", "NULL")),
    ## n_approx_obj_fun_points
    assertthat::is.number(n_approx_obj_fun_points),
    assertthat::noNA(n_approx_obj_fun_points),
    isTRUE(n_approx_obj_fun_points > 0),
    ## n_approx_replicates
    assertthat::is.count(n_approx_replicates),
    assertthat::noNA(n_approx_replicates),
    ## n_approx_states_per_replicate
    assertthat::is.count(n_approx_states_per_replicate),
    assertthat::noNA(n_approx_states_per_replicate),
    isTRUE(n_approx_states_per_replicate <=
           n_states(nrow(site_data), nrow(feature_data))),
    ## method_approx_states
    assertthat::is.string(method_approx_states),
    assertthat::noNA(method_approx_states),
    isTRUE(method_approx_states %in%
      c("uniform_with_replacement", "uniform_without_replacement",
        "weighted_with_replacement", "weighted_without_replacement")),
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

  # drop spatial information
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)

  # calculate prior matrix
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

  # main calculation
  withr::with_seed(seed, {
    out <- rcpp_approx_expected_value_of_decision_given_current_info_n_states(
      pij = pij,
      pu_costs = site_data[[site_management_cost_column]],
      pu_locked_in = site_management_locked_in,
      preweight = feature_data[[feature_preweight_column]],
      postweight = feature_data[[feature_postweight_column]],
      target = feature_data[[feature_target_column]],
      n_approx_obj_fun_points = n_approx_obj_fun_points,
      budget = total_budget,
      gap = optimality_gap,
      n_approx_replicates = n_approx_replicates,
      n_approx_states_per_replicate = n_approx_states_per_replicate,
      method_approx_states = method_approx_states)
  })

  # return result
  out
}
