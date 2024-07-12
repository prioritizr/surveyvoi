#' @include internal.R
NULL

#' Greedy heuristic prioritization
#'
#' Generate a prioritization for protected area establishment.
#'
#' @inheritParams evdci
#'
#' @details
#' The prioritization is generated using a greedy heuristic algorithm.
#' The objective function for this algorithm is calculated by:
#' (i) estimating the probability that each species meets its target,
#' and (ii) calculating the sum of these probabilities.
#' Note that this function underpins the value of information calculations
#' because it is used to assess a potential management decision
#' given updated information on the presence of particular species
#' in particular sites.
#'
#' @return
#' A `list` containing the following elements:
#' \describe{
#' \item{x}{`logical` vector indicating if each site is selected for
#' protection or not.}
#' \item{objval}{`numeric` value denoting the objective value for
#' the prioritization.}
#' }
#'
#' @examples
#' # set seeds for reproducibility
#' set.seed(123)
#'
#' # load example site data
#' data(sim_sites)
#' print(sim_sites)
#'
#' # load example feature data
#' data(sim_features)
#' print(sim_features)
#'
#' # set total budget for managing sites for conservation
#'  # (i.e. 50% of the cost of managing all sites)
#' total_budget <- sum(sim_sites$management_cost) * 0.5
#'
#' # generate reserve selection prioritization
#' results <- greedy_heuristic_prioritization(
#'   sim_sites, sim_features,
#'   c("p1", "p2", "p3"),
#'   "management_cost",
#'   "target",
#'   total_budget
#' )
#'
#' # print results
#' print(results)
#'
#' @export
greedy_heuristic_prioritization <- function(
  site_data,
  feature_data,
  site_probability_columns,
  site_management_cost_column,
  feature_target_column,
  total_budget,
  site_management_locked_in_column = NULL,
  site_management_locked_out_column = NULL
) {
  # assert arguments are valid
  assertthat::assert_that(
    ## site_data
    inherits(site_data, "sf"), ncol(site_data) > 0,
    nrow(site_data) > 0,
    ## feature_data
    inherits(feature_data, "data.frame"), ncol(feature_data) > 0,
    nrow(feature_data) > 0,
    ## feature_target_column
    assertthat::is.string(feature_target_column),
    all(assertthat::has_name(feature_data, feature_target_column)),
    is.numeric(feature_data[[feature_target_column]]),
    assertthat::noNA(feature_data[[feature_target_column]]),
    all(feature_data[[feature_target_column]] >= 0),
    ## site_management_cost_column
    assertthat::is.string(site_management_cost_column),
    all(assertthat::has_name(site_data, site_management_cost_column)),
    is.numeric(site_data[[site_management_cost_column]]),
    assertthat::noNA(site_data[[site_management_cost_column]]),
    ## total_budget
    assertthat::is.number(total_budget), assertthat::noNA(total_budget),
    isTRUE(total_budget > 0)
  )
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
  ## validate pij values
  validate_site_probability_data(site_data, site_probability_columns)

  # drop spatial information
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)

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

  ## validate it is possible to select enough planning units to meet
  ## the target for even a single feature,
  ## if not, then it's not possible to generate a single solution
  ## with an objective value > 0
  cheapest_sol <- rep(FALSE, nrow(site_data))
  cheapest_sol[site_management_locked_in] <- TRUE
  n_extra_needed <-
  if (min(feature_data[[feature_target_column]]) > sum(cheapest_sol)) {
    idx <- which(!site_management_locked_out & !site_management_locked_in)
    idx_costs <- site_data[[site_management_cost_column]]
    idx_costs <- idx_costs[idx]
    idx_order <- order(idx_costs)
    idx_add <- cumsum(idx_costs[idx_order]) <= total_budget
    sel_idx <- idx[idx_order[which(idx_add)]]
    cheapest_sol[sel_idx] <- TRUE
  }
  if (
    sum(cheapest_sol) <= min(feature_data[[feature_target_column]])
  ) {
    ### throw warning
    warning(
      paste(
        "it is not possible to select enough planning units to meet",
        "any targets for even a single feature",
        "(given the budget and locked out planning units)"
      ),
      call. = TRUE, immediate. = FALSE
    )
    ## return solution
    return(list(x = cheapest_sol, objval = 0))
  }

  # create prior data
  prior_data <-
    t(as.matrix(site_data[, site_probability_columns, drop = FALSE]))

  # main calculation
  out <- rcpp_greedy_heuristic_prioritization(
    rij = prior_data,
    pu_costs = site_data[[site_management_cost_column]],
    pu_locked_in = as.numeric(site_management_locked_in),
    pu_locked_out = as.numeric(site_management_locked_out),
    target = round(feature_data[[feature_target_column]]),
    budget = total_budget
  )

  # update objective value with exact calculation
  out$objval <- rcpp_expected_value_of_action(
    out$x,
    prior_data,
    round(feature_data[[feature_target_column]])
  )

  # return result
  out

}
