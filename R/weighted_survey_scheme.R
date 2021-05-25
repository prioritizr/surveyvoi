#' @include internal.R ilp.R env_div_survey_scheme.R
NULL

#' Weighted survey scheme
#'
#' Generate a survey scheme by selecting the set of sites with the greatest
#' overall weight value, a maximum budget for the survey scheme.
#'
#' @param weight_column `character` name of the column in the argument
#'   to `site_data` with the weights for each site.
#'
#' @inheritParams env_div_survey_scheme
#'
#' @inherit env_div_survey_scheme return
#'
#' @details
#'   Let \eqn{J} denote the set of sites (indexed by \eqn{j}), and let
#'   \eqn{b} denote the maximum budget available for surveying the sites.
#'   Next, let \eqn{c_j} represent the cost of surveying each site
#'   \eqn{j \in J}, and \eqn{w_j} denote the relative value (weight) for
#'   surveying each site \eqn{j \in J}. The set of sites with the greatest
#'   overall weight values, subject to a given budget can the be identified by
#'   solving the following integer programming problem. Here,
#'   \eqn{x_j} is the binary decision variable indicating each if site
#'   is selected in the survey scheme or not.
#'
#'   \deqn{\mathit{Maximize} \space \sum_{j \in J} x_j w_i \\
#'   \mathit{subject \space to} \\
#'   \sum_{j \in J} x_j c_j \leq b}{
#'   Minimize sum_j^J (xi * wi) subject to sum_j^J (xi * ci) <= b}
#'
#' @inheritSection env_div_survey_scheme Solver
#'
#' @examples
#' # set seed for reproducibility
#' set.seed(123)
#'
#' # simulate data
#' x <- sf::st_as_sf(
#'   tibble::tibble(x = rnorm(4), y = rnorm(4),
#'                  w = c(0.01, 10, 8, 1),
#'                  cost = c(1, 1, 1, 1)),
#'   coords = c("x", "y"))
#'
#' # plot site' locations and color by weight values
#' plot(x[, "w"], pch = 16, cex = 3)
#'
#' # generate scheme without any sites locked in
#' s <- weighted_survey_scheme(x, cost_column = "cost", survey_budget = 2,
#'                              weight_column = "w")
#'
#' # print solution
#' print(s)
#'
#' # plot solution
#' x$s <- c(s)
#' plot(x[, "s"], pch = 16, cex = 3)
#' @export
weighted_survey_scheme <- function(
  site_data, cost_column, survey_budget, weight_column, locked_in_column = NULL,
  locked_out_column = NULL, solver = "auto", verbose = FALSE) {
  # assert that arguments are valid
  assertthat::assert_that(
    ## site_data
    inherits(site_data, "sf"), nrow(site_data) > 0, ncol(site_data) > 0,
    ## cost_column
    assertthat::is.string(cost_column), assertthat::noNA(cost_column),
    all(assertthat::has_name(site_data, cost_column)),
    is.numeric(site_data[[cost_column]]),
    assertthat::noNA(site_data[[cost_column]]),
    all(site_data[[cost_column]] >= 0),
    ## survey_budget
    is.numeric(survey_budget), assertthat::noNA(survey_budget),
    all(survey_budget >= 0),
    ## weight_column
    assertthat::is.string(weight_column),
    all(assertthat::has_name(site_data, weight_column)),
    is.numeric(site_data[[weight_column]]),
    assertthat::noNA(site_data[[weight_column]]),
    all(site_data[[weight_column]] > 0))
  if (!is.null(locked_in_column)) {
    ## locked_in_column
    assertthat::assert_that(
      assertthat::is.string(locked_in_column),
      all(assertthat::has_name(site_data, locked_in_column)),
      is.logical(site_data[[locked_in_column]]),
      assertthat::noNA(site_data[[locked_in_column]]))
  }
  if (!is.null(locked_out_column)) {
    ## locked_out_column
    assertthat::assert_that(
      assertthat::is.string(locked_out_column),
      all(assertthat::has_name(site_data, locked_out_column)),
      is.logical(site_data[[locked_out_column]]),
      assertthat::noNA(site_data[[locked_out_column]]))
  }

  # set locked in sites
  if (is.null(locked_in_column)) {
    locked_in <- rep(FALSE, nrow(site_data))
  } else {
    locked_in <- site_data[[locked_in_column]]
  }

  # set locked out sites
  if (is.null(locked_out_column)) {
    locked_out <- rep(FALSE, nrow(site_data))
  } else {
    locked_out <- site_data[[locked_out_column]]
  }

  # return survey schemes
  weight_based_prioritizations(
    site_data[[weight_column]], survey_budget, site_data[[cost_column]],
    locked_in, locked_out,
    solver, verbose)
}
