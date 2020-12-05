#' @include internal.R ilp.R
NULL

#' Environmental diversity survey scheme
#'
#' Generate a survey scheme by maximizing the diversity of environmental
#' conditions that are surveyed.
#' Please note that this function requires the Gurobi optimization software
#' (<https://www.gurobi.com/>) and the \pkg{gurobi} R package
#' (installation instructions available for [Linux](https://www.gurobi.com/documentation/9.1/quickstart_linux/r_ins_the_r_package.html), [Windows](https://www.gurobi.com/documentation/9.1/quickstart_windows/r_ins_the_r_package.html), and [Mac OS](https://www.gurobi.com/documentation/9.1/quickstart_mac/r_ins_the_r_package.html)).
#'
#' @param site_data [sf::sf()] object containing the candidate survey
#'   sites.
#'
#' @param cost_column `character` name of the column in the argument to
#'   the argument to `site_data` that contains the cost for surveying each
#'   site. No missing (`NA`) values are permitted.
#'
#' @param survey_budget `numeric` vector of maximum budgets for the survey
#'   schemes. No missing (`NA`) values are permitted.
#'
#' @param env_vars_columns `character` vector names of the columns in
#'   the argument to `site_data` that contain `numeric` environmental
#'   variables. No missing (`NA`) values are permitted.
#'
#' @param method `character` name of the distance metric to use for
#'   calculating environmental dissimilarity scores. See
#'   [vegan::vegdist()] documentation the `method` parameter
#'   for other available distance metrics and more information.
#'   No missing (`NA`) values are permitted.
#'   Defaults to `"mahalanobis"` for Mahalanobis distances.
#'
#' @param locked_in_column `character` (optional) name of the column in
#'   the argument to `site_data` that contains `logical`
#'   (`TRUE`/ `FALSE`) values indicating if certain sites should be
#'   locked into the survey scheme.
#'   No missing (`NA`) values are permitted.
#'   Defaults to `NULL` such that no sites are locked in.
#'
#' @param locked_out_column `character` (optional) name of the column in
#'   the argument to `site_data` that contains `logical`
#'   (`TRUE`/ `FALSE`) values indicating if certain sites should be
#'   locked out of the survey scheme.
#'   No missing (`NA`) values are permitted.
#'   Defaults to `NULL` such that no sites are locked out.
#'
#' @param exclude_locked_out `logical` should locked out planning units
#'  be entirely excluded from the optimization process?
#'  Defaults to `FALSE`.
#'
#' @param verbose `logical` indicating if information should be
#'   printed while generating survey scheme(s). Defaults to `FALSE`.
#'
#' @details The integer programming formulation of the environmental diversity
#'   reserve selection problem (Faith & Walker 1996) is used to generate survey
#'   schemes.
#'
#' @references
#' Faith DP & Walker PA (1996) Environmental diversity: on the best-possible
#' use of surrogate data for assessing the relative biodiversity of sets of
#' areas. *Biodiversity & Conservation*, **5**, 399--415.
#'
#' @return `matrix` of `logical` (`TRUE`/ `FALSE`)
#'   values indicating if a site is selected in a scheme or not. Columns
#'   correspond to sites, and rows correspond to different schemes.
#'
#' @examples
#' \dontrun{
#' # set seed for reproducibility
#' set.seed(123)
#'
#' # simulate data
#'  x <- sf::st_as_sf(
#'    tibble::tibble(x = rnorm(4), y = rnorm(4),
#'                   v1 = c(0.1, 0.2, 0.3, 10), # environmental axis 1
#'                   v2 = c(0.1, 0.2, 0.3, 10), # environmental axis 2
#'                   cost = rep(1, 4)),
#'     coords = c("x", "y"))
#'
#' # plot the sites' environmental conditions
#' plot(x[, c("v1", "v2")], pch = 16, cex = 3)
#'
#' # generate scheme with a budget of 2
#' s <- env_div_survey_scheme(x, "cost", 2, c("v1", "v2"), "mahalanobis")
#'
#' # print scheme
#' print(s)
#'
#' # plot scheme
#' x$scheme <- c(s)
#' plot(x[, "scheme"], pch = 16, cex = 3)
#' }
#' @export
env_div_survey_scheme <- function(
  site_data, cost_column, survey_budget, env_vars_columns,
  method = "mahalanobis", locked_in_column = NULL, locked_out_column = NULL,
  exclude_locked_out = FALSE, verbose = FALSE) {
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
    ## env_vars_columns
    is.character(env_vars_columns), assertthat::noNA(env_vars_columns),
    length(env_vars_columns) > 0,
    all(assertthat::has_name(site_data, env_vars_columns)),
    all(sapply(cost_column, function(z) is.numeric(site_data[[z]]))),
    all(sapply(cost_column, function(z) assertthat::noNA(site_data[[z]]))),
    ## exclude_locked_out
    assertthat::is.flag(exclude_locked_out),
    assertthat::noNA(exclude_locked_out),
    ## method
    assertthat::is.string(method),
    ## survey_budget
    is.numeric(survey_budget), assertthat::noNA(survey_budget),
    all(survey_budget >= 0))
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

  # exclude locked out planning units if specified
  if (isTRUE(exclude_locked_out)) {
    cand_site_data <- site_data[!locked_out, , drop = FALSE]
    cand_locked_in <- locked_in[!locked_out]
    cand_locked_out <- locked_out[!locked_out]

  } else {
    cand_site_data <- site_data
    cand_locked_in <- locked_in
    cand_locked_out <- locked_out
  }

  # create environmental distance matrix
  env_dists <-
    as.matrix(vegan::vegdist(
      sf::st_drop_geometry(cand_site_data)[, env_vars_columns, drop = FALSE],
      method = method))

  # run optimization
  result <- distance_based_prioritizations(env_dists, survey_budget,
    cand_site_data[[cost_column]], cand_locked_in, cand_locked_out, verbose)

  # prepare output
  if (isTRUE(exclude_locked_out)) {
    out <- matrix(FALSE, nrow = nrow(result), ncol = nrow(site_data))
    idx <- which(!locked_out)
    for (i in seq_along(idx))
      out[, idx[i]] <- result[, i]
  } else {
    out <- result
  }

  # return output
  out
}
