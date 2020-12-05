#' @include internal.R ilp.R env_div_survey_scheme.R
NULL

#' Geographic coverage survey scheme
#'
#' Generate a survey scheme by maximizing the geographic coverage
#' of surveys.
#' Please note that this function requires the Gurobi optimization software
#' (<https://www.gurobi.com/>) and the \pkg{gurobi} R package
#' (installation instructions available for [Linux](https://www.gurobi.com/documentation/9.1/quickstart_linux/r_ins_the_r_package.html), [Windows](https://www.gurobi.com/documentation/9.1/quickstart_windows/r_ins_the_r_package.html), and [Mac OS](https://www.gurobi.com/documentation/9.1/quickstart_mac/r_ins_the_r_package.html)).
#'
#' @inheritParams env_div_survey_scheme
#'
#' @inherit env_div_survey_scheme return
#'
#' @details The integer programming formulation of the *p*-Median
#'   problem (Daskin & Maass 2015) is used to generate survey schemes.
#'
#' @references
#' Daskin MS & Maass KL (2015) The p-median problem. In *Location Science*
#' (pp. 21-45). Springer, Cham.
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
#' # plot the sites' locations
#' plot(st_geometry(x), pch = 16, cex = 3)
#'
#' # generate scheme with a budget of 2
#' s <- geo_cov_survey_scheme(x, "cost", 2)
#'
#' # print scheme
#' print(s)
#'
#' # plot scheme
#' x$scheme <- c(s)
#' plot(x[, "scheme"], pch = 16, cex = 3)
#' }
#' @export
geo_cov_survey_scheme <- function(
  site_data, cost_column, survey_budget, locked_in_column = NULL,
  locked_out_column = NULL, exclude_locked_out = FALSE, verbose = FALSE) {
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
    ## exclude_locked_out
    assertthat::is.flag(exclude_locked_out),
    assertthat::noNA(exclude_locked_out),
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

  # create geographic distance matrix
  if (all(sapply(sf::st_geometry(cand_site_data), inherits, "POINT"))) {
    geo_dists <- as.matrix(
      stats::dist(methods::as(cand_site_data, "Spatial")@coords,
                  method = "euclidean"))
  } else {
    geo_dists <- sf::st_distance(cand_site_data)
    geo_dists <- matrix(as.numeric(geo_dists[]),
                        ncol = ncol(geo_dists),
                        nrow = nrow(geo_dists))
  }

  # rescale distances to avoid issues with large values
  geo_dists[] <- scales::rescale(geo_dists[], to = c(0, 1e+4))

  # run optimization
  result <- distance_based_prioritizations(geo_dists, survey_budget,
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
