#' @include internal.R
NULL

#' Relative site richness scores
#'
#' Calculate relative site richness scores. Sites with greater scores
#' are predicted to be more likely to contain more species. Note that these
#' scores are relative to each other and scores calculated using different
#' matrices cannot be compared to each other.
#'
#' @inheritParams relative_site_uncertainty_scores
#'
#' @details
#' The relative site richness scores are calculated using the following
#' procedure:
#'
#' \enumerate{
#' \item Let \eqn{J} denote the set of sites (indexed by \eqn{j}),
#'   \eqn{I} denote the set of features (indexed by \eqn{i}), and
#'   \eqn{x_{ij}} denote the modelled probability of feature \eqn{i \in I}
#'   occurring in planning units \eqn{j \in J}.
#'
#' \item Next, we will sum the values for each site:
#'   \eqn{y_j = \sum_{i \in I} x_{ij}}.
#'
#' \item Finally, we will linearly rescale the \eqn{y_j} values between 0.01
#'   and 1 to produce the scores.
#'
#' }
#'
#' @return \code{numeric} \code{vector} of richness scores. Note that
#'  these values are automatically rescaled between 0.01 and 1.
#'
#' @examples
#' # set seed for reproducibility
#' set.seed(123)
#'
#' # simulate data for 3 features and 4 planning units
#' x <- tibble::tibble(x = rnorm(4), y = rnorm(4),
#'                      p1 = c(0.095, 0.032, 0.5, 0.924),
#'                      p2 = c(0.023, 0.014, 0.4, 0.919),
#'                      p3 = c(0.075, 0.046, 0.9, 0.977))
#' x <- sf::st_as_sf(x, coords = c("x", "y"))
#'
#' # print data,
#' # we can see that the fourth site has the highest modelled probabilities of
#' # occupancy across all species
#' print(x)
#'
#' # plot sites' occupancy probabilities
#' plot(x[, c("p1", "p2", "p3")], pch = 16, cex = 3)
#'
#' # calculate scores
#' s <- relative_site_richness_scores(x, c("p1", "p2", "p3"))
#'
#' # print scores,
#' # we can see that site 4 has the highest richness score
#' print(s)
#'
#' # plot sites' richness scores
#' x$s <- s
#' plot(x[, c("s")], pch = 16, cex = 3)
#'
#' @export
relative_site_richness_scores <- function(site_data, site_probability_columns) {
  # assert that arguments are valid
  assertthat::assert_that(
    inherits(site_data, "sf"),
    is.character(site_probability_columns),
    assertthat::noNA(site_probability_columns),
    all(assertthat::has_name(site_data, site_probability_columns)))
  ## validate pij values
  assertthat::assert_that(
    all(sapply(site_probability_columns,
               function(x) is.numeric(site_data[[x]]))),
    msg = "site_data values in site_probability_columns must be numeric")
  assertthat::assert_that(
    all(sapply(site_probability_columns,
               function(x) all(!is.na(site_data[[x]])))),
    msg = "site_data values in site_probability_columns must not be NA")
  assertthat::assert_that(
    all(sapply(site_probability_columns,
               function(x)
                 all(site_data[[x]] >= 0) && all(site_data[[x]] <= 1))),
    msg = paste("site_data values in site_probability_columns must be between",
                "0 and 1"))
  # calculate scores
  ## convert to matrix format
  x <- sf::st_drop_geometry(site_data)
  out <- t(as.matrix(x[, site_probability_columns, drop = FALSE]))
  ## set non-finite values to zero, e.g. because a species has exactly
  ## the same predicted probabilities of occupancy for all sites
  out[!is.finite(out)] <- 0
  ## sum scores for each site
  out <- colSums(out)
  ## rescale numbers to between 0.01 and 1
  out <- scales::rescale(out, to = c(0.01, 1))
  # return values
  out
}
