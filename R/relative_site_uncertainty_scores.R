#' @include internal.R
NULL

#' Relative site uncertainty scores
#'
#' Calculate scores to describe the overall uncertainty of modelled species'
#' occupancy predictions for each site. Sites with greater scores are associated
#' with greater uncertainty. Note that these scores are relative to each other
#' and uncertainty values calculated using different matrices cannot be
#' compared to each other.
#'
#' @inheritParams evdci
#'
#' @details
#' The relative site uncertainty scores are calculated as joint Shannon's
#' entropy statistics. Since we assume that species occur independently of each
#' other, we can calculate these statistics separately for each species in each
#' site and then sum together the statistics for species in the same site:
#'
#' \enumerate{
#' \item Let \eqn{J} denote the set of sites (indexed by \eqn{j}),
#'   \eqn{I} denote the set of features (indexed by \eqn{i}), and
#'   \eqn{x_{ij}} denote the modelled probability of feature \eqn{i \in I}
#'   occurring in planning units \eqn{j \in J}.
#'
#' \item Next, we will calculate the Shannon's entropy statistic for each
#'   species in each site:
#'   \eqn{y_{ij} = - \big( (x_ij \mathit{log}_2 x_{ij}) + (1 - x_ij \mathit{log}_2 1 - x_{ij}) \big) }
#'
#' \item Finally, we will sum the entropy statistics together for each site:
#'   \eqn{s_{j} = \sum_{i \in I} y_{ij}}
#'
#' }
#'
#' @return `numeric` `vector` of uncertainty scores. Note that
#'  these values are automatically rescaled between 0.01 and 1.
#'
#' @examples
#' # set seed for reproducibility
#' set.seed(123)
#'
#' # simulate data for 3 features and 5 sites
#' x <- tibble::tibble(x = rnorm(5), y = rnorm(5),
#'                     p1 = c(0.5, 0, 1, 0, 1),
#'                     p2 = c(0.5, 0.5, 1, 0, 1),
#'                     p3 = c(0.5, 0.5, 0.5, 0, 1))
#' x <- sf::st_as_sf(x, coords = c("x", "y"))
#'
#' # print data,
#' # we can see that site (row) 3 has the least certain predictions
#' # because it has many values close to 0.5
#' print(x)
#'
#' # plot sites' occupancy probabilities
#' plot(x[, c("p1", "p2", "p3")], pch = 16, cex = 3)
#'
#' # calculate scores
#' s <- relative_site_uncertainty_scores(x, c("p1", "p2", "p3"))
#'
#' # print scores,
#' # we can see that site 3 has the highest uncertainty score
#' print(s)
#'
#' # plot sites' uncertainty scores
#' x$s <- s
#' plot(x[, c("s")], pch = 16, cex = 3)
#'
#' @export
relative_site_uncertainty_scores <- function(
  site_data, site_probability_columns) {
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
  # calculate uncertainty values
  ## convert to matrix format
  x <- sf::st_drop_geometry(site_data)
  x <- t(as.matrix(x[, site_probability_columns, drop = FALSE]))
  ## calculate entropy scores for each site in each species
  x[] <- -1 * ((x * log2(x)) + ((1 - x) * log2(1 - x)))
  ## convert non-finite values to zeros, because log2(0) = NaN
  x[!is.finite(x)] <- 0
  ## calculate entropy score for each site
  colSums(x)
}
