#' Number of states
#'
#' Calculate the total number of presence/absence states for a given
#' number of sites and features.
#'
#' @param n_sites `numeric` number of sites.
#'
#' @param n_features `numeric` number of features.
#'
#' @return `numeric` value.
#'
#' @examples
#' # calculate number of states for 3 sites and 2 features
#' n_states(n_sites = 2, n_features = 3)
#'
#' @export
n_states <- function(n_sites, n_features) {
  # assert that arguments are valid
  assertthat::assert_that(
    assertthat::is.count(n_sites), assertthat::noNA(n_sites),
    assertthat::is.count(n_features), assertthat::noNA(n_features))
  # return result
  rcpp_n_states(n_sites * n_features) + 1
}
