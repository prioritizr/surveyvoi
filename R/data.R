#' @include internal.R
NULL

#' Simulated datasets
#'
#' Simulated data for prioritizing sites for ecological surveys.
#'
#' @usage data(sim_features)
#' @usage data(sim_sites)
#'
#' @format
#' \describe{
#'   \item{sim_sites}{[`sf::sf()`] object.}
#'   \item{sim_features}{[`tibble::tibble()`] object.}
#' }.
#'
#' @details
#' The simulated datasets provide data for six sites and three features. The
#' sites can potentially acquired for protected area establishment. However,
#' existing information on the spatial distribution of the features is
#' incomplete. Only some of the sites have existing ecological survey data.
#' To help inform management decisions, species distribution models have been
#' fitted to predict the probability of each species occupying each site.
#'
#' \describe{
#'   \item{`sim_sites`}{This object describes the sites and contains the
#'     following data: cost of surveying the sites (`survey_cost` column),
#'     cost of acquiring sites for conservation (`management_cost` column),
#'     results from previous ecological surveys (`f1`, `f2`, `f3` columns),
#'     previous survey effort (`n1`, `n2`, `n3` columns),
#'     environmental conditions of the sites (`e1`, `e2` columns),
#'     and modelled probability of the features occupying the sites
#'     (`p1`, `p2`, `p3` columns).}
#'
#'   \item{`sim_features`}{This object describes the features and contains
#'     the following data:
#'     the name of each feature (`name` column),
#'     whether each feature should be considered in future surveys
#'      (`survey` column),
#'     the sensitivity and specificity of the survey methodology for each
#      feature (`survey_sensitivity`, `survey_specificity` columns),
#'     the sensitivity and specificity of the species distribution model
#'     for each feature (`model_sensitivity`, `model_specificity` columns),
#'     and the representation target thresholds for each feature
#'     (`target` column).}
#' }
#'
#' @aliases sim_sites sim_features
#'
#' @keywords datasets
#'
#' @docType data
#'
#' @seealso These datasets were simulated using [`simulate_feature_data()`]
#'   and [`simulate_site_data()`].
#'
#' @examples
#' # load data
#' data(sim_sites, sim_features)
#'
#' # print feature data
#' print(sim_features, width = Inf)
#
#' # print site data
#' print(sim_sites, width = Inf)
#'
#' @name sim_data
NULL
