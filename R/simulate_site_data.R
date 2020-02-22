#' @include internal.R
NULL

#' Simulate site data
#'
#' Simulate site data for developing simulated survey schemes.
#'
#' @param n_sites \code{integer} number of sites.
#'
#' @param n_features \code{integer} number of features.
#'
#' @param proportion_of_sites_missing_data \code{numeric} proportion of sites
#'   that do not have existing presence/absence data. Values must be between
#'   zero and one.
#'
#' @param n_env_vars \code{integer} number of environmental variables for
#'   simulating feature distributions. Defaults to 3.
#'
#' @param survey_cost_intensity \code{numeric} intensity of the costs of
#'   surveying sites. Larger values correspond to larger costs on average.
#'   Defaults to 2.
#'
#' @param survey_cost_radius \code{numeric} value corresponding to
#'   the spatial homogeneity of the survey costs. Defaults to 0.1.
#'
#' @param management_cost_intensity \code{numeric} intensity of the costs of
#'   average cost of managing sites for conservation. Defaults to 100.
#'
#' @param management_cost_radius \code{numeric} value corresponding to
#'   the spatial homogeneity of the survey costs. Defaults to 0.1.
#'
#' @param output_probabilities \code{logical} value indicating if
#'   probability values of occupancy should be output or not. Defaults
#'   to \code{TRUE}.
#'
#' @details
#' The data are simulated using random fields
#' (\code{\link[RandomFields]{RFsimulate}}) to provide spatially auto-correlated
#' simulations.
#'
#' @return \code{\link[sf]{sf}} object with site data. Columns starting with
#'  (i) \code{"f"} (e.g. \code{f1}) contain presence/absence data,
#'  (ii) \code{"p"} (e.g. \code{p1}) contain prior probability data,
#'  (iii) \code{"e"} (e.g. \code{e1}) contain environmental data. Note
#'  that presence/absence and probability columns with the same integer suffix
#'  correspond to the same feature (e.g. \code{f1} and \code{p1} are the
#'  same feature).
#'
#' @examples
#' # set seed for reproducibility
#' set.seed(123)
#' RandomFields::RFoptions(seed = 123)
#'
#' # simulate data
#' d <- simulate_site_data(n_sites = 10, n_features = 4, prop = 0.5)
#'
#' # print data
#' print(d)
#'
#' # plot cost data
#' plot(d[, c("survey_cost", "management_cost")], axes = TRUE, pch = 16,
#'      cex = 2)
#'
#' # plot environmental data
#' plot(d[, c("e1", "e2", "e3")], axes = TRUE, pch = 16, cex = 2)
#'
#' # plot feature presence/absence data
#' plot(d[, c("f1", "f2", "f3", "f4")], axes = TRUE, pch = 16, cex = 2)
#;
#' # plot feature probability data
#' plot(d[, c("p1", "p2", "p3", "p4")], axes = TRUE, pch = 16, cex = 2)
#'
#' @seealso \code{\link{simulate_feature_data}}
#'
#' @export
simulate_site_data <- function(n_sites, n_features,
                               proportion_of_sites_missing_data,
                               n_env_vars = 3,
                               survey_cost_intensity = 10,
                               survey_cost_radius = 0.5,
                               management_cost_intensity = 100,
                               management_cost_radius = 0.5,
                               output_probabilities = TRUE) {
  # assert arguments are valid
  assertthat::assert_that(
    ## n_sites
    assertthat::is.count(n_sites), assertthat::noNA(n_sites),
    ## n_features
    assertthat::is.count(n_features), assertthat::noNA(n_features),
    ## n_env_vars
    assertthat::is.count(n_env_vars), assertthat::noNA(n_env_vars),
    ## proportion_of_sites_missing_data
    assertthat::is.number(proportion_of_sites_missing_data),
    assertthat::noNA(proportion_of_sites_missing_data),
    all(proportion_of_sites_missing_data >= 0),
    all(proportion_of_sites_missing_data <= 1),
    ## survey_cost_intensity
    assertthat::is.number(survey_cost_intensity),
    assertthat::noNA(survey_cost_intensity),
    all(survey_cost_intensity >= 0),
    ## survey_cost_radius
    assertthat::is.number(survey_cost_radius),
    assertthat::noNA(survey_cost_radius),
    all(survey_cost_radius >= 0),
    ## management_cost_intensity
    assertthat::is.number(management_cost_intensity),
    assertthat::noNA(management_cost_intensity),
    all(management_cost_intensity >= 0),
    ## management_cost_radius
    assertthat::is.number(management_cost_radius),
    assertthat::noNA(management_cost_radius),
    all(management_cost_radius >= 0),
    ## output_probabilities
    assertthat::is.flag(output_probabilities),
    assertthat::noNA(output_probabilities))
  # site locations
  site_data <- tibble::tibble(x = runif(n_sites), y = runif(n_sites))
  # cost data
  site_data$survey_cost <- suppressMessages(
    RandomFields::RFsimulate(
      x = site_data$x, y = site_data$y,
      model = RandomFields::RPpoisson(
        intensity = survey_cost_intensity,
        phi = RandomFields::RMtruncsupport(
          phi = RandomFields::RMgauss(),
          radius = survey_cost_radius)))[[1]])
  site_data$management_cost <- suppressMessages(
    RandomFields::RFsimulate(
      x = site_data$x, y = site_data$y,
      model = RandomFields::RPpoisson(
        intensity = management_cost_intensity,
        phi = RandomFields::RMtruncsupport(
          phi = RandomFields::RMgauss(),
          radius = management_cost_radius)))[[1]])
  # environmental data
  env_data <- as.matrix(suppressMessages(
    RandomFields::RFsimulate(model = RandomFields::RMgauss(),
                              x = site_data$x, y = site_data$y,
                              n = n_env_vars)@data))
  env_data <- apply(env_data, 2, function(x) scale(x))
  # feature model coefficients
  feature_coef <- matrix(rnorm(n_features * n_env_vars, sd = 5),
                         ncol = n_env_vars, nrow = n_features)
  # probability data
  pij <- matrix(NA, ncol = n_features, nrow = n_sites)
  for (i in seq_len(n_features))
    pij[, i] <- plogis(c(env_data  %*% feature_coef[i, ]))
  # presence/absence data
  rij <- pij
  rij[] <- rbinom(n_sites * n_features, 1, c(pij))
  # randomly specify sites missing data
  missing_sites <- sample.int(n_sites,
    round(proportion_of_sites_missing_data * n_sites), replace = FALSE)
  # set missing site data to NA
  rij[missing_sites, ] <- NA_real_
  # compile data
  ## env data
  for (i in seq_len(ncol(env_data)))
    site_data[[paste0("e", i)]] <- env_data[, i]
  ## rij data
  for (i in seq_len(ncol(rij)))
    site_data[[paste0("f", i)]] <- rij[, i]
  ## pij data
  if (output_probabilities) {
    for (i in seq_len(ncol(pij)))
      site_data[[paste0("p", i)]] <- pij[, i]
  }
  # return result
  st_as_sf(site_data, coords = c("x", "y"))
}
