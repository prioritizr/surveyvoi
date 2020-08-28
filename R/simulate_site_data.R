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
#' @param max_number_surveys_per_site \code{integer} maximum number of
#'   surveys per site in the simulated data. Defaults to 5.
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
#' @return \code{\link[sf]{sf}} object with site data. Columns that start with
#'  (i) \code{"f"} (e.g. \code{"d1"}) contain the proportion of
#'   times that each feature was detected in each site,
#'  (ii) \code{"n"} (e.g. \code{"n1"}) contain the number of
#'   of surveys for each feature within each site,
#'  (iii) \code{"p"} (e.g. \code{"p1"}) contain prior
#'  probability data,  and
#'  (iv) \code{"e"} (e.g. \code{"e1"}) contain environmental
#'  data. Note that columns that contain the same integer value (excepting
#'  environmental data columns) correspond to the same feature
#'  (e.g. \code{"d1"}, \code{"n1"}, \code{"p1"} contain data that correspond
#'  to the same feature).
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
#' # plot feature detection data
#' plot(d[, c("f1", "f2", "f3", "f4")], axes = TRUE, pch = 16, cex = 2)
#'
#' # plot feature survey effort
#' plot(d[, c("n1", "n2", "n3", "n4")], axes = TRUE, pch = 16, cex = 2)
#'
#' # plot feature prior probability data
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
                               max_number_surveys_per_site = 5,
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
    ## max_number_surveys_per_site
    assertthat::is.count(max_number_surveys_per_site),
    assertthat::noNA(max_number_surveys_per_site),
    ## output_probabilities
    assertthat::is.flag(output_probabilities),
    assertthat::noNA(output_probabilities))
  # site locations
  site_data <- tibble::tibble(x = stats::runif(n_sites),
                              y = stats::runif(n_sites))
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
  feature_coef <- matrix(stats::rnorm(n_features * n_env_vars, sd = 8),
                         ncol = n_env_vars, nrow = n_features)
  # probability data
  pij <- matrix(NA, ncol = n_features, nrow = n_sites)
  for (i in seq_len(n_features))
    pij[, i] <- stats::plogis(c(env_data  %*% feature_coef[i, ]))
  pij[] <- round(pij[], 3)
  # determine number of surveys per site
  n_surveys <- sample.int(max_number_surveys_per_site, n_sites, replace = TRUE)
  # simulate survey detections within each site
  p_det <- pij
  for (i in seq_len(n_features)) {
    for (j in seq_len(n_sites)) {
      p_det[j, i] <- mean(stats::rbinom(n_surveys[j], 1, pij[j, i]))
    }
  }
  # randomly specify sites missing data
  missing_sites <- sample.int(n_sites,
    round(proportion_of_sites_missing_data * n_sites), replace = FALSE)
  # set missing site data to zeros
  p_det[missing_sites, ] <- 0
  # compile data
  ## detection data
  for (i in seq_len(n_features))
    site_data[[paste0("f", i)]] <- p_det[, i]
  ## number of surveys data
  for (i in seq_len(n_features))
    site_data[[paste0("n", i)]] <- replace(n_surveys, missing_sites, 0)
  ## env data
  for (i in seq_len(n_env_vars))
    site_data[[paste0("e", i)]] <- env_data[, i]
  ## pij data
  if (output_probabilities) {
    for (i in seq_len(ncol(pij)))
      site_data[[paste0("p", i)]] <- pij[, i]
  }
  # return result
  sf::st_as_sf(site_data, coords = c("x", "y"))
}
