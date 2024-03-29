#' @include internal.R
NULL

#' Simulate site data
#'
#' Simulate site data for developing simulated survey schemes.
#'
#' @param n_sites `integer` number of sites.
#'
#' @param n_features `integer` number of features.
#'
#' @param proportion_of_sites_missing_data `numeric` proportion of sites
#'   that do not have existing presence/absence data. Values must be between
#'   zero and one.
#'
#' @param n_env_vars `integer` number of environmental variables for
#'   simulating feature distributions. Defaults to 3.
#'
#' @param survey_cost_intensity `numeric` intensity of the costs of
#'   surveying sites. Larger values correspond to larger costs on average.
#'   Defaults to 20.
#'
#' @param survey_cost_scale `numeric` value corresponding to
#'   the spatial homogeneity of the survey costs. Defaults to 5.
#'
#' @param management_cost_intensity `numeric` intensity of the costs of
#'   average cost of managing sites for conservation. Defaults to 100.
#'
#' @param management_cost_scale `numeric` value corresponding to
#'   the spatial homogeneity of the survey costs. Defaults to 30.
#'
#' @param max_number_surveys_per_site `integer` maximum number of
#'   surveys per site in the simulated data. Defaults to 5.
#'
#' @param output_probabilities `logical` value indicating if
#'   probability values of occupancy should be output or not. Defaults
#'   to `TRUE`.
#'
#' @return A [sf::sf()] object with site data.
#'  The `"management_cost"` column contains the site protection costs,
#'  and the `"survey_cost"` column contains the costs for surveying
#'  each site.
#'  Additionally, columns that start with
#'  (i) `"f"` (e.g. `"f1"`) contain the proportion of
#'   times that each feature was detected in each site,
#'  (ii) `"n"` (e.g. `"n1"`) contain the number of
#'   of surveys for each feature within each site,
#'  (iii) `"p"` (e.g. `"p1"`) contain prior
#'  probability data,  and
#'  (iv) `"e"` (e.g. `"e1"`) contain environmental
#'  data. Note that columns that contain the same integer value (excepting
#'  environmental data columns) correspond to the same feature
#'  (e.g. `"d1"`, `"n1"`, `"p1"` contain data that correspond
#'  to the same feature).
#'
#' @examples
#' # set seed for reproducibility
#' set.seed(123)
#'
#' # simulate data
#' d <- simulate_site_data(n_sites = 10, n_features = 4, prop = 0.5)
#'
#' # print data
#' print(d, width = Inf)
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
#' @seealso [simulate_feature_data()]
#'
#' @export
simulate_site_data <- function(n_sites, n_features,
                               proportion_of_sites_missing_data,
                               n_env_vars = 3,
                               survey_cost_intensity = 20,
                               survey_cost_scale = 5,
                               management_cost_intensity = 100,
                               management_cost_scale = 30,
                               max_number_surveys_per_site = 5,
                               output_probabilities = TRUE) {
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
    ## survey_cost_scale
    assertthat::is.number(survey_cost_scale),
    assertthat::noNA(survey_cost_scale),
    all(survey_cost_scale >= 0),
    ## management_cost_intensity
    assertthat::is.number(management_cost_intensity),
    assertthat::noNA(management_cost_intensity),
    all(management_cost_intensity >= 0),
    ## management_cost_scale
    assertthat::is.number(management_cost_scale),
    assertthat::noNA(management_cost_scale),
    all(management_cost_scale >= 0),
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
  site_data$survey_cost <- round(c(
    simulate_gaussian_random_field(
      n = 1,
      coords = matrix(c(site_data$x, site_data$y), ncol = 2),
      mu = survey_cost_intensity,
      sd = survey_cost_scale,
      scale = survey_cost_scale,
      truncate = 1
    )
  ))
  site_data$management_cost <- round(c(
    simulate_gaussian_random_field(
      n = 1,
      coords = matrix(c(site_data$x, site_data$y), ncol = 2),
      mu = management_cost_intensity,
      sd = management_cost_scale,
      scale = management_cost_scale,
      truncate = 1
    )
  ))
  # environmental data
  env_data <- simulate_gaussian_random_field(
    n = n_env_vars,
    coords = matrix(c(site_data$x, site_data$y), ncol = 2),
    mu = 0,
    sd = 1,
    scale = 1.0
  )
  env_data <- apply(env_data, 2, function(x) scale(x))
  # feature model coefficients
  feature_coef <- matrix(stats::rnorm(n_features * n_env_vars, sd = 8),
                         ncol = n_env_vars, nrow = n_features)
  # probability data
  pij <- matrix(NA, ncol = n_features, nrow = n_sites)
  for (i in seq_len(n_features))
    pij[, i] <- stats::plogis(c(env_data %*% feature_coef[i, ]))
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

#' Simulate a Gaussian random field
#'
#' @param n `integer` Number of simulations to generate.
#'
#' @param coords `matrix` Matrix containing coordinates for simulating data.
#'
#' @param mu `numeric` Parameter for simulations.
#'
#' @param sd `numeric` Parameter for simulations.
#'
#' @param scale `numeric` Parameter for simulations.
#'
#' @param truncate `numeric` Threshold lowest value in simulated values.
#' Defaults to negative infinity (`-Inf`) such that no data are truncated.
#"
#' @inherit simulate_random_field return
#'
#' @noRd
simulate_gaussian_random_field <- function(n, coords, mu, sd, scale,
                                           truncate = -Inf) {
  # assert valid arguments
  assertthat::assert_that(
    assertthat::is.count(n),
    assertthat::noNA(n),
    is.matrix(coords),
    assertthat::noNA(c(coords)),
    nrow(coords) >= 1,
    ncol(coords) == 2,
    assertthat::is.number(mu),
    assertthat::noNA(mu),
    assertthat::is.number(scale),
    assertthat::noNA(scale),
    assertthat::is.number(sd),
    assertthat::noNA(sd),
    assertthat::is.number(truncate),
    assertthat::noNA(truncate)
  )
  # main processing
  mu <- rep(mu, nrow(coords))
  p <- nrow(coords)
  chol_d <- chol(exp(-scale * as.matrix(stats::dist(coords))))
  out <- t(
    matrix(stats::rnorm(n = n * p, sd = sd), ncol = p) %*%
    chol_d + rep(mu, rep(n, p))
  )
  # ensure matrix output
  if (!is.matrix(out)) {
    out <- matrix(out, ncol = 1)
  }
  # truncate outputs
  if (any(out < truncate)) {
    out[which(out < truncate)] <- truncate
  }
  # return result
  out
}
