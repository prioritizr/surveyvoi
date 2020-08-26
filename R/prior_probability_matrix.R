#' Prior probability matrix
#'
#' Create prior probability matrix for the value of information analysis.
#'
#' @inheritParams evdci
#'
#' @return \code{matrix} object containing the prior probabilities of each
#'   feature occupying each site. Each row corresponds to a different
#'   feature and each column corresponds to a different site.
#'
#' @examples
#' # set seeds for reproducibility
#' library(RandomFields)
#' set.seed(123)
#' RFoptions(seed = 123)
#'
#' # simulate data
#' site_data <- simulate_site_data(n_sites = 5, n_features = 2,
#'                                 prop = 0.5)
#' feature_data <- simulate_feature_data(n_features = 2, prop = 1)
#'
#' # preview simulated data
#' print(site_data)
#' print(feature_data)
#'
#' # calculate prior probability matrix
#' prior_matrix <- prior_probability_matrix(
#'   site_data, feature_data, c("f1", "f2"), c("p1", "p2"),
#'   "survey_sensitivity", "survey_specificity", "model_sensitivity",
#'   "model_specificity")
#'
#' # preview prior probability matrix
#' print(prior_matrix)
#' @export
prior_probability_matrix <- function(
  site_data, feature_data,
  site_detection_columns, site_n_surveys_columns,
  site_probability_columns,
  feature_survey_sensitivity_column, feature_survey_specificity_column,
  feature_model_sensitivity_column, feature_model_specificity_column) {
  # assert arguments are valid
  assertthat::assert_that(
    ## site_data
    inherits(site_data, c("sf", "data.frame")),
    ncol(site_data) > 0, nrow(site_data) > 0,
    ## feature_data
    inherits(feature_data, "data.frame"), ncol(feature_data) > 0,
    nrow(feature_data) > 0,
    ## site_detection_columns
    is.character(site_detection_columns),
    identical(nrow(feature_data), length(site_detection_columns)),
    assertthat::noNA(site_detection_columns),
    all(assertthat::has_name(site_data, site_detection_columns)),
    ## site_n_surveys_columns
    is.character(site_n_surveys_columns),
    identical(nrow(feature_data), length(site_n_surveys_columns)),
    assertthat::noNA(site_n_surveys_columns),
    all(assertthat::has_name(site_data, site_n_surveys_columns)),
    ## site_probability_columns
    is.character(site_probability_columns),
    identical(nrow(feature_data), length(site_probability_columns)),
    assertthat::noNA(site_probability_columns),
    all(assertthat::has_name(site_data, site_probability_columns)),
    ## feature_survey_sensitivity_column
    assertthat::is.string(feature_survey_sensitivity_column),
    all(assertthat::has_name(feature_data, feature_survey_sensitivity_column)),
    is.numeric(feature_data[[feature_survey_sensitivity_column]]),
    assertthat::noNA(
      feature_data[[feature_survey_sensitivity_column]]),
    all(feature_data[[feature_survey_sensitivity_column]] >= 0),
    all(feature_data[[feature_survey_sensitivity_column]] <= 1),
    ## feature_survey_specificity_column
    assertthat::is.string(feature_survey_specificity_column),
    all(assertthat::has_name(feature_data, feature_survey_specificity_column)),
    is.numeric(feature_data[[feature_survey_specificity_column]]),
    assertthat::noNA(feature_data[[feature_survey_specificity_column]]),
    all(feature_data[[feature_survey_specificity_column]] >= 0),
    all(feature_data[[feature_survey_specificity_column]] <= 1),
    ## feature_model_sensitivity_column
    assertthat::is.string(feature_model_sensitivity_column),
    all(assertthat::has_name(feature_data, feature_model_sensitivity_column)),
    is.numeric(feature_data[[feature_model_sensitivity_column]]),
    assertthat::noNA(feature_data[[feature_model_sensitivity_column]]),
    all(feature_data[[feature_model_sensitivity_column]] >= 0),
    all(feature_data[[feature_model_sensitivity_column]] <= 1),
    ## feature_model_specificity_column
    assertthat::is.string(feature_model_specificity_column),
    all(assertthat::has_name(feature_data, feature_model_specificity_column)),
    is.numeric(feature_data[[feature_model_specificity_column]]),
    assertthat::noNA(feature_data[[feature_model_specificity_column]]),
    all(feature_data[[feature_model_specificity_column]] >= 0),
    all(feature_data[[feature_model_specificity_column]] <= 1))
  ## validate survey data
  validate_site_detection_data(site_data, site_detection_columns)
  validate_site_n_surveys_data(site_data, site_n_surveys_columns)

  # drop spatial data
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)

  # extract data
  dij <- t(as.matrix(site_data[, site_detection_columns, drop = FALSE]))
  nij <- t(as.matrix(site_data[, site_n_surveys_columns, drop = FALSE]))
  pij <- t(as.matrix(site_data[, site_probability_columns, drop = FALSE]))
  # calculate prior matrix
  ## initialize matrix
  prior <- matrix(0, ncol = ncol(dij), nrow = nrow(dij))
  ## calculate values for each feature within each site
  for (j in seq_len(ncol(dij))) {
    for (i in seq_len(nrow(dij))) {
      if (nij[i, j] > 0) {
        ### calculate values if has survey data
        prior[i, j] <- prior_probability_of_occupancy_given_survey_data(
          round(dij[i, j] * nij[i, j]), round((1 - dij[i, j]) * nij[i, j]),
          feature_data[[feature_survey_sensitivity_column]][i],
          feature_data[[feature_survey_specificity_column]][i],
          prior = 0.5, clamp = FALSE)
      } else if (pij[i, j] >= 0.5) {
        ### calculate values if no data available, and model predicts presence
        prior[i, j] <- feature_data[[feature_model_sensitivity_column]][i]
      } else {
        ### calculate values if no data available, and model predicts absence
        prior[i, j] <- 1 - feature_data[[feature_model_specificity_column]][i]
      }
    }
  }

  # clamp values to avoid issues with probabilities that are exactly zero or one
  prior[] <- pmin(prior[], 1 - 1e-10)
  prior[] <- pmax(prior[], 1e-10)

  # return result
  rownames(prior) <- site_detection_columns
  prior
}

#' Prior probability of occupancy given survey data
#'
#' Calculate the prior probability of a species occupying a site
#' given that a series of surveys which protected a number of
#' detections and non-detections.
#'
#' @param n_det \code{numeric} number of detections.
#'
#' @param n_nondet \code{numeric} number of non-detections.
#'
#' @param sensitivity \code{numeric} sensitivity of the surveys.
#'
#' @param specificity \code{numeric} specificity of the surveys.
#'
#' @param prior \code{numeric} initial prior probability of occupancy.
#'   Defaults to 0.5.
#'
#' @param clamp \code{logical} should values be clamped to between 1e-10
#'   and 1-1e10 to avoid issues with probabilities that are exactly equal to
#'   zero and one.
#'
#' @return \code{numeric} probability value.
#'
#' @noRd
prior_probability_of_occupancy_given_survey_data <- function(
  n_det, n_nondet, sensitivity, specificity, prior = 0.5, clamp = TRUE) {
  # assert that arguments are valid
  assertthat::assert_that(
    assertthat::is.number(n_det), assertthat::is.number(n_nondet),
    assertthat::is.number(sensitivity), assertthat::is.number(specificity),
    assertthat::is.number(prior),
    all(sensitivity >= 0), all(sensitivity <= 1),
    all(specificity >= 0), all(specificity <= 1),
    all(prior >= 0), all(prior <= 1),
    assertthat::is.flag(clamp), assertthat::noNA(clamp))
  # define local functions, based on case study 1 in
  # https://doi.org/10.1111/2041-210X.12423
  update_prior_prob_pres <- function(p, x, s1, s2) {
    if (x >= 0.5) {
      o <- (p * s1) / ((p * s1) + ((1 - p) * (1 - s2)))
    } else {
      o <- (p * (1 - s1)) / ((p * (1 - s1)) + ((1 - p) * s2))
    }
    o
  }
  # main calculations
  obs <- c(rep(1, n_det), rep(0, n_nondet))
  out <- prior
  for (i in obs) out <- update_prior_prob_pres(out, i, sensitivity, specificity)
  # return result
  if (isTRUE(clamp)) {
    out[out < 1e-10] <- 1e-10
    out[out > (1 - 1e-10)] <- (1 - 1e-10)
  }
  out
}
