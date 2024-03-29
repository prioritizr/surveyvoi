#' Prior probability matrix
#'
#' Create prior probability matrix for the value of information analysis.
#'
#' @inheritParams evdci
#'
#' @return A `matrix` object containing the prior probabilities of each
#'   feature occupying each site. Each row corresponds to a different
#'   feature and each column corresponds to a different site.
#'
#' @examples
#' # set seeds for reproducibility
#' set.seed(123)
#'
#' # load example site data
#' data(sim_sites)
#' print(sim_sites)
#'
#' # load example feature data
#' data(sim_features)
#' print(sim_features)
#'
#' # calculate prior probability matrix
#' prior_matrix <- prior_probability_matrix(
#'   sim_sites, sim_features,
#'   c("f1", "f2", "f3"), c("n1", "n2", "n3"), c("p1", "p2", "p3"),
#'   "survey_sensitivity", "survey_specificity",
#'   "model_sensitivity", "model_specificity")
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
  validate_site_detection_data(
    site_data, site_detection_columns, check_zeros = FALSE)
  validate_site_n_surveys_data(
    site_data, site_n_surveys_columns, check_zeros = FALSE)

  # calculate prior matrix
  prior <- internal_prior_probability_matrix(
    site_data, feature_data,
    site_detection_columns, site_n_surveys_columns,
    site_probability_columns,
    feature_survey_sensitivity_column, feature_survey_specificity_column,
    feature_model_sensitivity_column, feature_model_specificity_column,
    prefer_survey_data = FALSE)

  # return result
  rownames(prior) <- site_detection_columns
  prior
}


#' Internal prior probability matrix
#'
#' Internal function for generating a prior probability matrix.
#'
#' @inheritParams prior_probability_matrix
#'
#' @param prefer_survey_data `logical` should survey data be used
#'  preferentially instead of model predictions, even if the species
#'  distribution models outperform the survey methodology?
#'  Defaults to `FALSE`
#'
#' @inherit prior_probability_matrix return
#'
#' @noRd
internal_prior_probability_matrix <- function(
  site_data, feature_data,
  site_detection_columns, site_n_surveys_columns,
  site_probability_columns,
  feature_survey_sensitivity_column, feature_survey_specificity_column,
  feature_model_sensitivity_column, feature_model_specificity_column,
  prefer_survey_data = FALSE) {
  # drop spatial data
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)
  # extract data
  dij <- t(as.matrix(site_data[, site_detection_columns, drop = FALSE]))
  nij <- t(as.matrix(site_data[, site_n_surveys_columns, drop = FALSE]))
  pij <- t(as.matrix(site_data[, site_probability_columns, drop = FALSE]))
  survey_sensitivity <- feature_data[[feature_survey_sensitivity_column]]
  survey_specificity <- feature_data[[feature_survey_specificity_column]]
  model_sensitivity <- feature_data[[feature_model_sensitivity_column]]
  model_specificity <- feature_data[[feature_model_specificity_column]]
  # initialize matrix
  prior <- matrix(0, ncol = ncol(dij), nrow = nrow(dij))
  # determine if each site has survey data for each feature
  is_site_have_survey_data <- nij > 0.5
  # calculate model performance
  perf_model_tss <- model_sensitivity + model_specificity - 1
  # calculate survey methodology performance
  ## accounting for the number of surveys within planning unit j
  ## note we assume that each repeat survey is independent, so
  ## the sensitivity and specificities follow Bernoulli distribution
  perf_survey_tss <- dij
  perf_survey_tss[] <- -Inf
  idx <- which(is_site_have_survey_data, arr.ind = TRUE)
  for (ii in seq_len(nrow(idx))) {
    i <- idx[ii, 1]
    j <- idx[ii, 2]
    perf_survey_tss[i, j] <-
      (1 - prod(1 - rep(survey_sensitivity[i], nij[i, j]))) +
      (1 - prod(1 - rep(survey_specificity[i], nij[i, j]))) - 1
  }

  # calculate values for each feature within each site
  for (j in seq_len(ncol(dij))) {
    for (i in seq_len(nrow(dij))) {
      ## determine if the survey data provide a better understanding
      ## of whether feature i is in planning unit j
      is_survey_data_better_than_model <-
        perf_survey_tss[i, j] >= perf_model_tss[i]
      ## now calculate the posterior probabilities...
      if (is_site_have_survey_data[i, j] &&
          (is_survey_data_better_than_model || prefer_survey_data)) {
        ## calculate values if has survey data
        prior[i, j] <-
          prior_probability_of_occupancy(
            round(dij[i, j] * nij[i, j]), round((1 - dij[i, j]) * nij[i, j]),
            survey_sensitivity[i], survey_specificity[i],
            prior = 0.5, clamp = FALSE)
      } else if (pij[i, j] >= 0.5) {
        ## calculate values if no data available, and model predicts presence
        prior[i, j] <-
          prior_probability_of_occupancy(
            1, 0, model_sensitivity[i], model_specificity[i],
            prior = 0.5, clamp = FALSE)
      } else {
        ## calculate values if no data available, and model predicts absence
        prior[i, j] <-
          prior_probability_of_occupancy(
            0, 1, model_sensitivity[i], model_specificity[i],
            prior = 0.5, clamp = FALSE)
      }
    }
  }

  # clamp values to avoid issues with probabilities that are exactly zero or one
  prior[] <- pmin(prior[], 1 - 1e-10)
  prior[] <- pmax(prior[], 1e-10)

  # return result
  prior
}

#' Prior probability of occupancy
#'
#' Calculate the prior probability of a species occupying a site
#' given that a series of surveys which protected a number of
#' detections and non-detections.
#'
#' @param n_det `numeric` number of detections.
#'
#' @param n_nondet `numeric` number of non-detections.
#'
#' @param sensitivity `numeric` sensitivity of the surveys.
#'
#' @param specificity `numeric` specificity of the surveys.
#'
#' @param prior `numeric` initial prior probability of occupancy.
#'   Defaults to 0.5.
#'
#' @param clamp `logical` should values be clamped to between 1e-10
#'   and 1-1e10 to avoid issues with probabilities that are exactly equal to
#'   zero and one.
#'
#' @return A `numeric` probability value.
#'
#' @noRd
prior_probability_of_occupancy <- function(
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
