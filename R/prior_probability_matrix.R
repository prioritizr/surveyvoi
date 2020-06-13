#' Prior probability matrix
#'
#' Create prior probability matrix for the value of information analysis.
#'
#' @inheritParams evdci
#'
#' @details The prior matrix is constructed using a combination of previous
#'   survey results and modelled predictions for sites that have not been
#'   surveyed. Let \eqn{I} denote the set of features (indexed by \eqn{i}) and
#'   let \eqn{J} denote the set of sites (indexed by \eqn{j}).
#'   Let \eqn{D_j} indicate which sites have already been surveyed (using zeros
#'   and ones).
#'   Let \eqn{S_i} and \eqn{N_i} denote sensitivity and specificity
#'   of the survey method for features \eqn{i \in I} (respectively).
#'   Also let \eqn{H_{ij}} denote the presence or absence of each feature
#'   in each site if they have been surveyed (using zeros or ones),
#'   with unsurveyed sites denoted with -1 values
#'   Let \eqn{{H'}_{ij}} indicate the modelled probability of each feature
#'   occupying each site.
#'   Finally, let \eqn{{S'}_i} and \eqn{{N'}_i} denote
#'   sensitivity and specificity of the model for feature \eqn{i \in I}
#'   (respectively).
#'   Given such data, the (\eqn{P_{ij}}) prior probability of feature \eqn{i}
#'   occupying site \eqn{j} is:
#'
#' \deqn{
#' P_{ij} = \\
#' S_i, \text{ if } D_j = 1, H_{ij} = 1 \space (i \text{ detected in } j) \\
#' 1 - N_i, \text{ else if } D_j = 1, H_{ij} = 0 \space (i \text{ not detected in } j) \\
#' {S'}_i \times S_i, \text{ else if } D_j = 0, {H'}_{ij} \geq 0.5 \space (j \text{ not surveyed and } i \text{ predicted present in } j \text{)} \\
#' 1 - ({N'}_i \times N_i), \text{ else if } D_j = 0, {H'}_{ij} \geq 0.5 \space (j \text{ not surveyed and } i \text{ predicted absent in } j \text{)} \\
#' }
#'
#' Note that the prior probability calculations account for the fact that the
#' species distribution models were evaluated using the survey data. Since
#' the species distribution models were evaluated using survey data, this
#' means that sensitivity and specificity values of the models are
#' conditional on uncertainty present in the survey methodology. To
#' account for this, the prior probability of species \eqn{i} occurring within
#' planning unit \eqn{j} when relying on the species distribution model
#' predictions depends on the sensitivity and specificity of both the
#' species distribution models and the survey methodology.
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
  site_data, feature_data, site_occupancy_columns, site_probability_columns,
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
    ## site_occupancy_columns
    is.character(site_occupancy_columns),
    identical(nrow(feature_data), length(site_occupancy_columns)),
    assertthat::noNA(site_occupancy_columns),
    all(assertthat::has_name(site_data, site_occupancy_columns)),
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
  # drop spatial data
  if (inherits(site_data, "sf"))
    site_data <- st_drop_geometry(site_data)
  # extract data
  rij <- t(as.matrix(site_data[, site_occupancy_columns, drop = FALSE]))
  mij <- t(as.matrix(site_data[, site_probability_columns, drop = FALSE]))
  # calculate prior matrix
  ## initialize matrix
  prior <- matrix(0, ncol = ncol(rij), nrow = nrow(rij))
  ## add in prior probabilities for surveys that detected features
  for (f in seq_len(nrow(rij)))
    prior[f, which(rij[f, ] >= 0.5)] <-
      feature_data[[feature_survey_sensitivity_column]][f]
  ## add in prior probabilities for surveys that did not detect features
  for (f in seq_len(nrow(rij)))
    prior[f, which(rij[f, ] < 0.5)] <-
      1 - feature_data[[feature_survey_specificity_column]][f]
  ## add in prior probabilities for sites units without surveys
  for (f in seq_len(nrow(rij))) {
    ### add model sensitivity for predicted presence,
    ### and 1 - model specificity for predicted absence
    pos <- which(is.na(rij[f, ]))
    prior[f, pos] <-
      ((feature_data[[feature_model_sensitivity_column]][f] *
        feature_data[[feature_survey_sensitivity_column]][f]) *
       (mij[f, pos] >= 0.5)) +
      ((1 - (feature_data[[feature_model_specificity_column]][f] *
             feature_data[[feature_survey_specificity_column]][f])) *
       (mij[f, pos] < 0.5))
  }
  # clamp values that are exactly zero or one
  prior[] <- pmin(prior[], 1 - 1e-10)
  prior[] <- pmax(prior[], 1e-10)
  # return result
  rownames(prior) <- site_occupancy_columns
  prior
}
