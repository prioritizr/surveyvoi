#' Prior probability matrix
#'
#' Create prior probability matrix for the value of information analysis.
#'
#' @inheritParams expected_value_of_survey_scheme
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
#'   occupying each site. Finally, let \eqn{{S'}_i} and \eqn{{N'}_i} denote
#'   sensitivity and specificity of the model for feature \eqn{i \in I}.
#'   Given such data, the (\eqn{P_{ij}}) prior probability of feature \eqn{i}
#'   occupying site \eqn{j} is:
#'
#' \deqn{
#' P_{ij} = \\
#' S_i, \text{ if } D_j = 1, H_{ij} = 1 \space (i \text{ detected in } j) \\
#' 1 - N_i, \text{ else if } D_j = 1, H_{ij} = 0 \space (i \text{ not detected in } j) \\
#' {S'}_i {H'}_{ij}, \text{ else if } D_j = 0 \space (j \text{ not surveyed)} \\
#' }
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
#' # calculate prior probability matrix
#' prior_probability_matrix(site_data, feature_data, c("f1", "f2"),
#'                          c("p1", "p2"), "sensitivity", "specificity")
#'
#' @export
prior_probability_matrix <- function(
  site_data, feature_data, site_occupancy_columns, site_probability_columns,
  feature_survey_sensitivity_column, feature_survey_specificity_column,
  feature_model_sensitivity_column) {
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
    all(feature_data[[feature_model_sensitivity_column]] <= 1))
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
    ### add probabilities for modelled presences
    pos <- which(is.na(rij[f, ]))
    prior[f, pos] <-
      feature_data[[feature_model_sensitivity_column]][f] * mij[f, pos]
  }
  # return result
  rownames(prior) <- site_occupancy_columns
  prior
}
