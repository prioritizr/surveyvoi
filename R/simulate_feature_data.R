#' @include internal.R
NULL

#' Simulate feature data
#'
#' Simulate feature data for developing simulated survey schemes.
#'
#' @inheritParams simulate_site_data
#'
#' @param proportion_of_survey_features \code{numeric} proportion of features
#'   that will be examined in the new surveys. Values must be between
#'   zero and one.
#'
#' @return \code{\link[tibble]{tibble}} object. It contains the following
#' data:
#' \describe{
#' \item{\code{name}}{\code{character} name of each feature.}
#'
#' \item{\code{survey}}{\code{logical} (\code{TRUE} / \code{FALSE}) values
#'   indicating if each feature should be examined in surveys or not.}
#'
#' \item{\code{survey_sensitivity}}{\code{numeric} sensitivity (true positive
#'   rate) of the survey methodology for each features.}
#'
#' \item{\code{survey_specificity}}{\code{numeric} specificity (true negative
#'   rate) of the survey methodology for each features.}
#'
#' \item{\code{model_sensitivity}}{\code{numeric} specificity (true positive
#'   rate) of the occupancy models for each features.}
#'
#' \item{\code{alpha}}{\code{numeric} values used to parametrize
#'   the conservation benefit of managing of each feature.}
#'
#' \item{\code{gamma}}{\code{numeric} values used to parametrize
#'   the conservation benefit of managing of each feature.}
#'
#' }
#'
#' @seealso \code{\link{simulate_site_data}}
#'
#' @examples
#' # set seed for reproducibility
#' set.seed(123)
#'
#' # simulate data
#' d <- simulate_feature_data(n_features = 5,
#'                            proportion_of_survey_features = 0.5)
#' # print data
#' print(d)
#'
#' @export
simulate_feature_data <- function(n_features, proportion_of_survey_features) {
  # assert that arguments are valid
  assertthat::assert_that(
    assertthat::is.count(n_features), assertthat::noNA(n_features),
    assertthat::is.number(proportion_of_survey_features),
    assertthat::noNA(proportion_of_survey_features),
    isTRUE(proportion_of_survey_features > 0),
    isTRUE(proportion_of_survey_features <= 1))
  # create vector indicating which features should be surveyed
  survey_feature <- rep(FALSE, n_features)
  pos <- sample.int(n_features, ceiling(proportion_of_survey_features *
                                        n_features))
  survey_feature[pos] <- TRUE
  # return object
  tibble::tibble(
    name = paste0("f", seq_len(n_features)),
    survey = survey_feature,
    survey_sensitivity = runif(n_features, 0.95, 0.99),
    survey_specificity = runif(n_features, 0.8, 0.9),
    model_sensitivity = runif(n_features, 0.7, 0.8),
    alpha = abs(rnorm(n_features)) + 1,
    gamma = runif(n_features))
}
