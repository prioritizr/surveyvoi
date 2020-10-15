#' @include internal.R
NULL

#' Simulate feature data
#'
#' Simulate feature data for developing simulated survey schemes.
#'
#' @inheritParams simulate_site_data
#'
#' @param proportion_of_survey_features `numeric` proportion of features
#'   that will be examined in the new surveys. Values must be between
#'   zero and one. Defaults to 1 such that all features should be surveyed.
#'
#' @return \code{\link[tibble]{tibble}} object. It contains the following
#' data:
#' \describe{
#' \item{`name`}{`character` name of each feature.}
#'
#' \item{`survey`}{`logical` (`TRUE` / `FALSE`) values
#'   indicating if each feature should be examined in surveys or not.}
#'
#' \item{`survey_sensitivity`}{`numeric` sensitivity (true positive
#'   rate) of the survey methodology for each features.}
#'
#' \item{`survey_specificity`}{`numeric` specificity (true negative
#'   rate) of the survey methodology for each features.}
#'
#' \item{`model_sensitivity`}{`numeric` specificity (true positive
#'   rate) of the occupancy models for each features.}
#'
#' \item{`model_specificity`}{`numeric` specificity (true negative
#'   rate) of the occupancy models for each features.}
#'
#' \item{`target`}{`numeric` target values used to parametrize
#'   the conservation benefit of managing of each feature (defaults to 1).}
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
#' print(d, width = Inf)
#'
#' @export
simulate_feature_data <- function(
  n_features, proportion_of_survey_features = 1) {
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
    survey_sensitivity = stats::runif(n_features, 0.95, 0.99),
    survey_specificity = stats::runif(n_features, 0.8, 0.9),
    model_sensitivity = stats::runif(n_features, 0.7, 0.8),
    model_specificity = stats::runif(n_features, 0.8, 0.9),
    target = 1)
}
