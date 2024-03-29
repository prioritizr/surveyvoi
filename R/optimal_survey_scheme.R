#' @include internal.R evdsi.R
NULL

#' Optimal survey scheme
#'
#' Find the optimal survey scheme that maximizes value of information.
#' This function uses the exact method for
#' calculating the expected value of the decision given a survey scheme.
#'
#' @inheritParams evdsi
#'
#' @param n_threads `integer` number of threads to use for computation.
#'
#' @param survey_budget `numeric` maximum expenditure permitted
#'   for conducting surveys.
#'
#' @param site_survey_locked_out_column `character` name of the column
#'   in the argument to `site_data` that contains `logical`
#'   (`TRUE` / `FALSE`) values indicating which sites should
#'   be locked out (`TRUE`) from being selected for future surveys or
#'   (`FALSE`) not. No missing (`NA`) values are permitted in this
#'   column. This is useful if some sites will never be considered for future
#'   surveys (e.g. because they are too costly to survey, or have a
#'   low chance of containing the target species).
#'   Defaults to `NULL` such that no sites are locked out.
#'
#' @param verbose `logical` indicating if information should be
#'   printed during processing. Defaults to `FALSE`.
#'
#' @details
#' The optimal survey scheme is determined using a brute-force algorithm.
#' Initially, all feasible (valid) survey schemes are identified given the
#' survey costs and the survey budget (using
#' [feasible_survey_schemes()]. Next, the expected value of each and
#' every feasible survey scheme is computed
#' (using [evdsi()]).
#' Finally, the greatest expected value is identified, and all survey schemes
#' that share this greatest expected value are returned. Due to the nature of
#' this algorithm, it can take a very long time to complete.
#'
#' @inheritSection feasible_survey_schemes Dependencies
#'
#' @return A `matrix` of `logical` (`TRUE`/ `FALSE`)
#'   values indicating if a site is selected in the scheme or not. Columns
#'   correspond to sites, and rows correspond to different schemes. If
#'   there is only one optimal survey scheme then the `matrix` will only
#'   contain a single row. This matrix also has a `numeric` `"ev"`
#'   attribute that contains the expected value of each scheme.
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
#' # set total budget for managing sites for conservation
#' # (i.e. 50% of the cost of managing all sites)
#' total_budget <- sum(sim_sites$management_cost) * 0.5
#'
#' # set total budget for surveying sites for conservation
#' # (i.e. 40% of the cost of managing all sites)
#' survey_budget <- sum(sim_sites$survey_cost) * 0.4
#'
#' \dontrun{
#' # find optimal survey scheme using exact method
#' opt_survey <- optimal_survey_scheme(
#'   sim_sites, sim_features,
#'   c("f1", "f2", "f3"), c("n1", "n2", "n3"), c("p1", "p2", "p3"),
#'   "management_cost", "survey_cost",
#'   "survey", "survey_sensitivity", "survey_specificity",
#'   "model_sensitivity", "model_specificity",
#'   "target", total_budget, survey_budget)
#'
#' # print result
#' print(opt_survey)
#' }
#' @export
optimal_survey_scheme <- function(
  site_data, feature_data,
  site_detection_columns, site_n_surveys_columns, site_probability_columns,
  site_management_cost_column,
  site_survey_cost_column,
  feature_survey_column,
  feature_survey_sensitivity_column,
  feature_survey_specificity_column,
  feature_model_sensitivity_column,
  feature_model_specificity_column,
  feature_target_column,
  total_budget,
  survey_budget,
  site_management_locked_in_column = NULL,
  site_management_locked_out_column = NULL,
  site_survey_locked_out_column = NULL,
  prior_matrix = NULL,
  n_threads = 1,
  verbose = FALSE) {
  # assert arguments are valid
  assertthat::assert_that(
    ## site_data
    inherits(site_data, "sf"), ncol(site_data) > 0,
    nrow(site_data) > 0,
    ## feature_data
    inherits(feature_data, "data.frame"), ncol(feature_data) > 0,
    nrow(feature_data) > 0,
    ## site_detection_columns
    is.character(site_detection_columns),
    length(site_detection_columns) > 0,
    assertthat::noNA(site_detection_columns),
    all(assertthat::has_name(site_data, site_detection_columns)),
    length(site_detection_columns) == nrow(feature_data),
    ## site_n_surveys_columns
    is.character(site_n_surveys_columns),
    length(site_n_surveys_columns) > 0,
    assertthat::noNA(site_n_surveys_columns),
    all(assertthat::has_name(site_data, site_n_surveys_columns)),
    length(site_n_surveys_columns) == nrow(feature_data),
    ## site_probability_columns
    is.character(site_probability_columns),
    identical(nrow(feature_data), length(site_probability_columns)),
    assertthat::noNA(site_probability_columns),
    all(assertthat::has_name(site_data, site_probability_columns)),
    ## site_management_cost_column
    assertthat::is.string(site_management_cost_column),
    all(assertthat::has_name(site_data, site_management_cost_column)),
    is.numeric(site_data[[site_management_cost_column]]),
    assertthat::noNA(site_data[[site_management_cost_column]]),
    ## site_survey_cost_column
    assertthat::is.string(site_survey_cost_column),
    all(assertthat::has_name(site_data, site_survey_cost_column)),
    is.numeric(site_data[[site_survey_cost_column]]),
    assertthat::noNA(site_data[[site_survey_cost_column]]),
    ## feature_survey_column
    assertthat::is.string(feature_survey_column),
    all(assertthat::has_name(feature_data, feature_survey_column)),
    is.logical(feature_data[[feature_survey_column]]),
    assertthat::noNA(feature_data[[feature_survey_column]]),
    sum(feature_data[[feature_survey_column]]) >= 1,
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
    all(feature_data[[feature_model_specificity_column]] <= 1),
    ## feature_target_column
    assertthat::is.string(feature_target_column),
    all(assertthat::has_name(feature_data, feature_target_column)),
    is.numeric(feature_data[[feature_target_column]]),
    assertthat::noNA(feature_data[[feature_target_column]]),
    all(feature_data[[feature_target_column]] >= 0),
    ## total_budget
    assertthat::is.number(total_budget), assertthat::noNA(total_budget),
    isTRUE(total_budget > 0),
    ## survey_budget
    assertthat::is.number(survey_budget), assertthat::noNA(survey_budget),
    isTRUE(survey_budget > 0), isTRUE(survey_budget <= total_budget),
    isTRUE(survey_budget >= min(site_data[[site_survey_cost_column]])),
    ## prior_matrix
    inherits(prior_matrix, c("matrix", "NULL")),
    ## verbose
    assertthat::is.flag(verbose),
    assertthat::noNA(verbose),
    ## n_threads
    assertthat::is.count(n_threads),
    assertthat::noNA(n_threads))
  ## site_management_locked_in_column
  if (!is.null(site_management_locked_in_column)) {
    assertthat::assert_that(
      assertthat::is.string(site_management_locked_in_column),
      all(assertthat::has_name(site_data, site_management_locked_in_column)),
      is.logical(site_data[[site_management_locked_in_column]]),
      assertthat::noNA(site_data[[site_management_locked_in_column]]))
    assertthat::assert_that(
      sum(site_data[[site_management_locked_in_column]] *
          site_data[[site_management_cost_column]]) <=
      total_budget,
      msg = "cost of managing locked in sites exceeds total budget")
  }
  ## site_management_locked_out_column
  if (!is.null(site_management_locked_out_column)) {
    assertthat::assert_that(
      assertthat::is.string(site_management_locked_out_column),
      all(assertthat::has_name(site_data, site_management_locked_out_column)),
      is.logical(site_data[[site_management_locked_out_column]]),
      assertthat::noNA(site_data[[site_management_locked_out_column]]))
    if (all(site_data[[site_management_locked_out_column]]))
      warning("all sites locked out")
  }
  ## validate locked arguments if some locked in and some locked out
  if (!is.null(site_management_locked_in_column) &&
      !is.null(site_management_locked_out_column)) {
    assertthat::assert_that(
      all(site_data[[site_management_locked_in_column]] +
          site_data[[site_management_locked_out_column]] <= 1),
      msg = "at least one planning unit is locked in and locked out")
  }
  ## site_survey_locked_out_column
  if (!is.null(site_survey_locked_out_column)) {
    assertthat::assert_that(
      assertthat::is.string(site_survey_locked_out_column),
      all(assertthat::has_name(site_data, site_survey_locked_out_column)),
      is.logical(site_data[[site_survey_locked_out_column]]),
      assertthat::noNA(site_data[[site_survey_locked_out_column]]),
      !all(site_data[[site_survey_locked_out_column]]))
  }
  ## validate targets
  validate_target_data(feature_data, feature_target_column)
  ## validate survey data
  validate_site_detection_data(site_data, site_detection_columns)
  validate_site_n_surveys_data(site_data, site_n_surveys_columns)
  ## validate model probability values
  validate_site_probability_data(site_data, site_probability_columns)
  ## verify targets
  assertthat::assert_that(
    all(feature_data[[feature_target_column]] <= nrow(site_data)))
  if (!is.null(site_management_locked_out_column)) {
    assertthat::assert_that(
      all(feature_data[[feature_target_column]] <=
          sum(!site_data[[site_management_locked_out_column]])))
  }

  # prepare data for analysis
  ## drop spatial information
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)
  ## calculate prior matrix
  if (is.null(prior_matrix)) {
    pij <- prior_probability_matrix(
      site_data, feature_data, site_detection_columns,
      site_n_surveys_columns, site_probability_columns,
      feature_survey_sensitivity_column, feature_survey_specificity_column,
      feature_model_sensitivity_column, feature_model_specificity_column)
  } else {
    validate_prior_data(prior_matrix, nrow(site_data), nrow(feature_data))
    pij <- prior_matrix
  }
  ## prepare site management locked in data
  if (!is.null(site_management_locked_in_column)) {
    site_management_locked_in <- site_data[[site_management_locked_in_column]]
  } else {
    site_management_locked_in <- rep(FALSE, nrow(site_data))
  }
  ## prepare site survey locked out data
  if (!is.null(site_survey_locked_out_column)) {
    site_survey_locked_out <- site_data[[site_survey_locked_out_column]]
  } else {
    site_survey_locked_out <- rep(FALSE, nrow(site_data))
  }
  ## prepare locked out data
  if (!is.null(site_management_locked_out_column)) {
    site_management_locked_out <- site_data[[site_management_locked_out_column]]
  } else {
    site_management_locked_out <- rep(FALSE, nrow(site_data))
  }
  ## validate that targets are feasible given budget and locked out units
  sorted_costs <- sort(
    site_data[[site_management_cost_column]][!site_management_locked_out])
  sorted_costs <- sorted_costs[
    seq_len(max(feature_data[[feature_target_column]]))]
  assertthat::assert_that(
    sum(sorted_costs) <= total_budget,
    msg = paste("targets cannot be achieved given budget and locked out",
                "planning units"))
  ## extract site data
  nij <- t(as.matrix(site_data[, site_n_surveys_columns, drop = FALSE]))
  ## identify planning units that have been surveyed for all species
  site_survey_status <- colSums(nij < 0.5) == 0

  # generate survey schemes
  all_feasible_schemes <- feasible_survey_schemes(
    data.frame(x = site_data[[site_survey_cost_column]],
               out = site_survey_status | site_survey_locked_out),
    "x", survey_budget, NULL, "out", verbose = FALSE)

  # identify which schemes involve collecting new data vs. using existing data
  new_info_idx <- which(rowSums(all_feasible_schemes) > 0.5)
  current_idx <- which(rowSums(all_feasible_schemes) < 0.5)

  # calculate expected value of decision given scheme that does not survey sites
  evd_current <- rcpp_expected_value_of_decision_given_current_info(
    pij = pij,
    pu_costs = site_data[[site_management_cost_column]],
    pu_locked_in = site_management_locked_in,
    pu_locked_out = site_management_locked_out,
    target = round(feature_data[[feature_target_column]]),
    budget = total_budget)

  # calculate expected value of decision given schemes that survey sites
  ## initialize cluster
  if (n_threads > 1) {
    cl <- start_cluster(n_threads,
      c("pij", "new_info_idx", "all_feasible_schemes",
        "site_data", "feature_data",
        "site_management_locked_in",
        "site_management_locked_out",
        "feature_survey_column",
        "feature_survey_sensitivity_column",
        "feature_survey_specificity_column",
        "site_survey_cost_column",
        "site_management_cost_column",
        "feature_target_column",
        "total_budget",
        "rcpp_expected_value_of_decision_given_survey_scheme"))
    on.exit(try(stop_cluster(cl), silent = TRUE), add = TRUE)
  }
  ## run calculations
  evd_new_info <- plyr::laply(
    seq_along(new_info_idx),
    .parallel = n_threads > 1,
    .progress = ifelse(verbose && (n_threads == 1), "text", "none"),
    .paropts = list(.packages = "surveyvoi"),
    function(i) {
    ## extract i'th survey scheme
    site_survey_scheme <- all_feasible_schemes[new_info_idx[i], , drop = TRUE]
    ## calculate expected value of decision given survey scheme
    rcpp_expected_value_of_decision_given_survey_scheme(
      pij = pij,
      survey_features = feature_data[[feature_survey_column]],
      survey_sensitivity = feature_data[[feature_survey_sensitivity_column]],
      survey_specificity = feature_data[[feature_survey_specificity_column]],
      pu_survey_solution = site_survey_scheme,
      pu_survey_costs = site_data[[site_survey_cost_column]],
      pu_purchase_costs = site_data[[site_management_cost_column]],
      pu_purchase_locked_in = site_management_locked_in,
      pu_purchase_locked_out = site_management_locked_out,
      obj_fun_target = round(feature_data[[feature_target_column]]),
      total_budget = total_budget)
  })
  ## kill cluster
  if (n_threads > 1) {
    cl <- stop_cluster(cl)
  }

  # assemble expected values
  evd <- numeric(nrow(all_feasible_schemes))
  evd[new_info_idx] <- evd_new_info
  evd[current_idx] <- evd_current

  # find optimal solution(s)
  optimal_idx <- abs(max(evd) - evd) < 1e-15
  out <- all_feasible_schemes[optimal_idx, , drop = FALSE]
  attr(out, "ev") <- evd[optimal_idx]

  # return result
  out
}
