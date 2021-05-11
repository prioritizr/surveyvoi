#' @include internal.R evdsi.R
NULL

#' Approximately near optimal survey scheme
#'
#' Find a near optimal survey scheme that maximizes value of information.
#' This function uses the approximation method
#' for calculating the expected value of the decision given a survey scheme,
#' and a greedy heuristic algorithm to maximize this metric.
#'
#' @inheritParams approx_evdsi
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
#' Ideally, the brute-force algorithm would be used to identify the optimal
#' survey scheme. Unfortunately, it is not feasible to apply the brute-force
#' to large problems because it can take an incredibly long time to complete.
#' In such cases, it may be desirable to obtain a "relatively good" survey
#' scheme and the greedy heuristic algorithm is provided for such cases.
#' The greedy heuristic algorithm -- unlike the brute force algorithm --
#' is not guaranteed to identify an optimal solution -- or even a "relatively
#' good solution" for that matter -- though greedy heuristic algorithms tend to
#' deliver solutions that are 15\% from optimality. Specifically, this
#' greedy algorithms is implemented as:
#'
#' \enumerate{
#'
#' \item Initialize an empty *list of survey scheme solutions*, and an
#' empty *list of approximate expected values*.
#'
#' \item Calculate the expected value of current information.
#'
#' \item Add a survey scheme with no sites selected for surveying to the
#' *list of survey scheme solutions*, and add the expected value of current
#' information to the *list of approximate expected values*.
#'
#' \item Set the *current survey solution* as the survey scheme with no
#' sites selected for surveying.
#'
#' \item For each remaining candidate site that has not been selected for
#' a survey, generate a new candidate survey scheme with each candidate site
#' added to the current survey solution.
#'
#' \item Calculate the approximate expected value of each
#' new candidate survey scheme. If the cost of a given candidate survey scheme
#' exceeds the survey budget, then store a missing `NA value` instead.
#' Also if the the cost of a given candidate survey scheme plus the
#' management costs of locked in planning units exceeds the total budget,
#' then a store a missing value `NA` value too.
#'
#' \item If all of the new candidate survey schemes are associated with
#' missing `NA` values -- because they all exceed the survey budget -- then
#' go to step 12.
#'
#' \item Calculate the cost effectiveness of each new candidate survey
#' scheme. This calculated as the difference between the approximate expected
#' value of a given new candidate survey scheme and that of the
#' *current survey solution*, and dividing this difference by the the cost
#' of the newly selected candidate site.
#'
#' \item Find the new candidate survey scheme that is associated with the
#' highest cost-effectiveness value, ignoring any missing `NA` values.
#' This new candidate survey scheme is now set as the
#' *current survey scheme*.
#'
#' \item Store the *current survey scheme* in the
#' *list of survey scheme solutions* and store its approximate expected
#' value in the *list of approximate expected values*.
#'
#' \item Go to step 12.
#'
#' \item Find the solution in the *list of survey scheme solutions* that
#' has the highest expected value in the
#' *list of approximate expected values* and return this solution.
#'
#' }
#'
#' @return
#' `matrix` of `logical` (`TRUE`/ `FALSE`)
#' values indicating if a site is selected in the scheme or not. Columns
#' correspond to sites, and rows correspond to different schemes. If there
#' are no ties for the best identified solution, then the the `matrix`
#' will only contain a single row.
#'
#' @examples
#' \dontrun{
#' # set seeds for reproducibility
#' library(RandomFields)
#' set.seed(123)
#' RFoptions(seed = 201)
#'
#' # simulate data
#' site_data <- simulate_site_data(n_sites = 30, n_features = 2, prop = 0.1)
#' feature_data <- simulate_feature_data(n_features = 2, prop = 1)
#' feature_data$target <- c(10, 10)
#'
#' # preview simulated data
#' print(site_data)
#' print(feature_data)
#'
#' # set total budget for managing sites for conservation
#' # (i.e. 50% of the cost of managing all sites)
#' total_budget <- sum(site_data$management_cost) * 0.5
#'
#' # set total budget for surveying sites for conservation
#' # (i.e. 10% of the cost of managing all sites)
#' survey_budget <- sum(site_data$survey_cost) * 0.1
#'
#' # find survey scheme using approximate method and greedy heuristic algorithm
#' # (using 10 replicates so that this example completes relatively quickly)
#' approx_near_optimal_survey <- approx_near_optimal_survey_scheme(
#'   site_data, feature_data,
#'   c("f1", "f2"), c("n1", "n2"), c("p1", "p2"),
#'   "management_cost", "survey_cost",
#'   "survey", "survey_sensitivity", "survey_specificity",
#'   "model_sensitivity", "model_specificity",
#'   "target", total_budget, survey_budget)
#'
#' # print result
#' print(approx_near_optimal_survey)
#' }
#' @export
approx_near_optimal_survey_scheme <- function(
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
  n_approx_replicates = 100,
  n_approx_outcomes_per_replicate = 10000,
  seed = 500,
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
    ## n_threads
    assertthat::is.count(n_threads),
    assertthat::noNA(n_threads),
    ## n_approx_replicates
    assertthat::is.count(n_approx_replicates),
    assertthat::noNA(n_approx_replicates),
    ## n_approx_outcomes_per_replicate
    assertthat::is.count(n_approx_outcomes_per_replicate),
    assertthat::noNA(n_approx_outcomes_per_replicate),
    ## seed
    assertthat::is.number(seed),
    ## verbose
    assertthat::is.flag(verbose),
    assertthat::noNA(verbose))
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
  ## prepare locked out data
  if (!is.null(site_management_locked_out_column)) {
    site_management_locked_out <- site_data[[site_management_locked_out_column]]
  } else {
    site_management_locked_out <- rep(FALSE, nrow(site_data))
  }
  ## prepare site survey locked out data
  if (!is.null(site_survey_locked_out_column)) {
    site_survey_locked_out <- site_data[[site_survey_locked_out_column]]
  } else {
    site_survey_locked_out <- rep(FALSE, nrow(site_data))
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

  # calculate expected value of decision given scheme that does not survey sites
  evd_current <- withr::with_seed(seed, {
    rcpp_expected_value_of_decision_given_current_info(
      pij = pij,
      pu_costs = site_data[[site_management_cost_column]],
      pu_locked_in = site_management_locked_in,
      pu_locked_out = site_management_locked_out,
      target = round(feature_data[[feature_target_column]]),
      budget = total_budget)
  })

  # initialize looping variables
  candidate_sites <- !(site_survey_status | site_survey_locked_out)
  n_candidate_sites <- sum(candidate_sites)
  survey_solution_matrix <-
    matrix(FALSE, ncol = nrow(site_data), nrow = n_candidate_sites + 1)
  survey_solution_values <- numeric(n_candidate_sites + 1)
  survey_solution_values[1] <- mean(evd_current)

  # initialize progress bar
  if (isTRUE(verbose)) {
    pb <- progress::progress_bar$new(
      format = "  optimizing [:bar] :percent eta: :eta",
      total = n_candidate_sites + 1, clear = FALSE, width = 60,
      show_after = 0, force = TRUE)
    pb$tick(0)
  }

  # iterate over the the total number of available sites
  for (s in (1 + seq_len(n_candidate_sites))) {
    # extract previous solution
    prev_solution <- survey_solution_matrix[s - 1, ]
    # find candidate remaining sites that have not been selected yet
    curr_remaining_sites <- which(!prev_solution & candidate_sites)
    # calculate approx expected value of survey information when adding
    # each candidate remaining site to the previous solution
    ## initialize cluster
    if (n_threads > 1) {
      cl <- start_cluster(n_threads,
        c("pij", "prev_solution", "survey_budget", "seed",
          "site_data", "feature_data",
          "site_management_cost_column",
          "site_management_locked_in",
          "site_survey_cost_column",
          "feature_survey_column",
          "feature_survey_sensitivity_column",
          "site_management_locked_out",
          "feature_target_column",
          "total_budget",
          "n_approx_replicates",
          "n_approx_replicates",
          "n_approx_outcomes_per_replicate",
          "rcpp_approx_expected_value_of_decision_given_survey_scheme"))
    }
    ## run calculations
    curr_sites_approx_evsdi <- plyr::laply(
      curr_remaining_sites,
      .parallel = n_threads > 1,
      .paropts = list(.packages = "surveyvoi"),
      function(j) {
      ## generate solution
      curr_candidate_solution <- prev_solution
      curr_candidate_solution[j] <- TRUE
      ## cost of surveys exceeds survey budget then return NA
      curr_surv_cost <-
        sum(site_data[[site_survey_cost_column]] * curr_candidate_solution)
      if (curr_surv_cost > survey_budget) return(NA_real_)
      ## calculate cost of solution, if it exceeds the budget then return NA
      curr_total_cost <-
        curr_surv_cost +
        sum(site_data[[site_management_cost_column]] *
            site_management_locked_in)
      if (curr_total_cost > total_budget) return(NA_real_)
      ## calculate expected value of decision given survey scheme
      out <- withr::with_seed(seed, {
        rcpp_approx_expected_value_of_decision_given_survey_scheme(
          pij = pij,
          survey_features = feature_data[[feature_survey_column]],
          survey_sensitivity =
            feature_data[[feature_survey_sensitivity_column]],
          survey_specificity =
            feature_data[[feature_survey_specificity_column]],
          pu_survey_solution = curr_candidate_solution,
          pu_survey_costs = site_data[[site_survey_cost_column]],
          pu_purchase_costs = site_data[[site_management_cost_column]],
          pu_purchase_locked_in = site_management_locked_in,
          pu_purchase_locked_out = site_management_locked_out,
          obj_fun_target = round(feature_data[[feature_target_column]]),
          total_budget = total_budget,
          n_approx_replicates = n_approx_replicates,
          n_approx_outcomes_per_replicate = n_approx_outcomes_per_replicate,
          seed = seed)
      })
      ## return average expected value
      mean(out)
    })

    ## kill cluster
    if (n_threads > 1) {
      cl <- stop_cluster(cl)
    }

    # update progress bar
    if (isTRUE(verbose)) {
      pb$tick()
    }

    # check to see if main loop should be exited
    ## if all candidate solutions exceed the budget then exit loop
    if (all(is.na(curr_sites_approx_evsdi))) break()

    # penalise each objective value by the cost of the extra planning unit
    if (all(curr_sites_approx_evsdi < survey_solution_values[s])) {
      curr_eval_metrics <-
        (1 / (survey_solution_values[s] - curr_sites_approx_evsdi)) /
        site_data[[site_survey_cost_column]][curr_remaining_sites]
    } else {
      curr_eval_metrics <-
        (curr_sites_approx_evsdi - survey_solution_values[s]) /
        site_data[[site_survey_cost_column]][curr_remaining_sites]
    }

    # find the best site
    curr_best_idx <- which.max(curr_eval_metrics)
    curr_best_site <- curr_remaining_sites[curr_best_idx]

    # store the objective value
    survey_solution_values[s] <- curr_sites_approx_evsdi[curr_best_idx]

    # add the best site to the solution matrix
    survey_solution_matrix[s, ] <- prev_solution
    survey_solution_matrix[s, curr_best_site] <- TRUE
  }

  # find optimal solution(s)
  best_idx <- abs(max(survey_solution_values) - survey_solution_values) < 1e-15
  out <- survey_solution_matrix[best_idx,  , drop = FALSE]
  attr(out, "ev") <- survey_solution_values[best_idx]

  # return result
  out
}
