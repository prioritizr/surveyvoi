#' @include internal.R
NULL

#' Find all feasible survey schemes
#'
#' Generate a `matrix` representing all possible different
#' survey schemes given survey costs and a fixed budget.
#' Please note that this function requires the Gurobi optimization software
#' (<https://www.gurobi.com/>) and the \pkg{gurobi} R package
#' (installation instructions available for [Linux](https://www.gurobi.com/documentation/9.1/quickstart_linux/r_ins_the_r_package.html), [Windows](https://www.gurobi.com/documentation/9.1/quickstart_windows/r_ins_the_r_package.html), and [Mac OS](https://www.gurobi.com/documentation/9.1/quickstart_mac/r_ins_the_r_package.html)).
#'
#' @inheritParams env_div_survey_scheme
#'
#' @param survey_budget `numeric` the maximum possible expenditure
#'  permitted for conducting surveys.
#'
#' @param verbose `logical` indicating if information should be
#'   printed while searching for feasible schemes. Defaults to `FALSE`.
#'
#' @return [matrix()] where each row corresponds to a different
#'   survey scheme, and each column corresponds to a different planning unit.
#'   Cell values are `logical` (`TRUE` / `FALSE`) indicating
#'   if a given site is selected in a given survey scheme.
#'
#' @examples
#' \dontrun{
#' # set seed for reproducibility
#' set.seed(123)
#'
#' # simulate data
#' x <- sf::st_as_sf(tibble::tibble(x = rnorm(4), y = rnorm(4),
#'                                  cost = c(100, 200, 0.2, 1)),
#'                   coords = c("x", "y"))
#'
#' # print data
#' print(x)
#'
#' # plot site locations
#' plot(st_geometry(x), pch = 16, cex = 3)
#'
#' # generate all feasible schemes given a budget of 4
#' s <- feasible_survey_schemes(x, "cost", survey_budget = 4)
#'
#' # print schemes
#' print(s)
#'
#' # plot first scheme
#' x$scheme_1 <- s[1, ]
#' plot(x[, "scheme_1"], pch = 16, cex = 3)
#' }
#' @export
feasible_survey_schemes <- function(
  site_data, cost_column, survey_budget, locked_in_column = NULL,
  locked_out_column = NULL, verbose = FALSE) {
  # assert that arguments are valid
  assertthat::assert_that(
    ## site_data
    inherits(site_data, c("sf", "data.frame")),
    nrow(site_data) > 0, ncol(site_data) > 0,
    ## cost_column
    assertthat::is.string(cost_column), assertthat::noNA(cost_column),
    all(assertthat::has_name(site_data, cost_column)),
    is.numeric(site_data[[cost_column]]),
    assertthat::noNA(site_data[[cost_column]]),
    all(site_data[[cost_column]] >= 0),
    ## survey_budget
    is.numeric(survey_budget), assertthat::noNA(survey_budget),
    all(survey_budget >= 0),
    ## verbose
    assertthat::is.flag(verbose), assertthat::noNA(verbose))
  if (!is.null(locked_in_column)) {
    ## locked_in_column
    assertthat::assert_that(
      assertthat::is.string(locked_in_column),
      all(assertthat::has_name(site_data, locked_in_column)),
      is.logical(site_data[[locked_in_column]]),
      assertthat::noNA(site_data[[locked_in_column]]))
  }
  if (!is.null(locked_out_column)) {
    ## locked_out_column
    assertthat::assert_that(
      assertthat::is.string(locked_out_column),
      all(assertthat::has_name(site_data, locked_out_column)),
      is.logical(site_data[[locked_out_column]]),
      assertthat::noNA(site_data[[locked_out_column]]))
  }

  # set locked in sites
  locked_in <- rep(FALSE, nrow(site_data))
  if (!is.null(locked_in_column))
    locked_in <- site_data[[locked_in_column]]

  # set locked out sites
  locked_out <- rep(FALSE, nrow(site_data))
  if (!is.null(locked_out_column))
    locked_out <- site_data[[locked_out_column]]

  # generate feasible survey scheme
  if (abs(diff(range(site_data[[cost_column]]))) < 1e-10) {
    ## if all the sites have identical costs, no sites are locked in, and
    ## no sites are locked out, then generate the schemes using combinations
    out <- manual_feasible_survey_schemes(
      site_data[[cost_column]], survey_budget, locked_in, locked_out, verbose)
  } else {
    ## otherwise, if we have complicated criteria for generating the schemes,
    ## then use ILP to generate the schemes
    out <- gurobi_feasible_survey_schemes(
      site_data[[cost_column]], survey_budget, locked_in, locked_out, verbose)
  }

  # return result
  out
}

manual_feasible_survey_schemes <- function(cost, budget, locked_in,
                                           locked_out, verbose) {
  # initial calculations
  m <- length(cost)
  # find out which sites are not fixed in the schemes
  idx <- which(!locked_in & !locked_out)
  # if all sites are fixed then return a matrix with just the locked in sites
  if (length(idx) == 0)
    return(matrix(locked_in, nrow = 1))
  # if all sites except one are fixed then compute solution
  if (length(idx) == 1) {
    out <- matrix(FALSE, nrow = 2, ncol = m)
    out[, which(locked_in)] <- TRUE
    out[2, idx] <- TRUE
    return(out)
  }
  # find number of combinations
  remaining_budget <- budget - sum(cost[locked_in])
  k <- min(floor(remaining_budget / max(cost)), length(idx))
  # generate candidate schemes
  out <- plyr::llply(seq_len(k), .progress = ifelse(verbose, "text", "none"),
                     function(i) {
    cmb <- RcppAlgos::comboGeneral(idx, i)
    out <- matrix(FALSE, ncol = m, nrow = nrow(cmb))
    ind <- matrix(c(rep(seq_len(nrow(cmb)), times = i), c(cmb)), ncol = 2)
    out[ind] <- TRUE
    out
  })
  out <- do.call(rbind, out)
  # lock in sites
  out[, which(locked_in)] <- TRUE
  # add scheme with no sites selected if none locked in
  if (sum(locked_in) == 0)
    out <- rbind(FALSE, out)
  # remove duplicate schemes
  out <- as.matrix(dplyr::distinct(as.data.frame(out)))
  attr(out, "dimnames") <- NULL
  # return schemes
  out
}

gurobi_feasible_survey_schemes <- function(cost, budget, locked_in,
                                           locked_out, verbose) {

  # assign locked in status for ilp variables
  locked_in2 <- rep(FALSE, length(cost) * 2)
  raw <- matrix(FALSE, ncol = length(cost), nrow = 2)
  raw[1, which(locked_in)] <- TRUE
  pos <- which(raw)
  locked_in2[pos] <- TRUE

    # assign locked out status for ilp variables
  locked_out2 <- rep(FALSE, length(cost) * 2)
  raw <- matrix(FALSE, ncol = length(cost), nrow = 2)
  raw[1, which(locked_out)] <- TRUE
  pos <- which(raw)
  locked_out2[pos] <- TRUE

  # construct gurobi model matrix
  costs <- matrix(c(cost, rep(0, length(cost))),
                  nrow = 2, ncol = length(cost), byrow = TRUE)
  m <- list(
    obj = rep(1, length(costs)),
    A = methods::as(rcpp_feasible_actions_ilp_matrix(costs), "dgCMatrix"),
    sense = c(rep("=", ncol(costs)), "<="),
    rhs = c(rep(1, ncol(costs)), budget),
    lb = rep(0, length(costs)),
    ub = rep(1, length(costs)),
    vtype = rep("B", length(costs)), modelsense = "max")

  # set locked in and lock out sites
  m$lb[locked_in2] <- 1
  m$ub[locked_out2] <- 0
  assertthat::assert_that(
    all(m$ub >= m$lb),
    msg = "some sites are locked in and locked out")

  # generate solutions
  withr::with_locale(
    c(LC_CTYPE = "C"), {
    g <- suppressMessages(gurobi::gurobi(m, list(
      LogToConsole = as.integer(verbose), Presolve = 2, PoolSearchMode = 2,
      PoolSolutions = 1e+100)))
  })

  # verify that solutions are found
  if (identical(g$status, "INFEASIBLE"))
    stop("failed to generate solutions, please contact the package maintainer")

  # extract solutions
  s <- lapply(g$pool, `[[`, "xn")
  out <- vapply(s, function(x) max.col(matrix(x, byrow = TRUE,
                                              nrow = ncol(costs))),
              numeric(ncol(costs)))
  if (!is.matrix(out)) {
    out <- matrix(out, ncol = ncol(costs))
  } else {
    out <- t(out)
  }

  # convert values to logical
  out <- matrix(out == 1, ncol = ncol(out), nrow = nrow(out))

  # return out
  out
}
