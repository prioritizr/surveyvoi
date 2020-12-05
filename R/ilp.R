#' @include internal.R
NULL

#' Space based prioritizations
#'
#' This function generates space-based prioritizations by solving the
#' *p*-median problem to optimality.
#'
#' @param x `matrix` object containing a symmetric distance matrix that
#'   measures the dissimilarity between different planning units.
#'
#' @param budget `numeric` vector indicating the maximum expenditure for
#'   each solution. The values must be between one and
#'   the number of planning units (i.e. length of `x`).
#'
#' @param costs `numeric` cost values for each planning units. Defaults
#'   to a vector assigning a cost of one for each planning unit.
#'
#' @param locked_in `logical` vector indicating if each planning unit
#'   is (`TRUE`) locked into the solution or (`FALSE`) not. Defaults
#'   to a vector of `FALSE` for each planning unit.
#'
#' @param locked_out `logical` vector indicating if each planning unit
#'   is (`TRUE`) locked out of the solution or (`FALSE`) not. Defaults
#'   to a vector of `FALSE` for each planning unit.
#'
#' @param verbose `logical` indicating if information should be
#'   printed while generating prioritization. Defaults to `FALSE`.
#'
#' @return `matrix` containing the solution. Each column corresponds
#'   to a different planning unit and each row corresponds to a different
#'   solution. Cell values indicate if a given planning unit was selected in a
#'   given solution.
#'
#' @seealso For more information on the *p*-median problem
#'   see <http://www.hyuan.com/java/how.html>.
#'
#' @noRd
distance_based_prioritizations <- function(x, budget, costs, locked_in,
                                           locked_out, verbose) {
  # assert that data are valid
  assertthat::assert_that(
    ## x
    is.matrix(x), is.numeric(x), all(is.finite(c(x))), ncol(x) > 0, nrow(x) > 0,
    isSymmetric(x),
    ## budget
    all(sapply(budget, assertthat::is.number)), all(budget >= 0),
    ## costs
    is.numeric(costs), assertthat::noNA(costs), all(costs >= 0),
    identical(length(costs), ncol(x)),
    ## locked_in
    is.logical(locked_in), assertthat::noNA(locked_in),
    identical(length(locked_in), ncol(x)),
    ## locked_in
    is.logical(locked_out), assertthat::noNA(locked_out),
    identical(length(locked_out), ncol(x)),
    !any(locked_in & locked_out))
  # prepare data for optimization
  ## initialization
  n <- nrow(x)
  ## create objective function
  obj <- c(rep(0, n), c(x))
  ## create constraint matrix
  A <- rcpp_pmedian_constraint_matrix(x, costs)
  A <- Matrix::sparseMatrix(i = A$i, j = A$j, x = A$x, index1 = FALSE)
  ## create RHS vector
  rhs <- c(NA_real_, rep(1, n), rep(0, length(x)))
  ## create sense vector
  sense <- c("<=", rep("=", n), rep("<=", length(x)))
  ## create vtypes
  vtype <- rep("B", length(obj))
  ## create bounds vectors
  ub <- rep(1, length(obj))
  lb <- rep(0, length(obj))
  lb[which(locked_in)] <- 1 # lock in planning units
  ub[which(locked_out)] <- 0 # lock out planning units
  ## create gurobi input model object
  m <- list(modelsense = "min", obj = obj, sense = sense, rhs = rhs,
            lb = lb, ub = ub, vtype = vtype, A = A)
  ## create gurobi parameters list
  p <- list(Presolve = 2, MIPGap = 0, LogToConsole = as.integer(verbose))
  # prepare output matrix
  out <- matrix(NA, ncol = nrow(x), nrow = length(budget))
  # iterate over each budget and obtain the solution
  for (b in seq_along(budget)) {
    ## update budget
    m$rhs[1] <- budget[b]
    ## solve problem
    withr::with_locale(
      c(LC_CTYPE = "C"), {
      s <- gurobi::gurobi(m, p)$x[seq_len(n)]
    })
    assertthat::assert_that(
      isTRUE(sum(s * costs) <= budget[b]) ||
      isTRUE(abs((sum(s * costs) - budget[b])) <= 1e-5),
      msg = "solver returned infeasible solution")
    ## store solution
    out[b, ] <- as.logical(s)
  }
  # return matrix
  out
}

#' Weight based prioritization
#'
#' This is a convenience function to generate a series of prioritizations that
#' select a set of planning units with the highest sum of weights, subject
#' to a budget.
#'
#' @param x `numeric` vector containing the weights.
#'
#' @param budget `numeric` number indicating how many planning units
#'   should be selected in the solution. The value must be between one and
#'   the number of planning units (i.e. length of `x`).
#'
#' @param locked_in `logical` vector indicating if each planning unit
#'   is (`TRUE`) locked into the solution or (`FALSE`) not. Defaults
#'   to a vector of `FALSE` for each planning unit.
#'
#' @param locked_out `logical` vector indicating if each planning unit
#'   is (`TRUE`) locked out of the solution or (`FALSE`) not. Defaults
#'   to a vector of `FALSE` for each planning unit.
#'
#' @param verbose `logical` indicating if information should be
#'   printed while generating prioritization. Defaults to `FALSE`.
#'
#' @return `matrix` containing the solution. Each column corresponds
#'   to a different planning unit and each row corresponds to a different
#'   row. Cell values indicate if a given planning unit was selected in a
#'   given solution.
#'
#' @noRd
weight_based_prioritizations <- function(x, budget, costs, locked_in,
                                         locked_out, verbose) {
  # assert that data are valid
  assertthat::assert_that(
    ## x
    is.numeric(x), all(is.finite(c(x))), length(x) > 0,
    ## budget
    all(sapply(budget, assertthat::is.number)), all(budget >= 0),
    length(budget) > 0,
    ## costs
    is.numeric(costs), assertthat::noNA(costs), all(costs >= 0),
    identical(length(costs), length(x)),
    ## locked_in
    is.logical(locked_in), assertthat::noNA(locked_in),
    identical(length(locked_in), length(x)),
    ## locked_out
    is.logical(locked_out), assertthat::noNA(locked_out),
    identical(length(locked_out), length(x)),
    !any(locked_in & locked_out))
  assertthat::assert_that(
    all(sum(locked_in * costs) <= budget),
    msg = "locked in sites exceed one or more budgets")
  ## create objective function
  obj <- c(x)
  ## create constraints
  A <- methods::as(matrix(costs, nrow = 1, ncol = length(x)), "dgCMatrix")
  ## create RHS vector
  rhs <- c(NA_real_)
  ## create sense vector
  sense <- "<="
  ## create vtypes
  vtype <- rep("B", length(obj))
  ## create bounds vectors
  ub <- rep(1, length(obj))
  lb <- rep(0, length(obj))
  lb[which(locked_in)] <- 1 # lock in planning units
  ub[which(locked_out)] <- 0 # lock out planning units
  ## create gurobi input model object
  m <- list(modelsense = "max", obj = obj, sense = sense, rhs = rhs,
            lb = lb, ub = ub, vtype = vtype, A = A)
  ## create gurobi parameters list
  p <- list(Presolve = 2, MIPGap = 0, LogToConsole = as.integer(verbose))
  # prepare output matrix
  out <- matrix(NA, ncol = length(x), nrow = length(budget))
  # iterate over each budget and obtain the solution
  for (b in seq_along(budget)) {
    ## update budget
    m$rhs[1] <- budget[b]
    ## solve problem
    withr::with_locale(
      c(LC_CTYPE = "C"), {
      s <- gurobi::gurobi(m, p)$x
    })
    ## assert valid solution
    assertthat::assert_that(
      isTRUE(sum(s * costs) <= budget[b]) ||
      isTRUE(abs((sum(s * costs) - budget[b])) <= 1e-5),
      msg = "solver returned infeasible solution")
    ## store solution
    out[b, ] <- as.logical(s)
  }
  # return matrix
  out
}
