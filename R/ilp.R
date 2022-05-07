#' @include internal.R
NULL

#' Distance based prioritizations
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
#' @param costs `numeric` cost values for each planning units.
#'
#' @param locked_in `logical` vector indicating if each planning unit
#'   is (`TRUE`) locked into the solution or (`FALSE`) not.
#'
#' @param locked_out `logical` vector indicating if each planning unit
#'   is (`TRUE`) locked out of the solution or (`FALSE`) not.
#'
#' @param solver `character` name of the optimization solver to use
#'   for generating survey schemes.
#'   Available options include: `"Rsymphony"`, `"gurobi"` and `"auto"`.
#'   The `"auto"` method will use the Gurobi optimization software if
#'   it is available; otherwise, it will use the SYMPHONY software
#'   via the \pkg{Rsymphony} package.
#'
#' @param verbose `logical` indicating if information should be
#'   printed while generating prioritization.
#'
#' @return A `matrix` containing the solution. Each column corresponds
#'   to a different planning unit and each row corresponds to a different
#'   solution. Cell values indicate if a given planning unit was selected in a
#'   given solution.
#'
#' @seealso For more information on the *p*-median problem
#'   see <http://www.hyuan.com/java/how.html>.
#'
#' @noRd
distance_based_prioritizations <- function(
  x, budget, costs, locked_in, locked_out, solver, verbose) {
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
    !any(locked_in & locked_out),
    ## solver
    assertthat::is.string(solver),
    assertthat::noNA(solver))
  assertthat::assert_that(
    solver %in% c("auto", "Rsymphony", "gurobi"))

  # identify solver
  if (identical(solver, "auto")) {
    solver <- ifelse(
      requireNamespace("gurobi", quietly = TRUE), "gurobi", "Rsymphony")
  }

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
  ## if using Rsymphony, then convert to Rsymphony model format
  if (identical(solver, "Rsymphony")) {
    m <- rsymphony_model(m)
  }
  ## create parameters list
  if (identical(solver, "gurobi")) {
    p <- list(Presolve = 2, MIPGap = 0, LogToConsole = as.integer(verbose))
  } else {
    p <- list(gap_limit = -1, verbosity = ifelse(verbose, 1, -2))
  }

  # prepare output matrix
  out <- matrix(NA, ncol = nrow(x), nrow = length(budget))
  # iterate over each budget and obtain the solution
  for (b in seq_along(budget)) {
    ## update budget
    m$rhs[1] <- budget[b]
    ## solve problem
    if (identical(solver, "gurobi")) {
      withr::with_locale(
        c(LC_CTYPE = "C"), {
        s <- gurobi::gurobi(m, p)$x
      })
    } else {
      s <- rsymphony_solve(m, p)$x
    }
    ## extract relevant values from solution
    s <- round(s[seq_len(n)])
    ## verify that solution is feasible
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
#'   is (`TRUE`) locked into the solution or (`FALSE`) not.
#'
#' @param locked_out `logical` vector indicating if each planning unit
#'   is (`TRUE`) locked out of the solution or (`FALSE`) not.
#'
#' @param verbose `logical` indicating if information should be
#'   printed while generating prioritization.
#'
#' @param solver `character` name of the optimization solver to use
#'   for generating survey schemes.
#'   Available options include: `"Rsymphony"`, `"gurobi"` and `"auto"`.
#'   The `"auto"` method will use the Gurobi optimization software if
#'   it is available; otherwise, it will use the SYMPHONY software
#'   via the \pkg{Rsymphony} package.
#'
#' @return A `matrix` containing the solution. Each column corresponds
#'   to a different planning unit and each row corresponds to a different
#'   row. Cell values indicate if a given planning unit was selected in a
#'   given solution.
#'
#' @noRd
weight_based_prioritizations <- function(
  x, budget, costs, locked_in, locked_out, solver, verbose) {
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
    !any(locked_in & locked_out),
    ## solver
    assertthat::is.string(solver),
    assertthat::noNA(solver))
  assertthat::assert_that(
    all(sum(locked_in * costs) <= budget),
    msg = "locked in sites exceed one or more budgets")
  assertthat::assert_that(
    solver %in% c("auto", "Rsymphony", "gurobi"))
  # identify solver
  if (identical(solver, "auto")) {
    solver <- ifelse(
      requireNamespace("gurobi", quietly = TRUE), "gurobi", "Rsymphony")
  }
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
  ## if using Rsymphony, then convert to Rsymphony model format
  if (identical(solver, "Rsymphony")) {
    m <- rsymphony_model(m)
  }
  ## create parameters list
  if (identical(solver, "gurobi")) {
    p <- list(Presolve = 2, MIPGap = 0, LogToConsole = as.integer(verbose))
  } else {
    p <- list(gap_limit = 0, verbosity = ifelse(verbose, 1, -2))
  }
  # prepare output matrix
  out <- matrix(NA, ncol = length(x), nrow = length(budget))
  # iterate over each budget and obtain the solution
  for (b in seq_along(budget)) {
    ## update budget
    m$rhs[1] <- budget[b]
    ## solve problem
    if (solver == "gurobi") {
      withr::with_locale(
        c(LC_CTYPE = "C"), {
        s <- gurobi::gurobi(m, p)$x
      })
    } else {
      s <- rsymphony_solve(m, p)$x
    }
    ## verify that solution is feasible
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

#' Rsymphony model
#'
#' This function converts an optimization problem
#' prepared for the [Gurobi](https://wwww.gurobi.com) optimization software
#' so that it can be used with the \pkg{Rsymphony} package.
#'
#' @param model `list` object.
#'
#' @details
#' Note that only a small set of parameters and model specification
#' methods are supported.
#'
#' @return A `list` object.
#'
#' @noRd
rsymphony_model <- function(model) {
  # assert arguments are valid
  assertthat::assert_that(is.list(model))
  # define recognized model components
  mn <- c("modelsense", "obj", "sense", "rhs", "lb", "ub", "vtype", "A")
  assertthat::assert_that(
    setequal(names(model), mn),
    msg = "argument to model contains unsupported components")

  # prepare model
  ## constraint matrix
  names(model)[names(model) == "A"] <- "mat"
  ## variables types
  names(model)[names(model) == "vtype"] <- "types"
  model$types[model$types == "S"] <- "C"
  ## model sense
  model$max <- model$modelsense == "max"
  model$modelsense <- NULL
  ## constraint sense
  names(model)[names(model) == "sense"] <- "dir"
  model$dir[model$dir == "="] <- "=="
  ## variable bounds
  model$bounds <- list(
    lower = list(ind = seq_along(model$lb), val = model$lb),
    upper = list(ind = seq_along(model$ub), val = model$ub))
    model$lb <- NULL
    model$ub <- NULL

  # return result
  model
}

#' Rsymphony solve
#'
#' This function is a wrapper to solve an optimization problem using the
#' the \pkg{Rsymphony} package.
#'
#' @param model `list` object.
#'
#' @param parameters `list` object.
#'
#' @return A `list` object.
#'
#' @noRd
rsymphony_solve <- function(model, parameters) {
  # assert arguments are valid
  assertthat::assert_that(is.list(model), is.list(parameters))
  # return solution
  s <- do.call(Rsymphony::Rsymphony_solve_LP, append(model, parameters))
  list(x = s$solution, objval = s$objval)
}
