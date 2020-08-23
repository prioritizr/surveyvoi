#' Space based prioritizations
#'
#' This function generates space-based prioritizations by solving the
#' \emph{p}-median problem to optimality.
#'
#' @param x \code{matrix} object containing a symmetric distance matrix that
#'   measures the dissimilarity between different planning units.
#'
#' @param budget \code{numeric} vector indicating the maximum expenditure for
#'   each solution. The values must be between one and
#'   the number of planning units (i.e. length of \code{x}).
#'
#' @param costs \code{numeric} cost values for each planning units. Defaults
#'   to a vector assigning a cost of one for each planning unit.
#'
#' @param locked_in \code{logical} vector indicating if each planning unit
#'   is (\code{TRUE}) locked into the solution or (\code{FALSE}) not. Defaults
#'   to a vector of \code{FALSE} for each planning unit.
#'
#' @param locked_out \code{logical} vector indicating if each planning unit
#'   is (\code{TRUE}) locked out of the solution or (\code{FALSE}) not. Defaults
#'   to a vector of \code{FALSE} for each planning unit.
#'
#' @param verbose \code{logical} indicating if information should be
#'   printed while generating prioritization. Defaults to \code{FALSE}.
#'
#' @return \code{matrix} containing the solution. Each column corresponds
#'   to a different planning unit and each row corresponds to a different
#'   solution. Cell values indicate if a given planning unit was selected in a
#'   given solution.
#'
#' @seealso For more information on the \emph{p}-median problem
#'   see \url{http://www.hyuan.com/java/how.html}.
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
    s <- gurobi::gurobi(m, p)$x[seq_len(n)]
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
#' @param x \code{numeric} vector containing the weights.
#'
#' @param budget \code{numeric} number indicating how many planning units
#'   should be selected in the solution. The value must be between one and
#'   the number of planning units (i.e. length of \code{x}).
#'
#' @param locked_in \code{logical} vector indicating if each planning unit
#'   is (\code{TRUE}) locked into the solution or (\code{FALSE}) not. Defaults
#'   to a vector of \code{FALSE} for each planning unit.
#'
#' @param locked_out \code{logical} vector indicating if each planning unit
#'   is (\code{TRUE}) locked out of the solution or (\code{FALSE}) not. Defaults
#'   to a vector of \code{FALSE} for each planning unit.
#'
#' @param verbose \code{logical} indicating if information should be
#'   printed while generating prioritization. Defaults to \code{FALSE}.
#'
#' @return \code{matrix} containing the solution. Each column corresponds
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
    s <- gurobi::gurobi(m, p)$x
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

#' Create K-folds using site data.
#'
#' Create k-folds given survey data in multiple sites.
#'
#' @param prop_detected \code{integer} proportion of surveys for each site
#'   within which the species was detected. Each
#'   element corresponds to a different site, and values indicate the
#'   proportion of times a species was detected within a given site.
#'   If a site does not have any detections, then a value of zero should be
#'   used (not \code{NA}).
#'
#' @param n_total \code{integer} number of total surveys conducted within
#'   each site.
#'   Each element corresponds to a different site, and values indicate the
#'   number of surveys conducted within the given site.
#'   If a site does not have any non-detections, then a value of zero should be
#'   used (not \code{NA}).
#'
#' @param n \code{numeric} number of folds.
#'
#' @param index \code{integer} indices associated with each site.
#'   Defaults to a sequence ranging from 1 to the cardinality of the
#'   argument to \code{x} (i.e. \code{seq_along(x)}).
#'
#' @param seed \code{numeric} random number generated seed for generating
#'   folds. Defaults to 500.
#'
#' @details
#'  The sites will be stratified into folds will be stratified to ensure that
#'  each fold contains least one detection and one non-detection in the
#'  training and test datasets for subsequent model fitting. Note that
#'  sites with have zero detections and zero non-detections are
#'  randomly allocated to folds.
#'
#' @return \code{list} of \code{list} objects containing the
#'  indices excluded from each fold.
#'
#' @noRd
create_site_folds <- function(
  prop_detected, n_total, n, index = seq_along(x), seed = 500) {
  # assert arguments are valid
  assertthat::assert_that(
    is.numeric(prop_detected), length(prop_detected) > 0,
    assertthat::noNA(prop_detected),
    all(prop_detected >= 0), all(prop_detected <= 1),
    is.numeric(n_total), length(n_total) > 0,
    all(n_total >= 0), assertthat::noNA(n_total),
    identical(length(n_total), length(n_total)),
    assertthat::is.count(n),
    assertthat::noNA(n),
    is.numeric(index), assertthat::noNA(index),
    identical(length(prop_detected), length(index)),
    assertthat::is.count(seed))
  assertthat::assert_that(
    max(abs(round(n_total) - n_total)) < 1e-10,
    msg = "argument to n_total does not contain whole numbers")
  assertthat::assert_that(sum(round(prop_detected * n_total) > 0) >= n,
    msg = "not enough presences to create the specified number of folds")
  assertthat::assert_that(sum(round((1 - prop_detected) * n_total) > 0) >= n,
    msg = "not enough absence to create the specified number of folds")

  # initialization
  n_det <- round(prop_detected * n_total)
  n_nondet <- round((1 - prop_detected) * n_total)
  site_data <- tibble::tibble(
    idx = index, n_det = n_det, n_nondet = n_nondet, n_total = n_total)
  obs_y <- c(rep(rep(1, length(n_det)), n_det),
             rep(rep(0, length(n_nondet)), n_nondet))
  obs_index <- c(rep(index, n_det), rep(index, n_nondet))
  obs_data <- tibble::tibble(y = obs_y, idx = obs_index,
                             idf = factor(as.character(obs_index)))

  # organize site data with observations into folds
  withr::with_seed(seed, {
    # create folds
    obs_data2 <- groupdata2::fold(
      obs_data, num_col = "y", id_col = "idf", k = n, num_fold_cols = 5)
  })

  # find valid fold
  fold_columns <- setdiff(names(obs_data2), names(obs_data))
  found_valid <- FALSE
  for (f in fold_columns) {
    ## format fold data
    obs_data2$fold <- as.integer(as.character(obs_data2[[f]]))
    ## calculate statistics to determine if folding scheme is valid
    n_det_per_fold <- aggregate(obs_data2$y, by = list(obs_data2$fold), sum)$x
    n_nondet_per_fold <-
      aggregate(1 - obs_data2$y, by = list(obs_data2$fold), sum)$x
    ## if folding scheme is valid, then keep it
    if (all(n_det_per_fold > 0) && all(n_nondet_per_fold > 0)) {
      found_valid <- TRUE
      break()
    }
  }

  # throw error if no valid folding scheme was found, then throw error
  assertthat::assert_that(found_valid,
    msg = paste("could not find any valid folding schemes that have at least",
                "one detection and non-detection per fold, try again with a",
                "different seed."))

  # determine which fold each site belongs to
  site_data$fold <- vapply(site_data$idx, FUN.VALUE = integer(1), function(x) {
    if (x %in% obs_data2$idx) {
      out <- as.integer(obs_data2$fold[obs_data2$idx == x][[1]])
    } else {
      out <- NA_integer_
    }
  })

  # randomly allocate any sites that are missing fold values
  # (because they have no previous detections or non-detections)
  na_pos <- is.na(site_data$fold)
  if (any(na_pos)) {
    site_data$fold[na_pos] <-
      sample(seq_len(n), sum(na_pos), replace = sum(na_pos) > n)
  }

  # extract indices for folds
  site_data2 <- split(site_data, site_data$fold)
  train <- lapply(site_data2, function(i) setdiff(index, i$idx))
  test <- lapply(site_data2, function(i) i$idx)

  # return result
  list(train = train, test = test)
}
