r_prioritization <- function(
  rij, pu_costs, pu_locked_in, pu_locked_out, target, budget, gap = 0,
  file_path = "") {
  # assert that arguments are valid
  assertthat::assert_that(
    is.matrix(rij), ncol(rij) > 0, nrow(rij) > 0, assertthat::noNA(c(rij)),
    is.numeric(target), length(target) == nrow(rij), assertthat::noNA(target),
    length(unique(target)) == 1,
    is.numeric(pu_costs), length(pu_costs) == ncol(rij),
    ncol(rij) >= max(target),
    assertthat::noNA(pu_costs),
    is.numeric(pu_locked_in), length(pu_locked_in) == ncol(rij),
    assertthat::noNA(pu_locked_in),
    all(pu_locked_in %in% c(0, 1)),
    is.numeric(pu_locked_out), length(pu_locked_out) == ncol(rij),
    assertthat::noNA(pu_locked_out),
    all(pu_locked_out %in% c(0, 1)),
    assertthat::is.number(budget), assertthat::noNA(budget),
    assertthat::is.number(gap), assertthat::noNA(gap))
  assertthat::assert_that(
    sum(pu_costs[pu_locked_in > 0.5]) <= budget,
    msg = paste("locked in planning units exceed budget"))
  assertthat::assert_that(
    sum(sort(pu_costs[pu_locked_out < 0.5])[seq_len(max(target))]) <= budget,
    msg = paste("budget too low to obtain solution with a number of planning
                units >= max(targets)"))
 # init
  n_pu <- ncol(rij)
  n_f <- nrow(rij)
  # transform rij data
  pij <- rij
  pij[pij > (1 - 1e-10)] <- (1 - 1e-10)
  pij <- log(1 - pij) * 1e3
  # build problem
  p <- list()
  p$modelsense <- "min"
  p$ub <- c(1 - pu_locked_out, rep(0, n_f + 1))
  p$lb <- c(pu_locked_in, rep(-Inf, n_f + 1))
  p$obj <- c(rep(0, n_pu + n_f), 1)
  p$vtype <- c(rep("B", n_pu), rep("C", n_f + 1))
  p$rhs <- c(rep(0, n_f), rep(0, n_f), target[1], budget)
  p$sense <- c(rep("=", n_f), rep("<=", n_f), ">=", "<=")
  p$A <- rbind(
    cbind(pij, diag(n_f) * -1, 0),
    cbind(matrix(0, ncol = n_pu, nrow = n_f), diag(n_f), -1),
    c(rep(1, n_pu), rep(0, n_f + 1)),
    c(pu_costs, rep(0, n_f + 1)))
  # solve problem
  g <- gurobi::gurobi(p, list(MIPGap = gap, Presolve = -1, OutputFlag = 0,
                              Seed = 1, NumericFocus = 0, Threads = 1))
  assertthat::assert_that(identical(g$status, "OPTIMAL"))
  # write problem to disk if needed
  if (nchar(file_path) > 0)
    gurobi::gurobi_write(p, file_path)
  # return result
  list(x = as.logical(g$x[seq_len(n_pu)]), objval = g$objval, full_x = g$x)
}

brute_force_prioritization <- function(
  rij, pu_costs, pu_locked_in, pu_locked_out, target, budget) {
  # assert that arguments are valid
  assertthat::assert_that(
    is.matrix(rij), ncol(rij) > 0, nrow(rij) > 0, assertthat::noNA(c(rij)),
    is.numeric(target), length(target) == nrow(rij), assertthat::noNA(target),
    ncol(rij) >= max(target),
    is.numeric(pu_costs), length(pu_costs) == ncol(rij),
    assertthat::noNA(pu_costs),
    is.numeric(pu_locked_in), length(pu_locked_in) == ncol(rij),
    assertthat::noNA(pu_locked_in),
    all(pu_locked_in %in% c(0, 1)),
    is.numeric(pu_locked_out), length(pu_locked_out) == ncol(rij),
    assertthat::noNA(pu_locked_out),
    all(pu_locked_out %in% c(0, 1)),
    assertthat::is.number(budget), assertthat::noNA(budget))
  # set constants
  n_pu <- ncol(rij)
  n_f <- nrow(rij)
  m <- matrix(0, ncol = n_pu, nrow = 1)
  # brute force all combinations
  objs <- vapply(seq_len(rcpp_n_states(ncol(rij))), FUN.VALUE = numeric(1),
                 function(i) {
    # generate nth solution
    x <- c(rcpp_nth_state(i, m))
    # if number of selected planning units less than target return -Inf
    if (sum(x) < max(target))
      return(-Inf)
    # if any locked in planning units are not selected return -Inf
    if (any(x < pu_locked_in))
      return(-Inf)
    # if any locked out planning units are selected return -Inf
    if (any((x + pu_locked_out) == 2))
      return(-Inf)
    # if budget is exceeded then return -Inf
    if (sum(x * pu_costs) > budget)
      return(-Inf)
    # if feasible solution then return objective function
    rcpp_expected_value_of_action(as.logical(x), rij, target)
  })
  # find best solution
  list(x = as.logical(c(rcpp_nth_state(which.max(objs), m))),
       objval = max(objs))
}
