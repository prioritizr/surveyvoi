r_prioritization <- function(rij, pu_costs, pu_locked_in, preweight,
  postweight, target, n_approx_points, budget, gap, file_path) {
  # assert that arguments are valid
  assertthat::assert_that(
    is.matrix(rij), ncol(rij) > 0, nrow(rij) > 0, assertthat::noNA(c(rij)),
    is.numeric(preweight), is.numeric(postweight), is.numeric(target),
    length(preweight) == nrow(rij), length(postweight) == nrow(rij),
    length(target) == nrow(rij),
    is.numeric(pu_costs), length(pu_costs) == ncol(rij),
    assertthat::noNA(pu_costs),
    is.numeric(pu_locked_in), length(pu_locked_in) == ncol(rij),
    assertthat::noNA(pu_locked_in),
    all(pu_locked_in %in% c(0, 1)),
    assertthat::is.count(n_approx_points), assertthat::noNA(n_approx_points),
    assertthat::is.number(budget), assertthat::noNA(budget),
    assertthat::is.number(gap), assertthat::noNA(gap))
  # init
  n_pu <- ncol(rij)
  n_f <- nrow(rij)
  # build pwl components
  n_approx_points <- n_approx_points - 1
  obj_feature_held <- apply(rij, 1, function(x) {
    r1 <- min(x) * 0.99
    r2 <- sum(x) * 1.01
    seq(r1, r2, length.out = n_approx_points)
  })
  obj_feature_benefit <- obj_feature_held
  {for (i in seq_len(ncol(obj_feature_benefit))) {
    obj_feature_benefit[, i] <-
      r_conservation_benefit_amount(obj_feature_held[, i],
        preweight[i], postweight[i], target[i])
  }}
  # build problem
  p <- list()
  p$modelsense <- "max"
  p$ub <- c(rep(1, n_pu), rep(Inf, n_f))
  p$lb <- c(pu_locked_in, rep(0, n_f))
  p$obj <- rep(0, n_pu + n_f)
  p$rhs <- c(rep(0, n_f), budget)
  p$vtype <- c(rep("B", n_pu), rep("C", n_f))
  p$sense <- c(rep("=", n_f), "<=")
  p$A <- rbind(cbind(rij, diag(n_f) * -1), c(pu_costs, rep(0, n_f)))
  p$pwlobj <- lapply(seq_len(n_f), function(i) {
    if (abs(diff(range(obj_feature_benefit[, i]))) > 1e-5) {
      o <- list(var = n_pu + i,
                x = c(0, obj_feature_held[, i]),
                y = c(0, obj_feature_benefit[, i]))
    } else if (obj_feature_benefit[1, i] > 1e-5) {
      o <- list(var = n_pu + i,
                x = c(0, obj_feature_held[1, i]),
                y = c(0, obj_feature_benefit[1, i]))
    } else {
      o <- list(var = n_pu + i, x = c(0, 0.5, 1), y = c(0, 0.5, 1))
    }
    o
  })
  # solve problem
  g <- gurobi::gurobi(p, list(MIPGap = gap, Presolve = -1, OutputFlag = 0))
  assertthat::assert_that(identical(g$status, "OPTIMAL"))
  # write problem to disk if needed
  if (nchar(file_path) > 0)
    gurobi::gurobi_write(p, file_path)
  # return result
  list(x = as.logical(g$x[seq_len(n_pu)]), objval = g$objval, full_x = g$x)
}

brute_force_prioritization <- function(rij, preweight, postweight, target,
                                       pu_costs, pu_locked_in, budget) {
  # assert that arguments are valid
  assertthat::assert_that(
    is.matrix(rij), ncol(rij) > 0, nrow(rij) > 0, assertthat::noNA(c(rij)),
    is.numeric(preweight), is.numeric(postweight), is.numeric(target),
    is.numeric(pu_costs), length(pu_costs) == ncol(rij),
    assertthat::noNA(pu_costs),
    is.numeric(pu_locked_in), length(pu_locked_in) == ncol(rij),
    assertthat::noNA(pu_locked_in),
    all(pu_locked_in %in% c(0, 1)),
    assertthat::is.number(budget), assertthat::noNA(budget))
  # set constants
  n_pu <- ncol(rij)
  n_f <- nrow(rij)
  m <- matrix(0, ncol = n_pu, nrow = 1)
  # brute force all combinations
  objs <- vapply(seq_len(rcpp_n_states(ncol(rij))), FUN.VALUE = numeric(1),
                 function(i) {
    # generate nth solution
    x <- rcpp_nth_state(i, m)
    # if any locked in planning units are not selected return -Inf
    if (any(x < pu_locked_in))
      return(-Inf)
    # if budget is exceeded then return -Inf
    if (sum(x * pu_costs) > budget)
      return(-Inf)
    # if feasible solution then return objective function
    r_conservation_benefit_state(
      x[rep(1, n_f), ] * rij, preweight, postweight, target)
  })
  # find best solution
  list(x = as.logical(c(rcpp_nth_state(which.max(objs), m))),
       objval = max(objs))
}
