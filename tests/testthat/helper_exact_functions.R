r_expected_value_of_management_action <- function(
  solution, prior_data, alpha, gamma) {
  # initialization
  sub_prior_data <- prior_data[, solution, drop = FALSE]
  sub_prior_data_log <- log(sub_prior_data)
  # iterate over each state
  states <- seq(1, rcpp_n_states(length(sub_prior_data)))
  out <- sapply(states, function(i) {
    s <- rcpp_nth_state(i, sub_prior_data)
    v <- sum((alpha * rowSums(s)) ^ gamma)
    p <- sum(s[] * sub_prior_data_log[]) +
         sum((1 - s[]) * (1 - sub_prior_data_log[]))
    v * exp(p)
  })
  sum(out)
}

r_expected_value_of_management_decision_given_current_information <- function(
  prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
  budget, gap) {
  # find optimal solution
  solution <- rcpp_prioritization(
    prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
    budget, gap, "")$x
  # calculate expected value
  r_expected_value_of_management_action(solution, prior_data, alpha, gamma)
}

r_expected_value_of_management_decision_given_perfect_information <- function(
  prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
  budget, gap) {
  # calculate log of prior data
  prior_data_log <- log(prior_data)
  # calculate expected value for each state
  states <- seq(1, rcpp_n_states(length(prior_data)))
  out <- sapply(states, function(i) {
    s <- rcpp_nth_state(i, prior_data)
    solution <- r_prioritization(
      s, pu_costs, as.numeric(pu_locked_in), alpha, gamma,
      n_approx_obj_fun_points, budget, gap, "")$x
    v <- sum((alpha * rowSums(s[, solution, drop = FALSE])) ^ gamma)
    if (v < 1e-10) return(NA_real_)
    p <- sum(s[] * prior_data_log[]) +
         sum((1 - s[]) * (1 - prior_data_log[]))
    log(v) + p
  })
  out <- out[is.finite(out)]
  exp(rcpp_log_sum(out))
}
