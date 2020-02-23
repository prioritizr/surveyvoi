r_approx_expected_value_of_action <- function(
  solution, prior_data, alpha, gamma, states) {
  # initialization
  sub_prior_data <- prior_data[, solution, drop = FALSE]
  sub_prior_data_log <- log(sub_prior_data)
  # iterate over each state
  out <- sapply(states, function(i) {
    s <- rcpp_nth_state(i, sub_prior_data)
    v <- sum((alpha * rowSums(s)) ^ gamma)
    p <- sum(s[] * sub_prior_data_log[]) +
         sum((1 - s[]) * (1 - sub_prior_data_log[]))
    c(v, p)
  })
  out[1, ] <- log(out[1, ])
  out <- out[, is.finite(out[1, ]), drop = FALSE]
  out[2, ] <- out[2, ] / sum(out[2, ])
  exp(rcpp_log_sum(colSums(out)))
}

r_approx_expected_value_of_decision_given_current_info_fixed_states <- function(
  prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
  budget, gap, states) {
  # find optimal solution
  solution <- rcpp_prioritization(
    prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
    budget, gap, "")$x
  # calculate expected value
  r_approx_expected_value_of_action(solution, prior_data,
    alpha, gamma, states)
}

r_approx_expected_value_of_decision_given_current_info_n_states <- function(
  prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
  budget, gap, n) {
  # find optimal solution
  solution <- rcpp_prioritization(
    prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
    budget, gap, "")$x
  # calculate expected value
  r_approx_expected_value_of_action(solution, prior_data,
    alpha, gamma, rcpp_sample_k_nth_states(n, prior_data))
}
