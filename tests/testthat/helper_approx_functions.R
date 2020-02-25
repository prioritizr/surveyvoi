r_approx_expected_value_of_action <- function(
  solution, prior_data, alpha, gamma, states) {
  # initialization
  sub_prior_data <- prior_data[, solution, drop = FALSE]
  sub_prior_data_log <- log(sub_prior_data)
  sub_prior_data_log1m <- log(1 - sub_prior_data)
  # iterate over each state
  out <- sapply(states, function(i) {
    s <- rcpp_nth_state(i, sub_prior_data)
    v <- sum((alpha * rowSums(s)) ^ gamma)
    p <- sum(s[] * sub_prior_data_log[]) +
         sum((1 - s[]) * sub_prior_data_log1m[])
    c(v, p)
  })
  out[1, ] <- log(out[1, ])
  out <- out[, is.finite(out[1, ]), drop = FALSE]
  out[2, ] <- out[2, ] - rcpp_log_sum(out[2, ])
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
  budget, gap, n_replicates, n_states_per_replicate) {
  # find optimal solution
  solution <- rcpp_prioritization(
    prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
    budget, gap, "")$x
  # calculate expected value
  value <- sapply(seq_len(n_replicates), function(i) {
    r_approx_expected_value_of_action(
      solution, prior_data, alpha, gamma,
      rcpp_sample_k_uniform_no_replacement_nth_states(
        n_states_per_replicate, prior_data))
  })
  c(mean(value), se(value))
}

r_approx_expected_value_of_decision_given_perfect_info_fixed_states <- function(
  prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
  budget, gap, states) {
  # calculate log of prior data
  prior_data_log <- log(prior_data)
  prior_data_log1m <- log(1 - prior_data)
  # calculate expected value for each state
  out <- sapply(states, function(i) {
    s <- rcpp_nth_state(i, prior_data)
    solution <- r_prioritization(
      s, pu_costs, as.numeric(pu_locked_in), alpha, gamma,
      n_approx_obj_fun_points, budget, gap, "")$x
    v <- sum((alpha * rowSums(s[, solution, drop = FALSE])) ^ gamma)
    p <- sum(s[] * prior_data_log[]) +
         sum((1 - s[]) * prior_data_log1m[])
    c(v, p)
  })
  out[1, ] <- log(out[1, ])
  out <- out[, is.finite(out[1, ]), drop = FALSE]
  out[2, ] <- out[2, ] - rcpp_log_sum(out[2, ])
  exp(rcpp_log_sum(colSums(out)))
}

r_approx_expected_value_of_decision_given_perfect_info_n_states <- function(
  prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
  budget, gap, n_replicates, n_states_per_replicate) {
  value <- sapply(seq_len(n_replicates), function(i) {
    r_approx_expected_value_of_decision_given_perfect_info_fixed_states(
      prior_data, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points,
      budget, gap,
      rcpp_sample_k_uniform_no_replacement_nth_states(
        n_states_per_replicate, prior_data))
  })
  c(mean(value), se(value))
}

se <- function(x) sqrt(var(x) / length(x))
