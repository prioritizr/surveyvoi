r_approximate_expected_value_of_management_action <- function(
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
