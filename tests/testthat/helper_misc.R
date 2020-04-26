calculate_survey_features_rev_idx <- function(x) {
  n_f <- length(x)
  out <- rep(0, n_f);
  k <- 0
  for (i in seq_len(n_f)) {
    if (x[i]) {
      out[i] <- k;
      k <- k + 1
    }
  }
  out + 1
}

log_sum <- function(u, v) {
  max(u, v) + log(exp(u - max(u, v)) + exp(v -max(u, v)))
}

se <- function(x) sqrt(var(x) / length(x))

r_conservation_value_state <- function(
  state, preweight, postweight, target, total) {
  v <- rowSums(state)
  sum(ifelse(
    v < target,
    preweight * (v / target),
    preweight + (postweight * ((v - target) / (total - target)))))
}

r_conservation_value_amount <- function(amount, preweight, postweight,
  target, total) {
  assertthat::assert_that(is.numeric(amount),
  assertthat::is.number(preweight), assertthat::is.number(postweight),
  assertthat::is.number(target))
  ifelse(amount < target,
    preweight * (amount / target),
    preweight + (postweight * ((amount - target) / (total - target))))
}

wrap <- function(x) return(x)
