# https://stats.stackexchange.com/a/41263/103849
convolve_binomial <- function(p, threshold) {
  # p is a vector of probabilities of Bernoulli distributions.
  # The convolution of these distributions is returned as a vector
  # `z` where z[i] is the probability of i-1, i=1, 2, ..., length(p)+1.
  n <- length(p) + 1
  z <- c(1, rep(0, n - 1))
  for (p in p) {
    z <- (1 - p) * z + p * c(0, z[-n])
  }
  1 - sum(z[seq_len(round(threshold))])
}

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

# https://stackoverflow.com/a/778273/3483791
log_minus <- function(u, v) {
  u + log(1 + (-exp(v - u)))
}

se <- function(x) sqrt(var(x) / length(x))

r_conservation_value <- function(pij, target) {
  assertthat::assert_that(all(target <= ncol(pij)))
  out <- vapply(seq_along(target), FUN.VALUE = numeric(1), function(i) {
    convolve_binomial(pij[i, ], target[i])
  })
  prod(out)
}

r_approx_conservation_value <- function(pij, target) {
  assertthat::assert_that(all(target <= ncol(pij)))
  out <- vapply(seq_along(target), FUN.VALUE = numeric(1), function(i) {
    sum(PoissonBinomial::dpbinom(
      seq(target[i], ncol(pij)), pij[i, ], method = "Normal"))
  })
  prod(out)
}

r_proxy_conservation_value <- function(pij, target) {
  r_conservation_value(pij, rep(1, nrow(pij)))
}

r_conservation_value_amount <- function(amount, target) {
  assertthat::assert_that(is.numeric(amount), is.numeric(target))
  out <- rep(0, length(amount))
  not_met_idx <- which(amount < target)
  out[not_met_idx] <- 1 - (amount[not_met_idx] / target[not_met_idx])
  sum(out)
}

wrap <- function(x) return(x)
