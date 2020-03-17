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
