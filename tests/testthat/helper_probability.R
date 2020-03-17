total_probability_of_positive_result <- function(
  pij, survey_sensitivity, survey_specificity) {
  n_f <- nrow(pij)
  n_pu <- ncol(pij)
  out <- matrix(NA, nrow = n_f, ncol = n_pu)
  for (i in seq_len(n_pu)) {
    out[, i] <-
     (survey_sensitivity * pij[, i]) +
     ((1.0 - survey_specificity) * (1.0 - pij[, i]))
  }
  out
}

total_probability_of_negative_result <- function(
  pij, survey_sensitivity, survey_specificity) {
  n_f <- nrow(pij)
  n_pu <- ncol(pij)
  out <- matrix(NA, nrow = n_f, ncol = n_pu)
  for (i in seq_len(n_pu)) {
    out[, i] <-
     ((1 - survey_sensitivity) * pij[, i]) +
     (survey_specificity * (1.0 - pij[, i]))
  }
  out
}

total_probability_of_positive_model_result <- function(
  pij_subset, model_sensitivity, model_specificity, feature_outcome_idx) {
  n_f <- length(feature_outcome_idx)
  n_pu <- ncol(pij_subset)
  out <- matrix(NA, nrow = n_f, ncol = n_pu)
  for (j in seq_len(n_pu)) {
    for (i in seq_along(feature_outcome_idx)) {
    out[i, j] <-
     (model_sensitivity[i, feature_outcome_idx[i]] *
      pij_subset[i, j]) +
     ((1.0 - model_specificity[i, feature_outcome_idx[i]]) *
      (1.0 - pij_subset[i, j]))
  }}
  out
}

total_probability_of_negative_model_result <- function(
  pij_subset, model_sensitivity, model_specificity, feature_outcome_idx) {
  n_f <- length(feature_outcome_idx)
  n_pu <- ncol(pij_subset)
  out <- matrix(NA, nrow = n_f, ncol = n_pu)
  for (j in seq_len(n_pu)) {
    for (i in seq_along(feature_outcome_idx)) {
    out[i, j] <-
     ((1 - model_sensitivity[i, feature_outcome_idx[i]]) *
      pij_subset[i, j]) +
     ((model_specificity[i, feature_outcome_idx[i]]) *
      (1.0 - pij_subset[i, j]))
  }}
  out
}
