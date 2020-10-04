r_update_posterior_probability_matrix <- function(
  pij, oij,
  survey_features,
  survey_sensitivity, survey_specificity,
  pu_survey_solution) {
  # init
  n_pu <- ncol(pij)
  n_f <- nrow(pij)
  n_f_survey <- sum(survey_features)
  out <- matrix(NA, nrow = n_f, ncol = n_pu)

  # calculate total survey probabilities
  total_probability_of_survey_positive <-
    total_probability_of_positive_result(
      pij, survey_sensitivity, survey_specificity)
  total_probability_of_survey_negative <-
    total_probability_of_negative_result(
      pij, survey_sensitivity, survey_specificity)

  # main calculations
  out <- matrix(NA, ncol = n_pu, nrow = n_f)
  for (j in seq_len(n_pu)) {
    for (i in seq_len(n_f)) {
      if (!(survey_features[i] && pu_survey_solution[j])) {
        # if the species/planning unit is not being surveyed,
        # then use prior data
        out[i, j] <- pij[i, j]
      } else {
        ## use survey data
        if (oij[i, j] >= 0.5) {
          ## if survey gives detection
          out[i, j] <-
            (survey_sensitivity[i] * pij[i, j]) /
             total_probability_of_survey_positive[i, j]
        } else {
          ## if survey gives non-detection
          out[i, j] <-
            ((1 - survey_sensitivity[i]) * pij[i, j]) /
             total_probability_of_survey_negative[i, j]
        }
      }
    }
  }
  out[] <- pmax(out[], 1e-10)
  out[] <- pmin(out[], 1 - 1e-10)
  # return result
  out
}
