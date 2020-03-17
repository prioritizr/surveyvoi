r_posterior_probability_matrix <- function(
  rij, pij, oij,
  pu_survey_solution,
  survey_features,
  survey_sensitivity, survey_specificity,
  model_sensitivity, model_specificity) {
  # init
  n_pu <- ncol(rij)
  n_f <- nrow(rij)
  n_f_survey <- sum(survey_features)
  out <- matrix(NA, nrow = n_f, ncol = n_pu)
  survey_features_rev_idx <- calculate_survey_features_rev_idx(survey_features)

  # prepare feature data
  feature_outcome_idx <- rep(1, n_f_survey)
  model_sensitivity2 <- matrix(model_sensitivity[survey_features], ncol = 1)
  model_specificity2 <- matrix(model_specificity[survey_features], ncol = 1)
  pij_survey_species_subset <- pij[survey_features, , drop = FALSE ]

  # calculate total survey probabilities
  total_probability_of_survey_positive <-
    total_probability_of_positive_result(pij, survey_sensitivity,
                                         survey_specificity)
  total_probability_of_survey_negative <-
    total_probability_of_negative_result(pij, survey_sensitivity,
                                         survey_specificity)

  # calculate total model probabilities
  total_probability_of_model_positive <-
    total_probability_of_positive_model_result(
      pij_survey_species_subset, model_sensitivity2, model_specificity2,
      feature_outcome_idx)
  total_probability_of_model_negative <-
    total_probability_of_negative_model_result(
      pij_survey_species_subset, model_sensitivity2, model_specificity2,
      feature_outcome_idx)

  # main calculations
  for (j in seq_len(n_pu)) {
    for (i in seq_len(n_f)) {
      if (rij[i, j] > -0.5) {
        out[i, j] <- pij[i, j]
      } else if (!survey_features[i]) {
        out[i, j] <- pij[i, j]
      } else if ((pu_survey_solution[j] > 0.5) && (oij[i, j] >= 0.5)) {
        out[i, j] <-
          (survey_sensitivity[i] * pij[i, j]) /
           total_probability_of_survey_positive[i, j]
      } else if ((pu_survey_solution[j] > 0.5) && (oij[i, j] < 0.5)) {
        out[i, j] <-
          ((1 - survey_sensitivity[i]) * pij[i, j]) /
           total_probability_of_survey_negative[i, j]
      } else {
        sub_i <- survey_features_rev_idx[i]
        if (oij[i, j] >= 0.5) {
          out[i, j] <-
            (model_sensitivity2[sub_i, 1] * pij[i, j]) /
             total_probability_of_model_positive[sub_i, j]
        } else {
          out[i, j] <-
            ((1.0 - model_sensitivity2[sub_i, 1]) * pij[i, j]) /
             total_probability_of_model_negative[sub_i, j]
        }
      }
    }
  }
  # return result
  out
}
