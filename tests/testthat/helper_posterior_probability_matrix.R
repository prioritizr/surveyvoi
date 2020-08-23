r_posterior_probability_matrix <- function(
  nij, pij, oij,
  pu_survey_solution,
  survey_features,
  survey_sensitivity, survey_specificity,
  model_sensitivity, model_specificity) {
  # init
  n_pu <- ncol(pij)
  n_f <- nrow(pij)
  n_f_survey <- sum(survey_features)
  out <- matrix(NA, nrow = n_f, ncol = n_pu)
  survey_features_rev_idx <- calculate_survey_features_rev_idx(survey_features)

  # prepare feature data
  model_sensitivity2 <- model_sensitivity[survey_features]
  model_specificity2 <- model_specificity[survey_features]
  pij_survey_species_subset <- pij[survey_features, , drop = FALSE ]

  # calculate total survey probabilities
  total_probability_of_survey_positive <-
    total_probability_of_positive_result(
      pij, survey_sensitivity, survey_specificity)
  total_probability_of_survey_negative <-
    total_probability_of_negative_result(
      pij, survey_sensitivity, survey_specificity)

  # calculate total model probabilities
  total_probability_of_model_positive <-
    total_probability_of_positive_model_result(
      pij_survey_species_subset, model_sensitivity2, model_specificity2)
  total_probability_of_model_negative <-
    total_probability_of_negative_model_result(
      pij_survey_species_subset, model_sensitivity2, model_specificity2)

  # main calculations
  out <- matrix(NA, ncol = n_pu, nrow = n_f)
  for (j in seq_len(n_pu)) {
    for (i in seq_len(n_f)) {
      if (!survey_features[i] || (nij[i, j] > 0)) {
        # if the species is not being surveyed, then use prior data
        # or, if the planning unit already has survey data then use prior data
        out[i, j] <- pij[i, j]
      } else if ((pu_survey_solution[j] > 0.5) && (oij[i, j] >= 0.5)) {
        # if there is survey data for i'th pu and j'th species,
        # and the species was detected
        out[i, j] <-
          (survey_sensitivity[i] * pij[i, j]) /
           total_probability_of_survey_positive[i, j]
      } else if ((pu_survey_solution[j] > 0.5) && (oij[i, j] < 0.5)) {
        # if there is survey data for i'th pu and j'th species,
        # and the species was not detected
        out[i, j] <-
          ((1 - survey_sensitivity[i]) * pij[i, j]) /
           total_probability_of_survey_negative[i, j]
      } else {
        # if there is no survey data for i'th pu and j'th species,
        # then use model predictions
        sub_i <- survey_features_rev_idx[i]
        if (oij[i, j] >= 0.5) {
          # if the model predicts a presence
          out[i, j] <-
            (model_sensitivity2[sub_i] * pij[i, j]) /
             total_probability_of_model_positive[sub_i, j]
        } else {
          # if the model predicts an absence
          out[i, j] <-
            ((1.0 - model_sensitivity2[sub_i]) * pij[i, j]) /
             total_probability_of_model_negative[sub_i, j]
        }
      }
    }
  }
  # return result
  out
}
