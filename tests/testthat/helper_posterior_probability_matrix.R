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
      if (!survey_features[i]) {
        # if the species is not being surveyed, then use prior data
        out[i, j] <- pij[i, j]
      } else {
        # calculate indices
        sub_i <- survey_features_rev_idx[i]
        has_no_survey_data <- (nij[i, j] < 0.5) && !pu_survey_solution[j]
        # calculate model TSS
        curr_model_tss <-
          model_sensitivity2[sub_i] + model_specificity2[sub_i] - 1
        # calculate survey methodology TSS (account for multiple surveys)
        curr_survey_sens <-
          1 - (prod(1 - rep(survey_sensitivity[i], nij[i, j])))
        curr_survey_spec <-
          1 - (prod(1 - rep(survey_specificity[i], nij[i, j])))
        curr_survey_tss <- curr_survey_sens + curr_survey_spec - 1
        is_model_better_than_survey <- curr_model_tss > curr_survey_tss
        # generate posterior
        if (has_no_survey_data || is_model_better_than_survey) {
          ## use model for posterior if there is no survey data or
          ## the models are better
          if (oij[i, j] >= 0.5) {
            ## if model predicts presence
            out[i, j] <-
              (model_sensitivity2[sub_i] * pij[i, j]) /
               total_probability_of_model_positive[sub_i, j]
          } else {
            ## if model predicts absence
            out[i, j] <-
              ((1.0 - model_sensitivity2[sub_i]) * pij[i, j]) /
               total_probability_of_model_negative[sub_i, j]
          }
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
  }
  # return result
  out
}

r_find_rij_idx_based_on_models <- function(
  nij,
  pu_survey_solution,
  survey_features, survey_features_rev_idx,
  survey_sensitivity, survey_specificity,
  model_sensitivity, model_specificity) {
  out <- c()
  counter <- 1 # base 1 indexing
  for (j in seq_len(ncol(nij))) {
    for (i in seq_len(nrow(nij))) {
      if ((nij[i, j] < 0.5) &&
          (survey_features[i] > 0.5) &&
          (pu_survey_solution[j] < 0.5)) {
        out <- c(out, counter)
      }
      counter <- counter + 1
    }
  }
  out
}
