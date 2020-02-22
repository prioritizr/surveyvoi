context("posterior_probability_matrix")

# R implementation
r_posterior_probability_matrix <- function(
  rij, pij, oij, pu_survey_solution,
  survey_sensitivity, survey_specificity,
  survey_features, survey_features_rev_idx,
  model_sensitivity,
  total_probability_of_survey_positive, total_probability_of_survey_negative,
  total_probability_of_model_positive, o, out) {
  # init
  n_pu <- ncol(rij)
  n_f <- nrow(rij)
  survey_features_rev_idx <- survey_features_rev_idx + 1
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
          1 - ((survey_sensitivity[i] * pij[i, j]) /
               total_probability_of_survey_positive[i, j])
      } else {
        sub_i <- survey_features_rev_idx[i]
          out[i, j] <-
             (oij[i, j] * model_sensitivity[sub_i, o[sub_i]] * pij[i, j]) /
              total_probability_of_model_positive[sub_i, j]
      }
    }
  }
  # return result
  out
}

test_that("correct result", {
  # data
  ## constants
  set.seed(123)
  n_pu <- 10
  n_f <- 8
  n_f_survey <- 4
  n_o <- 5
  prop_missing_data <- 0.5
  pu_not_surveyed <- sort(sample.int(n_pu, ceiling(n_pu * prop_missing_data)))
  pu_cand_survey <- sort(sample(pu_not_surveyed,
                                ceiling(0.5 * length(pu_not_surveyed))))
  survey_features_idx <- sort(sample.int(n_f, n_f_survey))
  survey_features_rev_idx <- rep(0, n_f)
  survey_features_rev_idx[survey_features_idx] <-
    seq_along(survey_features_idx) - 1
  survey_features <- rep(FALSE, n_f)
  survey_features[survey_features_idx] <- TRUE
  ## prior probability data
  pij <- matrix(runif(n_pu * n_f), ncol = n_pu, nrow = n_f)
  ## sensitivity and specificity data
  survey_sensitivity <- runif(n_f, 0.9, 0.99)
  survey_specificity <- runif(n_f, 0.9, 0.99)
  model_sensitivity <- matrix(runif(n_f_survey * n_o, 0.2, 0.5), ncol = n_o,
                              nrow = n_f_survey)
  total_probability_of_survey_positive <-
    (survey_sensitivity - pij) + ((1 - survey_specificity) * (1 - pij))
  total_probability_of_survey_negative <-
    ((1 - survey_sensitivity) - pij) + (survey_specificity * (1 - pij))
  total_probability_of_model_positive <-
    matrix(runif(n_f_survey * n_pu), ncol = n_pu, nrow = n_f_survey)
  ## presence/absence survey data, -1s indicate planning unit not sampled
  rij <- pij
  rij[] <- round(pij)
  rij[, pu_not_surveyed] <- -1
  ### replace 1s with prior data for species that we are not surveying
  for (i in setdiff(seq_len(n_f), survey_features_idx))
  rij[i, pu_not_surveyed] <- pij[i, pu_not_surveyed]
  ## outcome data + modelled probability data
  oij <- rij
  oij[, pu_not_surveyed] <- runif(length(pu_not_surveyed) * n_f)
  oij[, pu_cand_survey] <- round(runif(length(pu_cand_survey) * n_f))
  ## vector indicating planning units selected for surveying
  pu_survey_solution <- rep(0, n_pu)
  pu_survey_solution[pu_cand_survey] <- 1
  ## vector indicating which planning units are selected for protection
  pu_purchase_solution <- sample(c(0, 1), n_pu, replace = TRUE)
  ## outcome positions
  o <- rep(3, n_f_survey)
  ## output matrix
  out <- matrix(-100, ncol = n_pu, nrow = n_f)
  # calculations
  r1 <- rcpp_posterior_probability_matrix_wrapper(
    rij, pij, oij, pu_survey_solution,
    survey_sensitivity, survey_specificity,
    survey_features, survey_features_rev_idx,
    model_sensitivity,
    total_probability_of_survey_positive, total_probability_of_survey_negative,
    total_probability_of_model_positive, o - 1, out)
  r2 <- r_posterior_probability_matrix(
    rij, pij, oij, pu_survey_solution,
    survey_sensitivity, survey_specificity,
    survey_features, survey_features_rev_idx,
    model_sensitivity,
    total_probability_of_survey_positive, total_probability_of_survey_negative,
    total_probability_of_model_positive, o, out)
  # tests
  expect_equal(r1, r2)
})
