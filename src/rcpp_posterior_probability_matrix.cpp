#include "rcpp_posterior_probability_matrix.h"

void posterior_probability_matrix(
  Eigen::MatrixXd &rij, // original rij data
  Eigen::MatrixXd &pij, // prior prob data
  Eigen::MatrixXd &oij, // outcome & modelled prediction data
  std::vector<bool> &pu_survey_solution, // is the planning unit being surveyed?
  std::vector<bool> &survey_features, // is the feature being surveyed?
  std::vector<std::size_t> &survey_features_rev_idx,
  Eigen::VectorXd &survey_sensitivity,
  Eigen::VectorXd &survey_specificity,
  Eigen::MatrixXd &total_probability_of_survey_positive,
  Eigen::MatrixXd &total_probability_of_survey_negative,
  Eigen::VectorXd &curr_model_sensitivity,
  Eigen::VectorXd &curr_model_specificity,
  Eigen::MatrixXd &total_probability_of_model_positive,
  Eigen::MatrixXd &total_probability_of_model_negative,
  Eigen::MatrixXd &out) {

  // initialization
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = rij.rows();
  std::size_t sub_i;

  // calculate the posterior probability for each feature in each planning unit
  for (std::size_t j = 0; j < n_pu; ++j) {
    for (std::size_t i = 0; i < n_f; ++i) {
      /// probability calculation depends on various factors...
      if (rij(i, j) > -0.5) {
        /// if the planning unit has already been surveyed,
        /// then use the prior
        out(i, j) = pij(i, j);
      } else if (!survey_features[i]) {
        /// if the planning unit has not already been surveyed,
        /// and the feature isn't relevant for future surveys,
        ///
        /// then the posterior probability of it occuring in the planning unit
        /// is simply the prior probability because we aren't learning any
        /// new information about it
        out(i, j) = pij(i, j);
      } else if ((pu_survey_solution[j]) && (oij(i, j) >= 0.5)) {
        /// if the planning unit has not already been surveyed,
        /// and is selected for surveying,
        /// and we're simulating what would happen if we detected the feature,
        ///
        /// then the posterior probability of it occuring in the planning unit
        /// is the probability that we will correctly detect the feature given
        /// our priors and the total probability that we will detect the
        /// feature accounting for false-presences
        out(i, j) =
          (survey_sensitivity(i) * pij(i, j)) /
          total_probability_of_survey_positive(i, j);
      } else if (pu_survey_solution[j]) {
        /// if the planning unit has not already been surveyed,
        /// and is selected for surveying,
        /// and we're simulating what would happen if we didn't detect the
        //  feature,
        //
        /// then the posterior probability of it not occuring in the planning
        /// unit is one minus the overall probability that we will correctly
        /// detect the feature given our priors and the total probability
        /// that we will not detect the feature acounting for false-negatives
        out(i, j) =
          ((1.0 - survey_sensitivity(i)) * pij(i, j)) /
          total_probability_of_survey_negative(i, j);
      } else {
        /// if the planning unit has not already been surveyed,
        /// and is not selected for surveying,
        /// and we're simulating the probability of the species being
        /// present based on the model,
        sub_i = survey_features_rev_idx[i];
        if (oij(i, j) >= 0.5) {
          /// if the model predicts a presence,
          /// then the posterior probability of the species being present
          /// is based on the sensitivity of the model, our prior probaiblity,
          /// and the overall probability of the model predicting a presence
          /// accounting for false-presences
          out(i, j) =
            (curr_model_sensitivity[sub_i] * pij(i, j)) /
            total_probability_of_model_positive(sub_i, j);
        } else {
          /// if the model predicts an absence,
          /// then the posterior probability of the species being present
          /// is based on the sensitivity of the model, our prior probaiblity,
          /// and the overall probability of the model predicting an absence
          /// accounting for false-absences
          out(i, j) =
            ((1.0 - curr_model_sensitivity[sub_i]) * pij(i, j)) /
            total_probability_of_model_negative(sub_i, j);
        }
        /// if posterior probability is NaN because the total
        /// probability of a model giving a positive result is zero,
        /// because it is an exceptionally poor model than that cannot
        /// make any correct presence predictions, then manually set
        /// the posterior probability to a really small non-zero number
        /// (since we have to take the log of this number later
        ///  and log(0) is infinity)
        if (std::isnan(out(i, j)))
          out(i, j) = 1.0e-15;
      }
    }
  }

  // return void
  return;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_posterior_probability_matrix(
  Eigen::MatrixXd rij,
  Eigen::MatrixXd pij,
  Eigen::MatrixXd oij,
  std::vector<bool> pu_survey_solution,
  std::vector<bool> survey_features,
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  Eigen::VectorXd model_sensitivity,
  Eigen::VectorXd model_specificity) {

  // initialize variables
  std::size_t n_f = rij.rows();
  std::size_t n_pu = rij.cols();
  std::size_t n_f_survey =
    std::accumulate(survey_features.begin(), survey_features.end(), 0);

  // subset prior probabilities to only include surveyed features
  Eigen::MatrixXd pij_survey_species_subset(n_f_survey, n_pu);
  {
    std::size_t k = 0;
    for (std::size_t i = 0; i < n_f; ++i) {
      if (survey_features[i]) {
        pij_survey_species_subset.row(k) = pij.row(i);
        ++k;
      }
    }
  }

  // subset model performance to only include surveyed features
  Eigen::VectorXd model_sensitivity2(n_f_survey);
  Eigen::VectorXd model_specificity2(n_f_survey);
  {
    std::size_t k = 0;
    for (std::size_t i = 0; i < n_f; ++i) {
      if (survey_features[i]) {
        model_sensitivity2[k] = model_sensitivity[i];
        model_specificity2[k] = model_specificity[i];
        ++k;
      }
    }
  }

  // calculate total survey probabilities
  Eigen::MatrixXd total_probability_of_survey_positive(n_f, n_pu);
  Eigen::MatrixXd total_probability_of_survey_negative(n_f, n_pu);
  total_probability_of_positive_result(
    pij, survey_sensitivity, survey_specificity,
    total_probability_of_survey_positive);
  total_probability_of_negative_result(
    pij, survey_sensitivity, survey_specificity,
    total_probability_of_survey_negative);

  // prepare feature ids
  std::vector<std::size_t> survey_features_rev_idx(n_f, 0);
  {
    std::size_t k = 0;
    for (std::size_t i = 0; i < n_f; ++i) {
      if (survey_features[i]) {
        survey_features_rev_idx[i] = k;
        ++k;
      }
    }
  }

  // calculate total model probabilities
  Eigen::MatrixXd total_probability_of_model_positive(n_f_survey, n_pu);
  Eigen::MatrixXd total_probability_of_model_negative(n_f_survey, n_pu);
  total_probability_of_positive_model_result(
    pij_survey_species_subset, model_sensitivity2, model_specificity2,
    total_probability_of_model_positive);
  total_probability_of_negative_model_result(
    pij_survey_species_subset, model_sensitivity2, model_specificity2,
    total_probability_of_model_negative);

  // calculate posterior matrix
  Eigen::MatrixXd out(n_f, n_pu);
  posterior_probability_matrix(
    rij, pij, oij,
    pu_survey_solution,
    survey_features, survey_features_rev_idx,
    survey_sensitivity, survey_specificity,
    total_probability_of_survey_positive,
    total_probability_of_survey_negative,
    model_sensitivity2, model_specificity2,
    total_probability_of_model_positive,
    total_probability_of_model_negative, out);
  return out;
}
