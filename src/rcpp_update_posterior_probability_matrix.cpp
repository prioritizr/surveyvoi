#include "rcpp_update_posterior_probability_matrix.h"

void update_posterior_probability_matrix(
  Eigen::MatrixXd &pij, // prior prob data with survey outcomes
  Eigen::MatrixXd &oij, // survey outcome matrix
  std::vector<std::size_t> &idx, // indices in rij matris to update
  Eigen::VectorXd &survey_sensitivity,
  Eigen::VectorXd &survey_specificity,
  Eigen::MatrixXd &total_probability_of_survey_positive,
  Eigen::MatrixXd &total_probability_of_survey_negative,
  Eigen::MatrixXd &out) {
  // initialization
  const std::size_t n_f = pij.rows();
  std::size_t i, j, ii;
  // calculate the posterior probability for each feature in each planning unit
  // that is being survyed
  for (auto itr = idx.cbegin(); itr != idx.cend(); ++itr) {
    // calculate indices
    // set contants
    ii = *itr;
    j = ii / n_f;
    i = ii - (j * n_f);
    // update posterior probabilities
    if (oij(ii) >= 0.5) {
      /// if the planning unit is selected for surveying,
      /// and we're simulating what would happen if we detected the feature,
      ///
      /// then the posterior probability of it occuring in the planning unit
      /// is the probability that we will correctly detect the feature given
      /// our priors and the total probability that we will detect the
      /// feature accounting for false-presences
      out(ii) =
        (survey_sensitivity(i) * pij(ii)) /
        total_probability_of_survey_positive(ii);
    } else {
      /// if the planning unit is selected for surveying,
      /// and we're simulating what would happen if we didn't detect the
      //  feature,
      //
      /// then the posterior probability of it not occuring in the planning
      /// unit is one minus the overall probability that we will correctly
      /// detect the feature given our priors and the total probability
      /// that we will not detect the feature acounting for false-negatives
      out(ii) =
        ((1.0 - survey_sensitivity(i)) * pij(ii)) /
        total_probability_of_survey_negative(ii);
    }
    /// clamp probability values to avoid numerical issues with probabilities
    /// that are eactly equal to zero or one
    out(ii) = std::max(out(ii), 1.0e-10);
    out(ii) = std::min(out(ii), 1.0 - 1.0e-10);
  }

  // return void
  return;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_update_posterior_probability_matrix(
  Eigen::MatrixXd pij, // prior prob data
  Eigen::MatrixXd oij, // survey outcomes matrix
  std::vector<bool> survey_features, // features that we want to survey, 0/1
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  std::vector<bool> pu_survey_solution // planning units to survey, 0/1
) {
  /// Initialization
  /// declare constants
  const std::size_t n_pu = pij.cols();
  const std::size_t n_f = pij.rows();
  const std::size_t n_f_survey =
    std::accumulate(survey_features.begin(), survey_features.end(), 0);
  const std::size_t n_pu_surveyed_in_scheme =
    std::accumulate(pu_survey_solution.begin(), pu_survey_solution.end(), 0);

  // Preliminary processing
  /// store indices for cells in the rij matrix that will be used for
  /// simulating different outcomes
  std::vector<std::size_t> idx;
  {
    std::size_t k = 0;
    idx.reserve(n_pu_surveyed_in_scheme * n_f_survey);
    for (std::size_t j = 0; j < n_pu; ++j) {
      for (std::size_t i = 0; i < n_f; ++i) {
        if (pu_survey_solution[j] && survey_features[i]) {
          idx.push_back(k);
        }
        ++k;
      }
    }
  }

  /// calculate the total probabilities of positive and negative outcomes
  /// from the surveys
  Eigen::MatrixXd total_probability_of_survey_positive(n_f, n_pu);
  Eigen::MatrixXd total_probability_of_survey_negative(n_f, n_pu);
  total_probability_of_positive_result(
    pij, survey_sensitivity, survey_specificity,
    total_probability_of_survey_positive);
  total_probability_of_negative_result(
    pij, survey_sensitivity, survey_specificity,
    total_probability_of_survey_negative);

  // Main processing
  /// create posterior matrix with most likely model outcomes
  Eigen::MatrixXd out = pij;
  update_posterior_probability_matrix(
    pij, oij, idx,
    survey_sensitivity, survey_specificity,
    total_probability_of_survey_positive,
    total_probability_of_survey_negative,
    out);

  // Exports
  /// return result
  return out;
}
