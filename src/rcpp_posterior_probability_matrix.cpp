#include "rcpp_posterior_probability_matrix.h"

void initialize_posterior_probability_matrix(
  Eigen::MatrixXd &pij, // prior prob data
  Eigen::MatrixXd &oij, // outcome & modelled prediction data
  std::vector<bool> &pu_survey_solution, // is the planning unit being surveyed?
  std::vector<bool> &survey_features, // is the feature being surveyed?
  std::vector<std::size_t> &survey_features_rev_idx,
  Eigen::VectorXd &survey_sensitivity,
  Eigen::VectorXd &survey_specificity,
  Eigen::MatrixXd &total_probability_of_survey_positive,
  Eigen::MatrixXd &total_probability_of_survey_negative,
  Eigen::MatrixXd &out) {
  // initialization
  const std::size_t n_pu = pij.cols();
  const std::size_t n_f = pij.rows();

  // calculate the posterior probability for each feature in each planning unit
  for (std::size_t j = 0; j < n_pu; ++j) {
    for (std::size_t i = 0; i < n_f; ++i) {
      /// probability calculation depends on various factors...
      if (!pu_survey_solution[j] || !survey_features[i]) {
        /// if the species is not being surveyed,
        /// or if the planning unit has no survey data,
        /// then use prior data
        out(i, j) = pij(i, j);
      } else if (pu_survey_solution[j] && (oij(i, j) >= 0.5)) {
        /// if the planning unit is selected for surveying,
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
        /// if the planning unit is selected for surveying,
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
      }
      /// clamp probability values to avoid numerical issues with probabilities
      /// that are eactly equal to zero or one
      out(i, j) = std::max(out(i, j), 1.0e-10);
      out(i, j) = std::min(out(i, j), 1.0 - 1.0e-10);
    }
  }

  // return void
  return;
}

void find_rij_idx_based_on_models(
  MatrixXd &model_pij, // model predictions and priors
  std::vector<std::vector<std::size_t>> &pu_model_prediction_idx,
  std::vector<bool> &survey_features, // is the feature being surveyed?
  std::vector<std::size_t> &survey_features_rev_idx, // reverse ids
  Eigen::MatrixXd &survey_tss, // TSS of survey scheme
  Eigen::VectorXd &model_sensitivity, // models sensitivity
  Eigen::VectorXd &model_specificity, // models specificity
  std::vector<std::size_t> &out) {
  // initialization
  const std::size_t n_pu = model_pij.cols();
  const std::size_t n_f = model_pij.rows();
  const std::size_t n_f_survey = survey_features_rev_idx.size();

  // declare looping variables
  std::size_t sub_i;

  // prepare output
  out.clear();
  out.reserve(n_f_survey * n_pu);

  // calculate model tss
  Eigen::VectorXd model_tss =
    model_sensitivity.array() + model_specificity.array();
  model_tss.array() -= 1.0;

  // determine which planning units will use modelled probabilities
  for (std::size_t i = 0; i < n_f_survey; ++i) {
    // extract feature idx
    sub_i = survey_features_rev_idx[i];
    for (auto itr = pu_model_prediction_idx.cbegin();
         itr != pu_model_prediction_idx.cend(); ++itr) {
      // extract planning unit idx
      j = *itr;
      if (pu_survey_solution[*itr]) {
        // if the planning unit is seleceted for additional surveys,
        // then we would only use the model predictions if they outperform
        // the survey methodology
        if (model_tss[i] >= survey_tss(i, j)) {
          curr_idx = (j * n_f) + sub_i;
        }
      } else {
        // if the planning unit is not selected for additional surveys,
        // then we must use the model predictions because no existing survey
        // data are available
        curr_idx = (j * n_f) + sub_i;
        out.push_back(curr_idx)
      }
    }
  }

  // return void
  out.shrink_to_fit();
  return;
}

void update_model_posterior_probabilities(
  std::vector<std::size_t> &rij_idx, // indices for model updates
  Eigen::MatrixXd &pij, // prior prob data
  Eigen::MatrixXd &oij, // outcome & modelled prediction data
  std::vector<std::size_t> &survey_features_rev_idx,
  Eigen::VectorXd &model_sensitivity,
  Eigen::VectorXd &model_specificity,
  Eigen::MatrixXd &total_probability_of_model_positive,
  Eigen::MatrixXd &total_probability_of_model_negative,
  Eigen::MatrixXd &out) {
  // initialization
  const std::size_t n_pu = oij.cols();
  const std::size_t n_f = oij.rows();
  std::size_t sub_i, i, j, ii;
  // calculate the posterior probability for each feature in each planning unit
  for (auto itr = rij_idx.cbegin(); itr != rij_idx.cend(); ++itr) {
    // determine species and planning unit id for ii'th index
    // (data are in column major format)
    ii = *itr;
    j = ii / n_f;
    i = ii - (j * n_f);
    // determine species reverse lookup id
    sub_i = survey_features_rev_idx[i];
    if (oij(ii) >= 0.5) {
      /// if the model predicts a presence,
      /// then the posterior probability of the species being present
      /// is based on the sensitivity of the model, our prior probaiblity,
      /// and the overall probability of the model predicting a presence
      /// accounting for false-presences
      out(ii) =
        (model_sensitivity[sub_i] * pij(ii)) /
        total_probability_of_model_positive(sub_i, j);
    } else {
      /// if the model predicts an absence,
      /// then the posterior probability of the species being present
      /// is based on the sensitivity of the model, our prior probaiblity,
      /// and the overall probability of the model predicting an absence
      /// accounting for false-absences
      out(ii) =
        ((1.0 - model_sensitivity[sub_i]) * pij(ii)) /
        total_probability_of_model_negative(sub_i, j);
    }
    /// clamp probability values to avoid numerical issues with /
    /// probabilities that are eactly equal to zero or one
    out(ii) = std::max(out(ii), 1.0e-10);
    out(ii) = std::min(out(ii), 1.0 - 1.0e-10);

    /// if posterior probability is NaN because the total
    /// probability of a model giving a positive result is zero,
    /// because it is an exceptionally poor model that cannot
    /// make any correct presence predictions, then manually set
    /// the posterior probability to a really small non-zero number
    /// (since we have to take the log of this number later
    ///  and log(0) is infinity)
    if (std::isnan(out(ii)))
      out(ii) = 1.0e-10;
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_posterior_probability_matrix(
  Eigen::MatrixXd nij,
  Eigen::MatrixXd pij,
  Eigen::MatrixXd oij,
  std::vector<bool> pu_survey_solution,
  std::vector<bool> survey_features,
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  Eigen::VectorXd model_sensitivity,
  Eigen::VectorXd model_specificity) {

  // initialize variables
  std::size_t n_f = nij.rows();
  std::size_t n_pu = nij.cols();
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

  // calculate total model probabilities
  Eigen::MatrixXd total_probability_of_model_positive(n_f_survey, n_pu);
  Eigen::MatrixXd total_probability_of_model_negative(n_f_survey, n_pu);
  total_probability_of_positive_model_result(
    pij_survey_species_subset, model_sensitivity2, model_specificity2,
    total_probability_of_model_positive);
  total_probability_of_negative_model_result(
    pij_survey_species_subset, model_sensitivity2, model_specificity2,
    total_probability_of_model_negative);

  // create default pu_model_prediction variable
  std::vector<bool> pu_model_prediction(n_pu);
  for (std::size_t i = 0; i < n_pu; ++i)
    pu_model_prediction[i] = !pu_survey_solution[i];

  // store indices of features that need surveying
  std::vector<std::size_t> survey_features_idx;
  survey_features_idx.reserve(n_f);
  for (std::size_t i = 0; i < n_f; ++i)
    if (survey_features[i])
      survey_features_idx.push_back(i);
  survey_features_idx.shrink_to_fit();

  // indices for features in sparse format for reverse lookup
  std::vector<std::size_t> survey_features_rev_idx(n_f, 1e5);
  create_reverse_lookup_id(survey_features, survey_features_rev_idx);

  // calculate TSS for survey data including new surveys
  Eigen::MatrixXd curr_survey_tss(n_f_survey, n_pu);
  calculate_survey_tss(
    nij,
    survey_features_idx,
    survey_sensitivity, survey_specificity,
    curr_survey_tss);

  // calculate posterior matrix
  Eigen::MatrixXd out(n_f, n_pu);
  initialize_posterior_probability_matrix(
    pij, oij,
    pu_survey_solution,
    survey_features, survey_features_rev_idx,
    survey_sensitivity, survey_specificity,
    total_probability_of_survey_positive,
    total_probability_of_survey_negative,
    out);

  // find rij indices to update with model estimats
  std::vector<std::size_t> curr_model_rij_idx;
  find_rij_idx_based_on_models(
    nij,
    survey_features, survey_features_rev_idx,
    pu_model_prediction,
    curr_survey_tss,
    model_sensitivity2, model_specificity2,
    curr_model_rij_idx);

  /// create posterior matrix with most likely model outcomes
  update_model_posterior_probabilities(
    curr_model_rij_idx,
    pij, oij,
    survey_features_rev_idx,
    model_sensitivity2, model_specificity2,
    total_probability_of_model_positive,
    total_probability_of_model_negative,
    out);

  // return result
  return out;
}


// [[Rcpp::export]]
Eigen::MatrixXd rcpp_initialize_posterior_probability_matrix(
  Eigen::MatrixXd nij,
  Eigen::MatrixXd pij,
  Eigen::MatrixXd oij,
  std::vector<bool> pu_survey_solution,
  std::vector<bool> survey_features,
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity) {

  // initialize variables
  std::size_t n_f = nij.rows();
  std::size_t n_pu = nij.cols();
  std::size_t n_f_survey =
    std::accumulate(survey_features.begin(), survey_features.end(), 0);

  // calculate total survey probabilities
  Eigen::MatrixXd total_probability_of_survey_positive(n_f, n_pu);
  Eigen::MatrixXd total_probability_of_survey_negative(n_f, n_pu);
  total_probability_of_positive_result(
    pij, survey_sensitivity, survey_specificity,
    total_probability_of_survey_positive);
  total_probability_of_negative_result(
    pij, survey_sensitivity, survey_specificity,
    total_probability_of_survey_negative);

  // store indices of features that need surveying
  std::vector<std::size_t> survey_features_idx;
  survey_features_idx.reserve(n_f);
  for (std::size_t i = 0; i < n_f; ++i)
    if (survey_features[i])
      survey_features_idx.push_back(i);
  survey_features_idx.shrink_to_fit();

  // indices for features in sparse format for reverse lookup
  std::vector<std::size_t> survey_features_rev_idx(n_f, 1e5);
  create_reverse_lookup_id(survey_features, survey_features_rev_idx);

  // calculate posterior matrix
  Eigen::MatrixXd out(n_f, n_pu);
  initialize_posterior_probability_matrix(
    pij, oij,
    pu_survey_solution,
    survey_features, survey_features_rev_idx,
    survey_sensitivity, survey_specificity,
    total_probability_of_survey_positive,
    total_probability_of_survey_negative,
    out);

  // return result
  return out;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_update_model_posterior_probabilities(
  Eigen::MatrixXd nij,
  Eigen::MatrixXd pij,
  Eigen::MatrixXd oij,
  std::vector<bool> pu_survey_solution,
  std::vector<bool> survey_features,
  Eigen::VectorXd survey_sensitivity,
  Eigen::VectorXd survey_specificity,
  Eigen::VectorXd model_sensitivity,
  Eigen::VectorXd model_specificity,
  Eigen::MatrixXd out) {

  // initialize variables
  std::size_t n_f = nij.rows();
  std::size_t n_pu = nij.cols();
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

  // store indices of features that need surveying
  std::vector<std::size_t> survey_features_idx;
  survey_features_idx.reserve(n_f);
  for (std::size_t i = 0; i < n_f; ++i)
    if (survey_features[i])
      survey_features_idx.push_back(i);
  survey_features_idx.shrink_to_fit();

  // indices for features in sparse format for reverse lookup
  std::vector<std::size_t> survey_features_rev_idx(n_f, 1e5);
  create_reverse_lookup_id(survey_features, survey_features_rev_idx);

  // calculate total model probabilities
  Eigen::MatrixXd total_probability_of_model_positive(n_f_survey, n_pu);
  Eigen::MatrixXd total_probability_of_model_negative(n_f_survey, n_pu);
  total_probability_of_positive_model_result(
    pij_survey_species_subset, model_sensitivity2, model_specificity2,
    total_probability_of_model_positive);
  total_probability_of_negative_model_result(
    pij_survey_species_subset, model_sensitivity2, model_specificity2,
    total_probability_of_model_negative);

  // calculate TSS for survey data including new surveys
  Eigen::MatrixXd curr_survey_tss(n_f_survey, n_pu);
  calculate_survey_tss(
    nij,
    survey_features_idx,
    survey_sensitivity, survey_specificity,
    curr_survey_tss);

  // create default pu_model_prediction variable
  std::vector<bool> pu_model_prediction(n_pu);
  for (std::size_t i = 0; i < n_pu; ++i)
    pu_model_prediction[i] = !pu_survey_solution[i];

  // find rij indices to update with model estimats
  std::vector<std::size_t> curr_model_rij_idx;
  find_rij_idx_based_on_models(
    nij,
    survey_features,
    survey_features_rev_idx,
    pu_model_prediction,
    curr_survey_tss,
    model_sensitivity2, model_specificity2,
    curr_model_rij_idx);

  /// create posterior matrix with most likely model outcomes
  update_model_posterior_probabilities(
    curr_model_rij_idx,
    pij, oij,
    survey_features_rev_idx,
    model_sensitivity2, model_specificity2,
    total_probability_of_model_positive,
    total_probability_of_model_negative,
    out);

  // return result
  return out;
}
