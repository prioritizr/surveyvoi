#include "rcpp_posterior_probability_matrix.h"

void initialize_posterior_probability_matrix(
  Eigen::MatrixXd &nij, // number of existing survey data
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
  const std::size_t n_pu = nij.cols();
  const std::size_t n_f = nij.rows();

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
  Eigen::MatrixXd &nij, // number of existing survey data
  std::vector<bool> &pu_survey_solution, // is the planning unit being surveyed?
  std::vector<bool> &survey_features, // is the feature being surveyed?
  std::vector<std::size_t> &survey_features_rev_idx,
  Eigen::VectorXd &survey_sensitivity,
  Eigen::VectorXd &survey_specificity,
  Eigen::VectorXd &model_sensitivity,
  Eigen::VectorXd &model_specificity,
  std::vector<std::size_t> &out) {
  // initialization
  const std::size_t n_pu = nij.cols();
  const std::size_t n_f = nij.rows();
  std::size_t counter = 0;
  out.clear();
  out.reserve(
    std::accumulate(survey_features.begin(), survey_features.end(), 0) * n_pu);

  // determine which planning units will use modelled probabilities
  for (std::size_t j = 0; j < n_pu; ++j) {
    for (std::size_t i = 0; i < n_f; ++i) {
      if ((nij(i, j) < 0.5) &&
          (survey_features[i]) &&
          (!pu_survey_solution[j])) {
        out.push_back(counter);
      }
      ++counter;
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
    ii = *itr;
    i = ii / n_pu;
    j = ii - (i * n_pu);
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

      print(wrap("pij"));
      print(wrap(pij));
      print(wrap("model_specificity"));
      print(wrap(model_specificity));
      print(wrap("model_sensitivity"));
      print(wrap(model_sensitivity));
      print(wrap("total_probability_of_model_positive"));
      print(wrap(total_probability_of_model_positive));
      print(wrap("total_probability_of_model_negative"));
      print(wrap(total_probability_of_model_negative));

      print(wrap("ii"));
      print(wrap(ii));
      print(wrap("i"));
      print(wrap(i));
      print(wrap("j"));
      print(wrap(j));



      out(ii) =
        ((1.0 - model_sensitivity[sub_i]) * pij(ii)) /
        total_probability_of_model_negative(sub_i, j);


      print(wrap("out(ii)"));
      print(wrap(out(ii)));

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
  initialize_posterior_probability_matrix(
    nij, pij, oij,
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
    pu_survey_solution,
    survey_features, survey_features_rev_idx,
    survey_sensitivity, survey_specificity,
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

  // calculate posterior matrix
  Eigen::MatrixXd out(n_f, n_pu);
  initialize_posterior_probability_matrix(
    nij, pij, oij,
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

  // find rij indices to update with model estimats
  std::vector<std::size_t> curr_model_rij_idx;
  find_rij_idx_based_on_models(
    nij,
    pu_survey_solution,
    survey_features, survey_features_rev_idx,
    survey_sensitivity, survey_specificity,
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
