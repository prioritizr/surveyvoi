#include "rcpp_probability.h"

void total_probability_of_positive_result(
  Eigen::MatrixXd &prior, Eigen::VectorXd &sensitivity,
  Eigen::VectorXd &specificity, Eigen::MatrixXd &out) {
  const std::size_t n_pu = prior.cols();
  for (std::size_t j = 0; j < n_pu; ++j)
    out.col(j) = (sensitivity.array() * prior.col(j).array()) +
                 ((1.0 - specificity.array()) * (1.0 - prior.col(j).array()));
  out = out.cwiseMax(1.0e-10);
  out = out.cwiseMin(1.0 - 1.0e-10);
  return;
}

void total_probability_of_negative_result(
  Eigen::MatrixXd &prior, Eigen::VectorXd &sensitivity,
  Eigen::VectorXd &specificity, Eigen::MatrixXd &out) {
  const std::size_t n_pu = prior.cols();
  for (std::size_t j = 0; j < n_pu; ++j)
    out.col(j) = (((1.0 - sensitivity.array())) * prior.col(j).array()) +
                 (specificity.array() * (1.0 - prior.col(j).array()));
  out = out.cwiseMax(1.0e-10);
  out = out.cwiseMin(1.0 - 1.0e-10);
  return;
}

void total_probability_of_positive_model_result(
  Eigen::MatrixXd &prior, // only contains rows for spp that need surveys
  Eigen::VectorXd &sensitivity, // only elements for spp needing surveys
  Eigen::VectorXd &specificity, // only elements for spp needing surveys
  Eigen::MatrixXd &out) {
  const std::size_t n_pu = prior.cols();
  const std::size_t n_f = sensitivity.size();
  for (std::size_t j = 0; j < n_pu; ++j) {
    for (std::size_t i = 0; i < n_f; ++i) {
      out(i, j) =
        (sensitivity[i] * prior(i, j)) +
        ((1.0 - specificity[i]) * (1.0 - prior(i, j)));
    }
  }
  out = out.cwiseMax(1.0e-10);
  out = out.cwiseMin(1.0 - 1.0e-10);
  return;
}

void total_probability_of_negative_model_result(
  Eigen::MatrixXd &prior, // only contains rows for spp that need surveys
  Eigen::VectorXd &sensitivity, // only elements for spp needing surveys
  Eigen::VectorXd &specificity, // only elements for spp needing surveys
  Eigen::MatrixXd &out) {
  const std::size_t n_pu = prior.cols();
  const std::size_t n_f = sensitivity.size();
  for (std::size_t j = 0; j < n_pu; ++j) {
    for (std::size_t i = 0; i < n_f; ++i) {
      out(i, j) =
        ((1.0 - sensitivity[i]) * prior(i, j)) +
        (specificity[i] * (1.0 - prior(i, j)));
    }
  }
  out = out.cwiseMax(1.0e-10);
  out = out.cwiseMin(1.0 - 1.0e-10);
  return;
}

double log_probability_of_outcome(
  Eigen::MatrixXd &oij,
  Eigen::MatrixXd &total_probability_of_survey_positive_log,
  Eigen::MatrixXd &total_probability_of_survey_negative_log,
  std::vector<std::size_t> &idx) {
  // init
  double out = 0.0;
  // main
  for (auto itr = idx.cbegin(); itr != idx.cend(); ++itr)
    out += (oij(*itr) * total_probability_of_survey_positive_log(*itr)) +
           ((1.0 - oij(*itr)) * total_probability_of_survey_negative_log(*itr));
  // return void
  return out;
}

double log_probability_of_model_outcome(
  Eigen::MatrixXd &oij,
  std::vector<std::size_t> &survey_features_rev_idx,
  Eigen::MatrixXd &total_probability_of_model_positive_log,
  Eigen::MatrixXd &total_probability_of_model_negative_log,
  std::vector<std::size_t> &idx) {
  // init
  double out = 0.0;
  const std::size_t n_pu = oij.cols();
  const std::size_t n_f = oij.rows();
  std::size_t i, j, ii, sub_i;
  // main
  for (auto itr = idx.cbegin(); itr != idx.cend(); ++itr) {
    // determine species and planning unit indices from rij index
    // (data are in column major format)
    ii = *itr;
    j = ii / n_f;
    i = ii - (j * n_f);
    // determine species reverse lookup id from species index
    sub_i = survey_features_rev_idx[i];
    // add probabilities
    out +=
      (oij(ii) * total_probability_of_model_positive_log(sub_i, j)) +
      ((1.0 - oij(ii)) * total_probability_of_model_negative_log(sub_i, j));
  }
  // return void
  return out;
}

double log_probability_of_state(
  Eigen::MatrixXd &sij,
  Eigen::MatrixXd &pij_log,
  Eigen::MatrixXd &pij_log1m) {
  return
    (sij.array() * pij_log.array()).sum() +
    ((1.0 - sij.array()) * pij_log1m.array()).sum();
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_total_probability_of_positive_result(
  Eigen::MatrixXd prior,
  Eigen::VectorXd sensitivity,
  Eigen::VectorXd specificity) {
  Eigen::MatrixXd out(prior.cols(), prior.rows());
  total_probability_of_positive_result(prior, sensitivity, specificity, out);
  return out;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_total_probability_of_negative_result(
  Eigen::MatrixXd prior,
  Eigen::VectorXd sensitivity,
  Eigen::VectorXd specificity) {
  Eigen::MatrixXd out(prior.cols(), prior.rows());
  total_probability_of_negative_result(prior, sensitivity, specificity, out);
  return out;
}

// [[Rcpp::export]]
double rcpp_probability_of_outcome(
  Eigen::MatrixXd oij,
  Eigen::MatrixXd total_probability_of_survey_positive,
  Eigen::MatrixXd total_probability_of_survey_negative,
  std::vector<std::size_t> idx) {
  // preapre data for calculations
  for(auto& i : idx)
    i -= 1;
  log_matrix(total_probability_of_survey_positive);
  log_matrix(total_probability_of_survey_negative);
  // calculate probability
  return std::exp(log_probability_of_outcome(
    oij, total_probability_of_survey_positive,
    total_probability_of_survey_negative, idx));
}

// [[Rcpp::export]]
double rcpp_probability_of_state(
  Eigen::MatrixXd sij,
  Eigen::MatrixXd pij) {
  Eigen::MatrixXd pij_log1m = pij;
  log_1m_matrix(pij_log1m);
  log_matrix(pij);
  return std::exp(log_probability_of_state(sij, pij, pij_log1m));
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_total_probability_of_positive_model_result(
  Eigen::MatrixXd prior,
  Eigen::VectorXd sensitivity,
  Eigen::VectorXd specificity) {
  Eigen::MatrixXd out(prior.cols(), prior.rows());
  total_probability_of_positive_model_result(
    prior, sensitivity, specificity, out);
  return out;
}

// [[Rcpp::export]]
Eigen::MatrixXd rcpp_total_probability_of_negative_model_result(
  Eigen::MatrixXd prior,
  Eigen::VectorXd sensitivity,
  Eigen::VectorXd specificity) {
  Eigen::MatrixXd out(prior.cols(), prior.rows());
  total_probability_of_negative_model_result(
    prior, sensitivity, specificity, out);
  return out;
}
