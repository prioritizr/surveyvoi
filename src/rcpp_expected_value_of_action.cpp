#include "rcpp_expected_value_of_action.h"

double expected_value_of_action(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij,
  Eigen::VectorXi &target) {
  // initialization
  const std::size_t n_pu = pij.cols();
  const std::size_t n_f = pij.rows();
  const std::size_t n_pu_solution = std::accumulate(solution.begin(),
                                                    solution.end(), 0);
  if (static_cast<int>(n_pu_solution) < target.maxCoeff()) {
    Rcpp::stop("prioritization contains fewer planning units than a target");
  }

  // prepare rij matrix containing only solution values
  std::size_t j = 0;
  Eigen::MatrixXd rij(n_f, n_pu_solution);
  for (std::size_t i = 0; i < n_pu; ++i) {
    if (solution[i]) {
      rij.col(j) = pij.col(i);
      ++j;
    }
  }

  // calculate probability of each species meeting the targets
  Eigen::VectorXd spp_prob(n_f);
  Rcpp::NumericVector curr_spp_probs;
  Rcpp::NumericVector curr_density_probs;
  Rcpp::IntegerVector curr_target_values((rij.cols() - target[0]) + 1);
  std::iota(curr_target_values.begin(), curr_target_values.end(), target[0]);
  for (std::size_t i = 0; i < n_f; ++i) {
    curr_spp_probs = wrap(rij.row(i));
    curr_density_probs = PoissonBinomial::dpb_dc(
      curr_target_values, curr_spp_probs);
    spp_prob[i] = Rcpp::sum(curr_density_probs);
  }

  // calculate probability of all species meeting targets
  return spp_prob.prod();
}

double approx_expected_value_of_action(
  Eigen::MatrixXd &pij,
  Rcpp::IntegerVector &target_values) {
  double out = 1.0;
  const std::size_t n_f = pij.rows();
  Rcpp::NumericVector curr_density_probs;
  Rcpp::NumericVector curr_spp_probs;
  for (std::size_t i = 0; i < n_f; ++i) {
    curr_spp_probs = wrap(pij.row(i));
    curr_density_probs = PoissonBinomial::dpb_na(
      target_values, curr_spp_probs, false);
    out *= Rcpp::sum(curr_density_probs);
  }
  return out;
}

double log_proxy_expected_value_of_action(
  Eigen::MatrixXd &log_1m_pij) {
  Eigen::VectorXd out = log_1m_pij.rowwise().sum();
  for (auto itr = out.data(); itr != out.data() + out.size(); ++itr)
    *itr = log_substract(0.0, *itr);
  return out.sum();
}

// [[Rcpp::export]]
double rcpp_expected_value_of_action(
  std::vector<bool> solution,
  Eigen::MatrixXd pij,
  Eigen::VectorXi target) {
  // return result
  return expected_value_of_action(
    solution, pij, target);
}
