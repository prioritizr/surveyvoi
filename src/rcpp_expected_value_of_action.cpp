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
  Eigen::VectorXd curr_spp_probs;
  for (std::size_t i = 0; i < n_f; ++i) {
    curr_spp_probs = rij.row(i);
    spp_prob[i] = convolve_binomial(curr_spp_probs, target[i]);
  }

  // calculate probability of all species meeting targets
  return spp_prob.prod();
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
