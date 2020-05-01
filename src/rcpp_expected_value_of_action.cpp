#include "rcpp_expected_value_of_action.h"

double expected_value_of_action(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &preweight,
  Eigen::VectorXd &postweight,
  Eigen::VectorXd &target) {

  // initialization
  const std::size_t n_pu = pij.cols();
  const double total = static_cast<double>(n_pu);

  // prepare rij matrix containing only solution values
  Eigen::MatrixXd rij = pij;
  for (std::size_t i = 0; i < n_pu; ++i)
    rij.col(i) *= solution[i];

  // return convervation value
  return conservation_value_state(rij, preweight, postweight, target, total);
}

// [[Rcpp::export]]
double rcpp_expected_value_of_action(
  std::vector<bool> solution,
  Eigen::MatrixXd pij,
  Eigen::VectorXd preweight,
  Eigen::VectorXd postweight,
  Eigen::VectorXd target) {
  // return result
  return expected_value_of_action(
    solution, pij, preweight, postweight, target);
}
