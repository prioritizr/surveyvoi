#include "package.h"

// [[Rcpp::export]]
Rcpp::List rcpp_pmedian_constraint_matrix(
  Eigen::MatrixXd x, Rcpp::NumericVector costs) {
  // initialization
  std::size_t n = x.rows();
  std::size_t n_vars = n + (n * n);
  std::size_t n_reserve = n + (n * n) + (n * n) + (n * n);
  // create constraints matrix
  std::vector<double> A_i;
  std::vector<double> A_j;
  std::vector<double> A_x;
  A_i.reserve(n_reserve);
  A_j.reserve(n_reserve);
  A_x.reserve(n_reserve);

  // add constraints to ensure that selected sites are within budget
  for (std::size_t i = 0; i < n; ++i) {
    A_i.push_back(0);
    A_j.push_back(i);
    A_x.push_back(costs[i]);
  }
  // add constraints to ensure each point is assigned to a single selected point
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      A_i.push_back(i + 1);
      A_j.push_back(n + (i * n) + j);
      A_x.push_back(1.0);
    }
  }
  // add constraints to ensure that points can only be assigned to selected
  // points
  std::size_t counter = n + 1;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      A_i.push_back(counter + j);
      A_j.push_back(i);
      A_x.push_back(-1.0);
    }
    for (std::size_t j = 0; j < n; ++j) {
      A_i.push_back(counter + j);
      A_j.push_back(n + (j * n) + i);
      A_x.push_back(1.0);
    }
    counter += n;
  }
  // return result
  return Rcpp::List::create(Rcpp::Named("i") = A_i,
                            Rcpp::Named("j") = A_j,
                            Rcpp::Named("x") = A_x);
}
