#include "rcpp_prioritization.h"

// [[Rcpp::export]]
Rcpp::List rcpp_prioritization(
  Eigen::MatrixXd rij,
  Eigen::VectorXd pu_costs,
  Eigen::VectorXd pu_locked_in,
  Eigen::VectorXd preweight,
  Eigen::VectorXd postweight,
  Eigen::VectorXd target,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap,
  std::string file_path) {
  // initialization
  std::vector<bool> solution(rij.cols());
  // main processing
  Prioritization p(rij.cols(), rij.rows(), pu_costs, pu_locked_in,
                   preweight, postweight, target, n_approx_obj_fun_points,
                   budget, gap);
  p.add_rij_data(rij);
  p.solve();
  p.get_solution(solution);
  if (file_path.size() > 0)
    p.save_problem(file_path);
  // export
  return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(solution),
                            Rcpp::Named("objval") = p.get_obj_value());
}
