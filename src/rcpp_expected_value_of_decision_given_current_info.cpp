#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_expected_value_of_action.h"

// [[Rcpp::export]]
double rcpp_expected_value_of_decision_given_current_info(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap) {

  // find optimal management action using prior data
  std::vector<bool> solution(pij.cols());
  Prioritization p(pij.cols(), pij.rows(), pu_costs, pu_locked_in,
                   alpha, gamma, n_approx_obj_fun_points, budget, gap);
  p.add_rij_data(pij);
  p.solve();
  p.get_solution(solution);

  // calculate log prior probabilities
  Eigen::MatrixXd pij_log1m = pij;
  log_1m_matrix(pij_log1m);
  log_matrix(pij);

  // calculate expected value of management action
  return expected_value_of_action(solution, pij, pij_log1m, alpha, gamma);
}
