#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_heuristic_prioritization.h"
#include "rcpp_expected_value_of_action.h"

// [[Rcpp::export]]
double rcpp_expected_value_of_decision_given_current_info(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &pu_locked_out,
  Eigen::VectorXi &target,
  double budget) {

  // find optimal management action using prior data
  std::vector<bool> solution(pij.cols());
  greedy_heuristic_prioritization(
    pij, pu_costs, pu_locked_in, pu_locked_out, target, budget, solution);

  // calculate expected value of management action
  return expected_value_of_action(solution, pij, target);
}
