#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_approx_expected_value_of_action.h"

// [[Rcpp::export]]
double rcpp_approx_expected_value_of_decision_given_current_info(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &pu_locked_out,
  Eigen::VectorXd &preweight,
  Eigen::VectorXd &postweight,
  Eigen::VectorXd &target,
  double budget,
  double gap,
  std::size_t n_approx_states) {

  /// clamp number of approximation states to total number of states
  mpz_class n_states_total;
  mpz_class n_approx_states2 = n_approx_states;
  if ((pij.cols() * pij.rows()) < 30) {
    n_states(pij.cols() * pij.rows(), n_states_total);
    n_states_total = n_states_total + 1; // increment to include all states
    if (cmp(n_approx_states2, n_states_total) > 0)
      n_approx_states = n_states_total.get_ui();
  }

  // find optimal management action using prior data
  std::vector<bool> solution(pij.cols());
  Prioritization p(pij.cols(), pij.rows(), pu_costs,
                   pu_locked_in, pu_locked_out,
                   preweight, postweight, target, budget, gap);

  // generate solution
  p.add_rij_data(pij);
  p.solve();
  p.get_solution(solution);

  // calculate log matrices
  Eigen::MatrixXd pij_log = pij;
  Eigen::MatrixXd pij_log1m = pij;
  log_matrix(pij_log);
  log_1m_matrix(pij_log1m);

  // calculate expected value of management action
  return approx_expected_value_of_action(
    solution, pij, pij_log, pij_log1m, preweight, postweight, target,
    n_approx_states);
}
