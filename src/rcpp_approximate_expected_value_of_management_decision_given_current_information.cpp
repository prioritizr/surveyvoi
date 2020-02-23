#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_approximate_expected_value_of_management_action.h"

double approximate_expected_value_of_management_decision_given_current_information(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap,
  std::vector<mpz_t> &states) {
  // find optimal management action using prior data
  std::vector<bool> solution(pij.cols());
  Prioritization p(pij.cols(), pij.rows(), pu_costs, pu_locked_in,
                   alpha, gamma, n_approx_obj_fun_points, budget, gap);
  p.add_rij_data(pij);
  p.solve();
  p.get_solution(solution);

  // calculate log prior probabilities
  pij.array() = pij.array().log();

  // calculate expected value of management action
  double out = approximate_expected_value_of_management_action(
    solution, pij, alpha, gamma, states);

  // return result
  return out;
}

// [[Rcpp::export]]
double rcpp_approximate_expected_value_of_management_decision_given_current_information_n_states(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap,
  std::size_t n_approx_states) {

  // generate states
  std::vector<mpz_t> states(n_approx_states);
  sample_k_nth_states(n_approx_states, pij, states);

  // calculate result
  double out = approximate_expected_value_of_management_decision_given_current_information(
    pij, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points, budget,
    gap, states);

  // clear memory
  for (std::size_t i = 0; i < n_approx_states; ++i)
    mpz_clear(states[i]);

  // return result
  return out;
}

// [[Rcpp::export]]
double rcpp_approximate_expected_value_of_management_decision_given_current_information_fixed_states(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap,
  std::vector<std::size_t> states) {
  // convert state indices from std::size_t to mpz_t
  const std::size_t n = states.size();
  std::vector<mpz_t> states2(n);
  for (std::size_t i = 0; i < n; ++i) {
    mpz_init(states2[i]);
    mpz_set_ui(states2[i], states[i]);
  }

  // calculate result
  double out = approximate_expected_value_of_management_decision_given_current_information(
    pij, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points, budget,
    gap, states2);

  // clear memory
  for (std::size_t i = 0; i < n; ++i)
    mpz_clear(states2[i]);

  // return result
  return out;
}
