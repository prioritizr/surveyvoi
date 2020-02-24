#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_approx_expected_value_of_action.h"

double approx_expected_value_of_decision_given_current_info(
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
  Eigen::MatrixXd pij_log1m = pij;
  pij_log1m.array() = (1.0 - pij_log1m.array()).array().log();
  pij.array() = pij.array().log();

  // calculate expected value of management action
  double out = approx_expected_value_of_action(
    solution, pij, pij_log1m, alpha, gamma, states);

  // return result
  return out;
}

// [[Rcpp::export]]
double rcpp_approx_expected_value_of_decision_given_current_info_n_states(
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
  sample_k_uniform_nth_states(n_approx_states, pij, states);

  // calculate result
  double out = approx_expected_value_of_decision_given_current_info(
    pij, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points, budget,
    gap, states);

  // clear memory
  for (std::size_t i = 0; i < n_approx_states; ++i)
    mpz_clear(states[i]);

  // return result
  return out;
}

// [[Rcpp::export]]
double rcpp_approx_expected_value_of_decision_given_current_info_fixed_states(
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
  double out = approx_expected_value_of_decision_given_current_info(
    pij, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points, budget,
    gap, states2);

  // clear memory
  for (std::size_t i = 0; i < n; ++i)
    mpz_clear(states2[i]);

  // return result
  return out;
}
