#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_approx_expected_value_of_action.h"

// [[Rcpp::export]]
Rcpp::NumericVector
  rcpp_approx_expected_value_of_decision_given_current_info_n_states(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap,
  std::size_t n_approx_replicates,
  std::size_t n_approx_states_per_replicate) {

  // initialize
  Rcpp::NumericVector out(n_approx_replicates);
  std::vector<mpz_t> states(n_approx_states_per_replicate);
  for (std::size_t j = 0; j < n_approx_states_per_replicate; ++j)
    mpz_init(states[j]);

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

  // main processing
  for (std::size_t i = 0; i < n_approx_replicates; ++i) {
    /// generate states
    sample_k_uniform_no_replacement_nth_states(
      n_approx_states_per_replicate, pij, states);
    /// calculate result
    out[i] = approx_expected_value_of_action(
        solution, pij, pij_log1m, alpha, gamma, states);
  }

  // clear memory
  for (std::size_t i = 0; i < n_approx_states_per_replicate; ++i)
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
  // initialize states
  const std::size_t n = states.size();
  std::vector<mpz_t> states2(n);
  for (std::size_t i = 0; i < n; ++i) {
    mpz_init(states2[i]);
    mpz_set_ui(states2[i], states[i]);
  }

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

  // calculate result
  double out = approx_expected_value_of_action(
    solution, pij, pij_log1m, alpha, gamma, states2);

  // clear memory
  for (std::size_t i = 0; i < n; ++i)
    mpz_clear(states2[i]);

  // return result
  return out;
}
