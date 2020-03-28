#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"

// [[Rcpp::export]]
double rcpp_expected_value_of_decision_given_perfect_info(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap) {

  // initialization
  /// initialize loop variables
  const std::size_t n_pu = pij.cols();
  double out = std::numeric_limits<double>::infinity();
  double curr_value_given_state_occurring;
  double curr_probability_of_state_occurring;
  double curr_expected_value_given_state;
  Eigen::MatrixXd curr_state(pij.rows(), pij.cols());
  Eigen::MatrixXd curr_rij(pij.rows(), pij.cols());

  /// create log version of probabilities
  Eigen::MatrixXd pij_log = pij;
  Eigen::MatrixXd pij_log1m = pij;
  log_matrix(pij_log);
  log_1m_matrix(pij_log1m);

  /// initialize prioritization
  std::vector<bool> solution(n_pu);
  Prioritization p(pij.cols(), pij.rows(), pu_costs, pu_locked_in,
                   alpha, gamma, n_approx_obj_fun_points, budget, gap);

  /// determine number of states
  mpz_class n;
  n_states(curr_state.size(), n);
  n = n + 1;

  /// initialize loop iterator
  mpz_class i = 1;

  // main processing
  while (cmp(i, n) < 0) {

    /// generate the i'th state
    nth_state(i, curr_state);
    /// generate solution for state
    p.add_rij_data(curr_state);
    p.solve();
    p.get_solution(solution);
    /// create matrix only containing feature data for selected planning units
    curr_rij = curr_state;
    for (std::size_t j = 0; j < n_pu; ++j)
      curr_rij.col(j).array() *= static_cast<double>(solution[j]);
    /// calculate the value of the prioritization given the state
    curr_value_given_state_occurring =
      alpha.cwiseProduct(curr_rij.rowwise().sum()).array().
        pow(gamma.array()).sum();
    /// if prioritization has a non-zero value then proceed with remaining
    /// calculations for this state
    if (curr_value_given_state_occurring > 1.0e-10) {
      /// calculate probability of the state occurring
      curr_probability_of_state_occurring =
        log_probability_of_state(curr_state, pij_log, pij_log1m);
      /// add the value of the prioritization given the state,
      /// weighted by the probability of the state occuring
      curr_expected_value_given_state =
        std::log(curr_value_given_state_occurring) +
        curr_probability_of_state_occurring;
      /// calculate expected value of action
      if (std::isinf(out)) {
        out = curr_expected_value_given_state;
      } else {
        out = log_sum(out, curr_expected_value_given_state);
      }
    }

    /// increment loop variable
    i = i + 1;
  }

  // return result
  return std::exp(out);
}
