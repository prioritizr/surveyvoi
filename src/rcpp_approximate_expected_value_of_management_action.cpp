#include "rcpp_approximate_expected_value_of_management_action.h"

double approximate_expected_value_of_prioritization(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij_log,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::vector<mpz_t> &states) {
  // initialization
  const std::size_t n_pu = pij_log.cols();
  const std::size_t solution_size =
    std::accumulate(solution.begin(), solution.end(), 0);
  Eigen::MatrixXd sub_pij_log(pij_log.rows(), solution_size);

  // preliminary processing
  /// subset planning units selected in the solution
  {
    std::size_t i = 0;
    for (std::size_t j = 0; j < n_pu; ++j) {
      if (solution[j]) {
        sub_pij_log.col(i) = pij_log.col(j);
        ++i;
      }
    }
  }

  // main processing
  /// determine number of states that affect the solution
  const std::size_t n_approx_states = states.size();
  /// initialize loop variables
  Eigen::VectorXd curr_value_given_state_occurring(n_approx_states);
  Eigen::VectorXd curr_probability_of_state_occurring(n_approx_states);
  Eigen::MatrixXd curr_state(sub_pij_log.rows(), sub_pij_log.cols());
  /// iterate over each state
  for (std::size_t i = 0; i < n_approx_states; ++i) {
    //// generate the i'th state
    nth_state(states[i], curr_state);
    //// calculate the value of the prioritization given the state
    curr_value_given_state_occurring[i] =
      std::log(alpha.cwiseProduct(curr_state.rowwise().sum()).array().
        pow(gamma.array()).sum());
    /// calculate probability of the state occurring
    curr_probability_of_state_occurring[i] =
      log_probability_of_state(curr_state, sub_pij_log);
  }

  // rescale probabilities
  curr_probability_of_state_occurring.array() /=
    curr_probability_of_state_occurring.sum();

  // calculate values weighted by probabilities
  curr_probability_of_state_occurring.array() +=
    curr_value_given_state_occurring.array();

  // return result
  return std::exp(log_sum(curr_probability_of_state_occurring));
}

// [[Rcpp::export]]
double rcpp_approximate_expected_value_of_prioritization(
  std::vector<bool> solution,
  Eigen::MatrixXd pij,
  Eigen::VectorXd alpha,
  Eigen::VectorXd gamma,
  std::vector<std::size_t> states) {
  // calculate log pij
  pij.array() = pij.array().log();

  // convert state indices from std::size_t to mpz_t
  const std::size_t n = states.size();
  std::vector<mpz_t> states2(n);
  for (std::size_t i = 0; i < n; ++i) {
    mpz_init(states2[i]);
    mpz_set_ui(states2[i], states[i]);
  }

  // calculate result
  double out = approximate_expected_value_of_prioritization(
    solution, pij, alpha, gamma, states2);

  // clear memory
  for (std::size_t i = 0; i < n; ++i)
    mpz_clear(states2[i]);

  // return result
  return out;
}
