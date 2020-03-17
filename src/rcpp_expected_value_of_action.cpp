#include "rcpp_expected_value_of_action.h"

double expected_value_of_action(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij_log,
  Eigen::MatrixXd &pij_log1m,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma) {

  // initialization
  const std::size_t n_pu = pij_log.cols();
  const std::size_t solution_size =
    std::accumulate(solution.begin(), solution.end(), 0);
  Eigen::MatrixXd sub_pij_log(pij_log.rows(), solution_size);
  Eigen::MatrixXd sub_pij_log1m(pij_log.rows(), solution_size);

  // preliminary processing
  /// subset planning units selected in the solution
  {
    std::size_t i = 0;
    for (std::size_t j = 0; j < n_pu; ++j) {
      if (solution[j]) {
        sub_pij_log.col(i) = pij_log.col(j);
        sub_pij_log1m.col(i) = pij_log1m.col(j);
        ++i;
      }
    }
  }

  // main processing
  /// initialize loop variables
  double out = std::numeric_limits<double>::infinity();
  double curr_value_given_state_occurring;
  double curr_probability_of_state_occurring;
  double curr_expected_value_given_state;
  Eigen::MatrixXd curr_state(sub_pij_log.rows(), sub_pij_log.cols());
  /// determine number of states that affect the solution
  mpz_t n;
  mpz_init(n);
  n_states(curr_state.size(), n);
  mpz_add_ui(n, n, 1);
  /// initialize loop iterator
  mpz_t i;
  mpz_init(i);
  mpz_set_ui(i, 1);
  /// iterate over each state
  while (mpz_cmp(i, n) < 0) {
    //// generate the i'th state
    nth_state(i, curr_state);
    //// caculculate the value of the prioritization given the state
    curr_value_given_state_occurring =
      alpha.cwiseProduct(curr_state.rowwise().sum()).array().
        pow(gamma.array()).sum();
    /// if prioritization has a non-zero value then proceed with remaining
    /// calculations for this state
    if (curr_value_given_state_occurring > 1.0e-10) {
      /// calculate probability of the state occurring
      curr_probability_of_state_occurring =
        log_probability_of_state(curr_state, sub_pij_log, sub_pij_log1m);
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
    mpz_add_ui(i, i, 1);
  }

  // clear memory
  mpz_clear(i);
  mpz_clear(n);

  // return result
  return std::exp(out);
}

double expected_value_of_action(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma) {
  // calculate log pij
  Eigen::MatrixXd pij1m = pij;
  log_1m_matrix(pij1m);
  log_matrix(pij);
  // return result
  return expected_value_of_action(solution, pij, pij1m, alpha, gamma);
}

// [[Rcpp::export]]
double rcpp_expected_value_of_action(
  std::vector<bool> solution,
  Eigen::MatrixXd pij,
  Eigen::VectorXd alpha,
  Eigen::VectorXd gamma) {
  // calculate log pij
  Eigen::MatrixXd pij1m = pij;
  log_1m_matrix(pij1m);
  log_matrix(pij);
  // return result
  return expected_value_of_action(solution, pij, pij1m, alpha, gamma);
}
