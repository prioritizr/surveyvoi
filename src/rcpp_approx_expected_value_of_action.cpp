#include "rcpp_approx_expected_value_of_action.h"

double approx_expected_value_of_action(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij_log,
  Eigen::MatrixXd &pij_log1m,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::vector<mpz_t> &states) {
  // initialization
  const std::size_t n_pu = pij_log.cols();
  const std::size_t solution_size =
    std::accumulate(solution.begin(), solution.end(), 0);
  Eigen::MatrixXd sub_pij_log(pij_log.rows(), solution_size);
  Eigen::MatrixXd sub_pij_log1m(pij_log1m.rows(), solution_size);

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
  /// determine number of states that affect the solution
  const std::size_t n_approx_states = states.size();
  /// initialize loop variables
  double v;
  std::size_t k = 0;
  Eigen::MatrixXd curr_state(sub_pij_log.rows(), sub_pij_log.cols());
  std::vector<double> value_given_state_occurring;
  std::vector<double> prob_of_state_occurring;
  value_given_state_occurring.reserve(n_approx_states);
  prob_of_state_occurring.reserve(n_approx_states);

  /// iterate over each state
  for (std::size_t i = 0; i < n_approx_states; ++i) {
    //// generate the i'th state
    nth_state(states[i], curr_state);
    //// calculate the value of the prioritization given the state
    v = alpha.cwiseProduct(curr_state.rowwise().sum()).array().
        pow(gamma.array()).sum();
    if (v > 1.0e-10) {
      /// store value if the priorititization has a non-zero benefit
      value_given_state_occurring.push_back(v);
      /// calculate probability of the state occurring
      prob_of_state_occurring.push_back(
        log_probability_of_state(curr_state, sub_pij_log, sub_pij_log1m));
      /// increment counter
      ++k;
    }
  }

  // check that at least one state had a non-zero value
  assert_gt_value(k, (std::size_t) 0,
    "all states have zero value, try increasing argument to n_approx_states");

  // create Eigen maps of data
  Eigen::VectorXd value_given_state_occurring2 =
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
      value_given_state_occurring.data(), k);
  Eigen::VectorXd prob_of_state_occurring2 =
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
      prob_of_state_occurring.data(), k);

  // rescale probabilities
  prob_of_state_occurring2.array() -=
    log_sum(prob_of_state_occurring2);

  // calculate values weighted by probabilities
  prob_of_state_occurring2.array() +=
    value_given_state_occurring2.array().log();

  // return result
  return std::exp(log_sum(prob_of_state_occurring2));
}

// [[Rcpp::export]]
double rcpp_approx_expected_value_of_action(
  std::vector<bool> solution,
  Eigen::MatrixXd pij,
  Eigen::VectorXd alpha,
  Eigen::VectorXd gamma,
  std::vector<std::size_t> states) {
  // calculate log pij
  Eigen::MatrixXd pij_log1m = pij;
  pij_log1m.array() = (1.0 - pij_log1m.array()).array().log();
  pij.array() = pij.array().log();

  // convert state indices from std::size_t to mpz_t
  const std::size_t n = states.size();
  std::vector<mpz_t> states2(n);
  for (std::size_t i = 0; i < n; ++i) {
    mpz_init(states2[i]);
    mpz_set_ui(states2[i], states[i]);
  }

  // calculate result
  double out = approx_expected_value_of_action(
    solution, pij, pij_log1m, alpha, gamma, states2);

  // clear memory
  for (std::size_t i = 0; i < n; ++i)
    mpz_clear(states2[i]);

  // return result
  return out;
}
