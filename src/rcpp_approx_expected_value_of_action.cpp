#include "rcpp_approx_expected_value_of_action.h"

double approx_expected_value_of_action(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij_log,
  Eigen::MatrixXd &pij_log1m,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::vector<mpz_class> &states) {
  // initialization
  const std::size_t n_pu = pij_log.cols();
  const std::size_t n_approx_states = states.size();

  // main processing
  /// initialize loop variables
  double v, p;
  std::size_t k = 0;
  Eigen::MatrixXd curr_state(pij_log.rows(), pij_log.cols());
  Eigen::MatrixXd curr_rij(pij_log.rows(), pij_log.cols());
  std::vector<double> value_given_state_occurring;
  std::vector<double> prob_of_state_occurring;
  std::vector<double> all_prob_of_state_occurring;
  value_given_state_occurring.reserve(n_approx_states);
  prob_of_state_occurring.reserve(n_approx_states);
  all_prob_of_state_occurring.reserve(n_approx_states);

  /// iterate over each state
  for (std::size_t i = 0; i < n_approx_states; ++i) {
    //// generate the i'th state
    nth_state(states[i], curr_state);
    /// create matrix only containing feature data for selected planning units
    curr_rij = curr_state;
    for (std::size_t j = 0; j < n_pu; ++j)
      curr_rij.col(j) *= solution[j];
    //// calculate the value of the prioritization given the state
    v = alpha.cwiseProduct(curr_rij.rowwise().sum()).array().
        pow(gamma.array()).sum();
    /// calculate probability of the state occurring
    p = log_probability_of_state(curr_state, pij_log, pij_log1m);
    /// store probability of state occurring
    all_prob_of_state_occurring.push_back(p);
    /// store values if non-zero benefit
    if (v > 1.0e-10) {
      value_given_state_occurring.push_back(v);
      prob_of_state_occurring.push_back(p);
      ++k;
    }
  }

  // check that at least one state had a non-zero value
  assert_gt_value(k, (std::size_t) 0,
    "all states have zero value, try increasing argument to n_approx_states_per_replicate");

  // create Eigen maps of data
  Eigen::VectorXd value_given_state_occurring2 =
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
      value_given_state_occurring.data(), k);
  Eigen::VectorXd prob_of_state_occurring2 =
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
      prob_of_state_occurring.data(), k);
  Eigen::VectorXd all_prob_of_state_occurring2 =
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
      all_prob_of_state_occurring.data(), n_approx_states);

  // rescale probabilities
  prob_of_state_occurring2.array() -= log_sum(all_prob_of_state_occurring2);

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
  log_1m_matrix(pij_log1m);
  log_matrix(pij);

  // convert state indices from std::size_t to mpz_class
  const std::size_t n = states.size();
  std::vector<mpz_class> states2(n);
  for (std::size_t i = 0; i < n; ++i)
    states2[i] = states[i];

  // calculate result
  double out = approx_expected_value_of_action(
    solution, pij, pij_log1m, alpha, gamma, states2);

  // return result
  return out;
}
