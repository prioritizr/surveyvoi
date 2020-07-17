#include "rcpp_approx_expected_value_of_action.h"

double approx_expected_value_of_action(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij,
  Eigen::MatrixXd &pij_log,
  Eigen::MatrixXd &pij_log1m,
  Eigen::VectorXd &preweight,
  Eigen::VectorXd &postweight,
  Eigen::VectorXd &target,
  std::size_t n_approx_states) {

  // initialization
  const std::size_t n_pu = pij_log.cols();
  const double total = static_cast<double>(n_pu);
  const std::size_t solution_size =
    std::accumulate(solution.begin(), solution.end(), 0);

  // preliminary processing
  /// generate states
  std::vector<mpz_class> states(n_approx_states);
  sample_n_weighted_states_without_replacement(n_approx_states, pij, states);

  // Rcout << "states = ";
  // for (std::size_t i = 0; i < n_approx_states; ++i)
  //   Rcout << states[i].get_str() << ", ";
  // Rcout << std::endl;

  // main processing
  /// initialize loop variables
  double v, p;
  std::size_t k = 0;
  Eigen::MatrixXd curr_state(pij_log.rows(), pij_log.cols());
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
    /// calculate probability of the state occurring
    p = log_probability_of_state(curr_state, pij_log, pij_log1m);
    /// create matrix only containing feature data for selected planning units
    for (std::size_t j = 0; j < n_pu; ++j)
      curr_state.col(j) *= solution[j];
    //// calculate the value of the prioritization given the state
    v = conservation_value_state(
      curr_state, preweight, postweight, target, total);
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

  // print(wrap(value_given_state_occurring2));
  // print(wrap(prob_of_state_occurring2));

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
  Eigen::VectorXd preweight,
  Eigen::VectorXd postweight,
  Eigen::VectorXd target,
  std::size_t n_approx_states) {
  // calculate log matrices
  Eigen::MatrixXd pij_log = pij;
  Eigen::MatrixXd pij_log1m = pij;
  log_matrix(pij_log);
  log_1m_matrix(pij_log1m);

  // return result
  return approx_expected_value_of_action(
    solution, pij, pij_log, pij_log1m, preweight, postweight, target,
    n_approx_states);
}
