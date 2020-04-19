#include "rcpp_approx_expected_value_of_action.h"

void approx_expected_value_of_action_values(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij_log,
  Eigen::MatrixXd &pij_log1m,
  Eigen::VectorXd &preweight,
  Eigen::VectorXd &postweight,
  Eigen::VectorXd &target,
  std::vector<mpz_class> &states,
  std::vector<double> &value_given_state_occurring,
  std::vector<double> &prob_of_state_occurring,
  std::vector<double> &extra_prob_of_state_occurring,
  double extra_prob = 0) // used if state probabilities are conditional on
                        // another event (e.g. survey outcomes)
  {
  // initialization
  const std::size_t n_pu = pij_log.cols();
  const double total = static_cast<double>(n_pu);
  const std::size_t n_approx_states = states.size();

  // main processing
  /// initialize loop variables
  double v, p;
  Eigen::MatrixXd curr_state(pij_log.rows(), pij_log.cols());

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
    v = conservation_benefit_state(
      curr_state, preweight, postweight, target, total);
    // store data
    if (v > 1.0e-10) {
      value_given_state_occurring.push_back(v);
      prob_of_state_occurring.push_back(p + extra_prob);
    } else {
      extra_prob_of_state_occurring.push_back(p + extra_prob);
    }
  }

  // return void
  return;
}

double approx_expected_value_of_action(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij_log,
  Eigen::MatrixXd &pij_log1m,
  Eigen::VectorXd &preweight,
  Eigen::VectorXd &postweight,
  Eigen::VectorXd &target,
  std::vector<mpz_class> &states) {

  // initialization
  const std::size_t n_approx_states = states.size();
  std::vector<double> value_given_state_occurring;
  std::vector<double> prob_of_state_occurring;
  std::vector<double> extra_prob_of_state_occurring;
  value_given_state_occurring.reserve(n_approx_states);
  prob_of_state_occurring.reserve(n_approx_states);
  extra_prob_of_state_occurring.reserve(n_approx_states);

  // main processing
  approx_expected_value_of_action_values(
    solution, pij_log, pij_log1m, preweight, postweight, target, states,
    value_given_state_occurring, prob_of_state_occurring,
    extra_prob_of_state_occurring);

  // check that at least one state had a non-zero value
  const std::size_t k = value_given_state_occurring.size();
  const std::size_t k2 = extra_prob_of_state_occurring.size();
  assert_gt_value(k, (std::size_t) 0,
    "all states have zero value, try increasing argument to n_approx_states_per_replicate");

 // create Eigen maps of data
  Eigen::VectorXd value_given_state_occurring2 =
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
      value_given_state_occurring.data(), k);
  Eigen::VectorXd prob_of_state_occurring2 =
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
      prob_of_state_occurring.data(), k);
  Eigen::VectorXd extra_prob_of_state_occurring2 =
    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
      extra_prob_of_state_occurring.data(), k2);

  // calculate total probaiblity
  double total_prob = log_sum(prob_of_state_occurring2);
  if (k2 > 0)
    total_prob = log_sum(total_prob, log_sum(extra_prob_of_state_occurring2));

  // caclulate corrected values
  prob_of_state_occurring2.array() -= total_prob;
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
    solution, pij, pij_log1m, preweight, postweight, target, states2);

  // return result
  return out;
}
