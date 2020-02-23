#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_approx_expected_value_of_action.h"

double approx_expected_value_of_decision_given_perfect_info(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap,
  std::vector<mpz_t> &states) {

  // initialization
  // constant variables
  const std::size_t n_pu = pij.cols();
  const std::size_t n_approx_states = states.size();

  /// create log version of probabilities
  Eigen::MatrixXd pij_log(pij.cols(), pij.rows());
  pij_log = pij.array().log();

  /// initialize prioritization
  std::vector<bool> solution(n_pu);
  Prioritization p(pij.cols(), pij.rows(), pu_costs, pu_locked_in,
                   alpha, gamma, n_approx_obj_fun_points, budget, gap);

  /// initialize loop variables
  double curr_value_given_state_occurring;
  std::vector<double> value_given_state_occurring;
  std::vector<double> prob_of_state_occurring;
  value_given_state_occurring.reserve(n_approx_states);
  prob_of_state_occurring.reserve(n_approx_states);
  Eigen::MatrixXd curr_state(pij.rows(), pij.cols());
  Eigen::MatrixXd curr_rij(pij.rows(), pij.cols());
  std::size_t k = 0;

  // main processing
  for (std::size_t i = 0; i < n_approx_states; ++i) {
    /// generate the i'th state
    nth_state(states[i], curr_state);
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
      /// store value given state
      value_given_state_occurring.push_back(curr_value_given_state_occurring);
      /// store probability of state occurring
      prob_of_state_occurring.push_back(
        log_probability_of_state(curr_state, pij_log));
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
  prob_of_state_occurring2.array() /=
    prob_of_state_occurring2.sum();

  // calculate values weighted by probabilities
  prob_of_state_occurring2.array() +=
    value_given_state_occurring2.array().log();

  // return result
  return std::exp(log_sum(prob_of_state_occurring2));
}

// [[Rcpp::export]]
double rcpp_approx_expected_value_of_decision_given_perfect_info_n_states(
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
  double out = approx_expected_value_of_decision_given_perfect_info(
    pij, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points, budget,
    gap, states);

  // clear memory
  for (std::size_t i = 0; i < n_approx_states; ++i)
    mpz_clear(states[i]);

  // return result
  return out;
}

// [[Rcpp::export]]
double rcpp_approx_expected_value_of_decision_given_perfect_info_fixed_states(
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
  double out = approx_expected_value_of_decision_given_perfect_info(
    pij, pu_costs, pu_locked_in, alpha, gamma, n_approx_obj_fun_points, budget,
    gap, states2);

  // clear memory
  for (std::size_t i = 0; i < n; ++i)
    mpz_clear(states2[i]);

  // return result
  return out;
}
