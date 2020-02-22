#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"

// [[Rcpp::export]]
double rcpp_approximate_expected_value_of_management_decision_given_perfect_information(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap,
  const std::size_t n_approx_states) {

  // initialization
  /// create log version of probabilities
  Eigen::MatrixXd pij_log(pij.cols(), pij.rows());
  pij_log = pij.array().log();

  /// initialize prioritization
  std::vector<bool> solution(pij.cols());
  Prioritization p(pij.cols(), pij.rows(), pu_costs, pu_locked_in,
                   alpha, gamma, n_approx_obj_fun_points, budget, gap);

  // calculate total number of states
  mpz_t n_states_total;
  mpz_init(n_states_total);
  n_states(pij.size(), n_states_total);

  // generate random seed for generating states
  unsigned long int max_uint = std::numeric_limits<unsigned long int>::max();
  unsigned long int rng_seed = Rcpp::sample(max_uint, 1)[0];

  // generate states
  gmp_randstate_t rng_state;
  gmp_randinit_default(rng_state);
  gmp_randseed_ui(rng_state, rng_seed);
  std::vector<mpz_t> states(n_approx_states);
  for (std::size_t i = 0; i < n_approx_states; ++i) {
    mpz_init(states[i]);
    mpz_urandomm(states[i], rng_state, n_states_total);
  }

  /// initialize loop variables
  Eigen::VectorXd curr_value_given_state_occurring(n_approx_states);
  Eigen::VectorXd curr_probability_of_state_occurring(n_approx_states);
  Eigen::MatrixXd curr_state(pij.rows(), pij.cols());

  // main processing
  for (std::size_t i = 0; i < n_approx_states; ++i) {
    /// generate the i'th state
    nth_state(states[i], curr_state);
    /// generate solution for state
    p.add_rij_data(curr_state);
    p.solve();
    p.get_solution(solution);
    /// calculate the value of the prioritization given the state
    curr_value_given_state_occurring[i] =
      std::log(alpha.cwiseProduct(curr_state.rowwise().sum()).array().
        pow(gamma.array()).sum());
    /// calculate probability of the state occurring
    curr_probability_of_state_occurring[i] =
      log_probability_of_state(curr_state, pij_log);
  }

  // rescale probabilities
  curr_probability_of_state_occurring.array() /=
    curr_probability_of_state_occurring.sum();

  // calculate values weighted by probabilities
  curr_probability_of_state_occurring.array() +=
    curr_value_given_state_occurring.array();

  // clear memory
  gmp_randclear(rng_state);
  mpz_clear(n_states_total);
  for (std::size_t i = 0; i < n_approx_states; ++i)
    mpz_clear(states[i]);

  // return result
  return std::exp(log_sum(curr_probability_of_state_occurring));
}
