#include "package.h"
#include "functions.h"
#include "rcpp_states.h"
#include "rcpp_probability.h"
#include "rcpp_prioritization.h"
#include "rcpp_approximate_expected_value_of_management_action.h"

// [[Rcpp::export]]
double rcpp_appproximate_expected_value_of_management_decision_given_current_information(
  Eigen::MatrixXd &pij,
  Eigen::VectorXd &pu_costs,
  Eigen::VectorXd &pu_locked_in,
  Eigen::VectorXd &alpha,
  Eigen::VectorXd &gamma,
  std::size_t n_approx_obj_fun_points,
  double budget,
  double gap,
  const std::size_t n_approx_states) {

  // find optimal management action using prior data
  std::vector<bool> solution(pij.cols());
  Prioritization p(pij.cols(), pij.rows(), pu_costs, pu_locked_in,
                   alpha, gamma, n_approx_obj_fun_points, budget, gap);
  p.add_rij_data(pij);
  p.solve();
  p.get_solution(solution);

  // calculate log prior probabilities
  pij.array() = pij.array().log();

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

  // calculate expected value of management action
  double out = approximate_expected_value_of_management_action(
    solution, pij, alpha, gamma, states);

  // clear memory
  gmp_randclear(rng_state);
  mpz_clear(n_states_total);
  for (std::size_t i = 0; i < n_approx_states; ++i)
    mpz_clear(states[i]);

  // return result
  return out;
}
