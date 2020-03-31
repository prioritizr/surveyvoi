#include "rcpp_sample_states.h"

void sample_n_weighted_states_with_replacement(
  std::size_t k, Eigen::MatrixXd &pij, std::vector<mpz_class> &out) {
  // init
  const std::size_t n_v = pij.size();
  Eigen::MatrixXd states(pij.cols(), pij.rows());
  // main
  for (std::size_t i = 0; i < k; ++i) {
    // generate i'th state
    for (std::size_t j = 0; j < n_v; ++j)
      states(j) = Rcpp::rbinom(n_v, 1, pij(j))[0];
    // identify state number
    which_state(states, out[i]);
  }
  // return void
  return;
}

void sample_n_uniform_states_with_replacement(
  std::size_t k, Eigen::MatrixXd &pij, std::vector<mpz_class> &out) {
  // init
  std::size_t seed = static_cast<std::size_t>(Rcpp::sample(1e+6, 1, true)[0]);
  mpz_t n, tmp;
  mpz_init(n);
  mpz_init(tmp);
  n_states(pij.size(), n);
  mpz_add_ui(n, n, 1);
  gmp_randstate_t rng;
  gmp_randinit_default(rng);
  gmp_randseed_ui(rng, seed);
  // main
  for (std::size_t i = 0; i < k; ++i) {
    mpz_urandomm(tmp, rng, n);
    out[i] = mpz_class(tmp);
  }
  // clean up
  mpz_clear(n);
  mpz_clear(tmp);
  gmp_randclear(rng);
  // return void
  return;
}

// based on https://stackoverflow.com/a/25179870/3483791
void sample_n_uniform_states_without_replacement(
  std::size_t k, Eigen::MatrixXd &pij, std::vector<mpz_class> &out) {

  // initialize random number generator
  std::size_t seed = static_cast<std::size_t>(Rcpp::sample(1e+6, 1, true)[0]);
  gmp_randstate_t rng;
  gmp_randinit_default(rng);
  gmp_randseed_ui(rng, seed);

  // initialize uint gmp numbers
  mpz_t t_mpz, N_mpz;
  mpz_init_set_ui(t_mpz, 0);
  mpz_init(N_mpz);
  n_states(pij.size(), N_mpz);
  mpz_add_ui(N_mpz, N_mpz, 1);

  // initialize floating point gmp numbers
  mpf_t N_mpf, t_mpf, u_mpf, s_mpf, z_mpf;
  mpf_init(u_mpf);
  mpf_init(s_mpf);
  mpf_init(z_mpf);
  mpf_init_set_ui(t_mpf, 0);
  mpf_init(N_mpf);
  mpf_set_z(N_mpf, N_mpz);

  // initialize uint numbers
  std::size_t m = 0;

  // main
  while (m < k) {
    /// generate random float between zero and one
    mpf_urandomb(u_mpf, rng, 23);
    /// randomly determine if the t should be appended to the result
    mpf_sub(s_mpf, N_mpf, t_mpf);
    mpf_mul(z_mpf, s_mpf, u_mpf);
    if (mpf_cmp_ui(z_mpf, k - m) >= 0) {
      mpf_add_ui(t_mpf, t_mpf, 1);
      mpz_add_ui(t_mpz, t_mpz, 1);
    } else {
      out[m] = mpz_class(t_mpz);
      ++m;
      mpf_add_ui(t_mpf, t_mpf, 1);
      mpz_add_ui(t_mpz, t_mpz, 1);
    }
  }
  // clean up
  mpz_clear(t_mpz);
  mpz_clear(N_mpz);
  mpf_clear(u_mpf);
  mpf_clear(N_mpf);
  mpf_clear(t_mpf);
  mpf_clear(s_mpf);
  mpf_clear(z_mpf);
  gmp_randclear(rng);
  // return void
  return;
}

namespace std {
template <> struct hash<mpz_class>
{
  size_t operator()(const mpz_class& x) const
  {
    return x.get_mpz_t()[0]._mp_size != 0 ?
      static_cast<size_t>(x.get_mpz_t()[0]._mp_d[0]) : 0;
  }
};
}

void sample_n_weighted_states_without_replacement(
  std::size_t k, Eigen::MatrixXd &pij, std::vector<mpz_class> &out) {
  // init
  const std::size_t n_v = pij.size();
  Eigen::MatrixXd states(pij.cols(), pij.rows());
  std::unordered_set<mpz_class> state_set;
  state_set.reserve(k);
  mpz_class state_id;
  std::size_t n_unique_states = 0;

  // generate states
  while(n_unique_states < k) {
    // generate new state
    for (std::size_t j = 0; j < n_v; ++j)
      states(j) = Rcpp::rbinom(n_v, 1, pij(j))[0];
    // identify state index
    state_id = 0;
    which_state(states, state_id);
    // add state to set and increment number of unique states
    if (state_set.insert(state_id).second) {
      ++n_unique_states;
    }
  }

  // extract states
  std::vector<mpz_class>::iterator itr = out.begin();
  for (std::unordered_set<mpz_class>::iterator itr2 = state_set.begin();
       itr2 != state_set.end(); ++itr2, ++itr)
    *itr = *itr2;

  // return void
  return;
}

// [[Rcpp::export]]
std::vector<std::size_t> rcpp_sample_n_weighted_states_with_replacement(
  std::size_t k, Eigen::MatrixXd &pij) {
  // init
  std::vector<mpz_class> s(k);
  std::vector<std::size_t> o(k);
  // generate states
  sample_n_weighted_states_with_replacement(k, pij, s);
  // extract values
  for (std::size_t i = 0; i < k; ++i)
    o[i] = s[i].get_ui();
  // return result
  return o;
}

// [[Rcpp::export]]
std::vector<std::size_t> rcpp_sample_n_uniform_states_with_replacement(
  std::size_t k, Eigen::MatrixXd &pij) {
  // init
  std::vector<mpz_class> s(k);
  std::vector<std::size_t> o(k);
  // generate states
  sample_n_uniform_states_with_replacement(k, pij, s);
  // extract values
  for (std::size_t i = 0; i < k; ++i)
    o[i] = s[i].get_ui();
  // return result
  return o;
}

// [[Rcpp::export]]
std::vector<std::size_t> rcpp_sample_n_uniform_states_without_replacement(
  std::size_t k, Eigen::MatrixXd &pij) {
  // init
  std::vector<mpz_class> s(k);
  std::vector<std::size_t> o(k);
  // generate states
  sample_n_uniform_states_without_replacement(k, pij, s);
  // extract values
  for (std::size_t i = 0; i < k; ++i)
    o[i] = s[i].get_ui();
  // return result
  return o;
}

// [[Rcpp::export]]
std::vector<std::size_t> rcpp_sample_n_weighted_states_without_replacement(
  std::size_t k, Eigen::MatrixXd &pij) {
  // init
  std::vector<mpz_class> s(k);
  std::vector<std::size_t> o(k);
  // generate states
  sample_n_weighted_states_without_replacement(k, pij, s);
  // extract values
  for (std::size_t i = 0; i < k; ++i)
    o[i] = s[i].get_ui();
  // return result
  return o;
}
