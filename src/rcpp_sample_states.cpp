#include "rcpp_sample_states.h"

namespace std {
template <> struct hash<mpz_class>
{
  size_t operator()(const mpz_class& x) const
  {
    return std::hash<std::string>{}(x.get_str());
  }
};
}

void sample_n_states(
  std::size_t k, Eigen::MatrixXd &pij, std::string &method,
  std::vector<mpz_class> &out) {
  // sample states according to specififed method
  if (method == "weighted_without_replacement") {
    sample_n_weighted_states_without_replacement(k, pij, out);
  } else if (method == "weighted_with_replacement") {
    sample_n_weighted_states_with_replacement(k, pij, out);
  } else if (method == "uniform_with_replacement") {
    sample_n_uniform_states_with_replacement(k, pij, out);
  } else if (method == "uniform_without_replacement") {
    sample_n_uniform_states_without_replacement(k, pij, out);
  } else {
    Rcpp::stop("method not recognized");
  }
  // return void
  return;
}

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

void sample_n_uniform_states_without_replacement(
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
    std::transform(pij.data(), pij.data() + pij.size(), states.data(),
                   [=](double p){ return R::rbinom(1, 0.5); });
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
    std::transform(pij.data(), pij.data() + pij.size(), states.data(),
                   [=](double p){ return R::rbinom(1, p); });
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
