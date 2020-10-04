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
  std::size_t k, Eigen::MatrixXd &pij, int seed, std::vector<mpz_class> &out) {
  // set seed
  set_seed(seed);
  // sample states
  sample_n_uniform_states_without_replacement(k, pij, out);
  // return void
  return;
}

void sample_n_weighted_states_with_replacement(
  std::size_t k, Eigen::MatrixXd &pij, std::vector<mpz_class> &out) {
  // init
  out.resize(k);
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
  out.resize(k);
  const std::size_t n_v = pij.size();
  Eigen::MatrixXd states(pij.cols(), pij.rows());
  // main
  for (std::size_t i = 0; i < k; ++i) {
    // generate i'th state
    for (std::size_t j = 0; j < n_v; ++j)
      states(j) = Rcpp::rbinom(n_v, 1, 0.5)[0];
    // identify state number
    which_state(states, out[i]);
  }
  // return void
  return;
}

void sample_n_uniform_states_without_replacement(
  std::size_t k, Eigen::MatrixXd &pij, std::vector<mpz_class> &out) {
  // init
  const std::size_t n_v = pij.size();
  const double threshold = std::log(1.0 - 1.0e-5);
  Eigen::MatrixXd states(pij.cols(), pij.rows());
  std::unordered_set<mpz_class> state_set;
  state_set.reserve(k);
  mpz_class state_id;
  std::size_t n_unique_states = 0;
  double sampled_prob = std::numeric_limits<double>::infinity();
  Eigen::MatrixXd log_pij = pij;
  Eigen::MatrixXd log_pij_1m = pij;
  log_matrix(log_pij);
  log_1m_matrix(log_pij_1m);
  double state_prob;
  // generate states
  while(n_unique_states < k) {
    // generate new state
    std::transform(pij.data(), pij.data() + pij.size(), states.data(),
                   [=](double p){ return R::rbinom(1, 0.5); });
    // identify state index
    state_id = 0;
    which_state(states, state_id);
    // add the state to the set
    if (state_set.insert(state_id).second) {
      // if the newly sampled state hasn't been generated before,
      // then increment number of unique states
      ++n_unique_states;
      // calculate probability of the state occurring
      state_prob = log_probability_of_outcome(states, log_pij, log_pij_1m);
      // add probability of the state to the running total of all sampled states
      if (std::isinf(sampled_prob)) {
        sampled_prob = state_prob;
      } else {
        sampled_prob = log_sum(sampled_prob, state_prob);
      }
      // stop looking for new states if the remaining stats have a very low
      // probability of occuring
      if (sampled_prob >= threshold)
        break;
    }
  }
  // extract states
  out.resize(state_set.size());
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
  const double threshold = std::log(1.0 - 1.0e-5);
  Eigen::MatrixXd states(pij.cols(), pij.rows());
  std::unordered_set<mpz_class> state_set;
  state_set.reserve(k);
  mpz_class state_id;
  std::size_t n_unique_states = 0;
  double sampled_prob = std::numeric_limits<double>::infinity();
  Eigen::MatrixXd log_pij = pij;
  Eigen::MatrixXd log_pij_1m = pij;
  log_matrix(log_pij);
  log_1m_matrix(log_pij_1m);
  double state_prob;
  // generate states
  while(n_unique_states < k) {
    // generate new state
    std::transform(pij.data(), pij.data() + pij.size(), states.data(),
                   [=](double p){ return R::rbinom(1, p); });
    // identify state index
    state_id = 0;
    which_state(states, state_id);
    // add the state to the set
    if (state_set.insert(state_id).second) {
      // if the newly sampled state hasn't been generated before,
      // then increment number of unique states
      ++n_unique_states;
      // calculate probability of the state occurring
      state_prob = log_probability_of_outcome(states, log_pij, log_pij_1m);
      // add probability of the state to the running total of all sampled states
      if (std::isinf(sampled_prob)) {
        sampled_prob = state_prob;
      } else {
        sampled_prob = log_sum(sampled_prob, state_prob);
      }
      // stop looking for new states if the remaining stats have a very low
      // probability of occuring
      if (sampled_prob >= threshold)
        break;
    }
  }
  // extract states
  out.resize(state_set.size());
  std::vector<mpz_class>::iterator itr = out.begin();
  for (std::unordered_set<mpz_class>::iterator itr2 = state_set.begin();
       itr2 != state_set.end(); ++itr2, ++itr)
    *itr = *itr2;
  // return void
  return;
}

// [[Rcpp::export]]
std::vector<std::size_t> rcpp_sample_n_weighted_states_with_replacement(
  std::size_t k, Eigen::MatrixXd &pij, int seed) {
  // init
  std::vector<mpz_class> s;
  set_seed(seed);
  // generate states
  sample_n_weighted_states_with_replacement(k, pij, s);
  // extract values
  const std::size_t n = s.size();
  std::vector<std::size_t> o(n);
  for (std::size_t i = 0; i < n; ++i)
    o[i] = s[i].get_ui();
  // return result
  return o;
}

// [[Rcpp::export]]
std::vector<std::size_t> rcpp_sample_n_uniform_states_with_replacement(
  std::size_t k, Eigen::MatrixXd &pij, int seed) {
  // init
  std::vector<mpz_class> s;
  set_seed(seed);
  // generate states
  sample_n_uniform_states_with_replacement(k, pij, s);
  // extract values
  const std::size_t n = s.size();
  std::vector<std::size_t> o(n);
  for (std::size_t i = 0; i < n; ++i)
    o[i] = s[i].get_ui();
  // return result
  return o;
}

// [[Rcpp::export]]
std::vector<std::size_t> rcpp_sample_n_uniform_states_without_replacement(
  std::size_t k, Eigen::MatrixXd &pij, int seed) {
  // init
  std::vector<mpz_class> s;
  set_seed(seed);
  // generate states
  sample_n_uniform_states_without_replacement(k, pij, s);
  // extract values
  const std::size_t n = s.size();
  std::vector<std::size_t> o(n);
  for (std::size_t i = 0; i < n; ++i)
    o[i] = s[i].get_ui();
  // return result
  return o;
}

// [[Rcpp::export]]
std::vector<std::size_t> rcpp_sample_n_weighted_states_without_replacement(
  std::size_t k, Eigen::MatrixXd &pij, int seed) {
  // init
  std::vector<mpz_class> s;
  set_seed(seed);
  // generate states
  sample_n_weighted_states_without_replacement(k, pij, s);
  // extract values
  const std::size_t n = s.size();
  std::vector<std::size_t> o(n);
  for (std::size_t i = 0; i < n; ++i)
    o[i] = s[i].get_ui();
  // return result
  return o;
}
