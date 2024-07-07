#include "rcpp_expected_value_of_action.h"

double expected_value_of_action(
  std::vector<bool> &solution,
  Eigen::MatrixXd &pij,
  Eigen::VectorXi &target) {
  // initialization
  const std::size_t n_pu = pij.cols();
  const std::size_t n_f = pij.rows();
  const std::size_t n_pu_solution = std::accumulate(solution.begin(),
                                                    solution.end(), 0);
  if (static_cast<int>(n_pu_solution) < target.maxCoeff()) {
    Rcpp::stop("prioritization contains fewer planning units than a target");
  }

  // prepare rij matrix containing only solution values
  std::size_t j = 0;
  Eigen::MatrixXd rij(n_f, n_pu_solution);
  for (std::size_t i = 0; i < n_pu; ++i) {
    if (solution[i]) {
      rij.col(j) = pij.col(i);
      ++j;
    }
  }

  // calculate probability of each species meeting the targets
  Eigen::VectorXd spp_prob(n_f);
  Rcpp::NumericVector curr_spp_probs;
  Rcpp::NumericVector curr_density_probs;
  std::vector<int> curr_target_values;
  for (std::size_t i = 0; i < n_f; ++i) {
    curr_spp_probs = Rcpp::wrap(rij.row(i));
    curr_target_values.resize((n_pu_solution - target[i]) + 1);
    std::iota(
      curr_target_values.begin(),
      curr_target_values.end(),
      target[i]
    );
    curr_density_probs = PoissonBinomial::dpb_dc(
      Rcpp::wrap(curr_target_values), curr_spp_probs);
    spp_prob[i] = Rcpp::sum(curr_density_probs);
  }

  // calculate sum of probabilities of species meeting targets
  return spp_prob.sum();
}

double exact_expected_value_of_action(
  Eigen::MatrixXd &pij,
  Rcpp::IntegerVector &target_values) {
  if (static_cast<int>(pij.cols()) < Rcpp::max(target_values)) {
    Rcpp::stop("prioritization contains fewer planning units than a target");
  }
  double out = 0.0;
  const std::size_t n_f = pij.rows();
  Rcpp::NumericVector curr_density_probs;
  Rcpp::NumericVector curr_spp_probs;
  std::vector<int> curr_target_values;
  for (std::size_t i = 0; i < n_f; ++i) {
    curr_spp_probs = wrap(pij.row(i));
    curr_target_values.resize((pij.cols() - target_values[i]) + 1);
    std::iota(
      curr_target_values.begin(),
      curr_target_values.end(),
      target_values[i]
    );
    curr_density_probs = PoissonBinomial::dpb_dc(
      Rcpp::wrap(curr_target_values), curr_spp_probs);
    out += Rcpp::sum(curr_density_probs);
  }
  return out;
}

double approx_expected_value_of_action(
  Eigen::MatrixXd &pij,
  Rcpp::IntegerVector &target_values) {
  if (static_cast<int>(pij.cols()) < Rcpp::max(target_values)) {
    Rcpp::stop("prioritization contains fewer planning units than a target");
  }
  double out = 0.0;
  const std::size_t n_f = pij.rows();
  const std::size_t n_pu = pij.cols();
  Rcpp::NumericVector curr_density_probs;
  std::vector<double> curr_spp_probs;
  std::vector<int> curr_target_values;
  for (std::size_t i = 0; i < n_f; ++i) {
    // find non-zero species probabilities
    // this is because the approximation methods in PoissonBinomial
    // require non-zero values
    curr_spp_probs.reserve(n_pu);
    for (std::size_t j = 0; j < n_pu; ++j) {
      if (pij(i, j) > 0.0) {
        curr_spp_probs.push_back(pij(i, j));
      }
    }
    // if there are no non-zero probabilities, then skip rest of loop
    if (curr_spp_probs.size() == 0) {
      continue;
    }
    // prepare target values
    curr_target_values.resize((curr_spp_probs.size() - target_values[i]) + 1);
    std::iota(
      curr_target_values.begin(),
      curr_target_values.end(),
      target_values[i]
    );
    // calculate probability that target is met
    if (
      (curr_spp_probs.size() == 1) &&
      (target_values[i] == 1.0)
    ) {
      /// if running calculations for only one probability value and
      /// the target is equal to one, then manually calculate the probability
      /// since PoissonBinomial library gets it wrong
      out += curr_spp_probs[0];
    } else {
      /// othweise use the PoissonBinomial library to run the calculations
      curr_density_probs = PoissonBinomial::dpb_na(
        Rcpp::wrap(curr_target_values),
        Rcpp::wrap(curr_spp_probs),
        false
      );
      // add the prob. of the species' target being met to the running total
      out += Rcpp::sum(curr_density_probs);
    }
    // clean up species probability data
    curr_spp_probs.clear();
  }
  return out;
}

double log_proxy_expected_value_of_action(
  Eigen::MatrixXd &log_1m_pij) {
  Eigen::VectorXd out = log_1m_pij.rowwise().sum();
  // note we that we use log_substract(0, ..) because 0 is log(1.0)
  for (auto itr = out.data(); itr != out.data() + out.size(); ++itr)
    *itr = log_subtract(0.0, *itr);
  return log_sum(out);
}

// [[Rcpp::export]]
double rcpp_expected_value_of_action(
  std::vector<bool> solution,
  Eigen::MatrixXd pij,
  Eigen::VectorXi target) {
  // return result
  return expected_value_of_action(
    solution, pij, target);
}
