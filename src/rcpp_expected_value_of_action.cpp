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
  if (static_cast<int>(n_pu_solution) < target.minCoeff()) {
    Rcpp::warning(
      "prioritization contains fewer planning units than minimum target"
    );
    return 0.0;
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
    if (curr_spp_probs.size() < target[i]) {
      spp_prob[i] = 0.0;
    } else {
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
  }

  // calculate sum of probabilities of species meeting targets
  return spp_prob.sum();
}

// this function is used for internal debugging when we need to calculate
// the expected value of actions using methods that produce correct
// (not approximation) results

// # nocov start
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
// # nocov end

double approx_expected_value_of_action(
  std::vector<std::vector<double>> &pij,
  Rcpp::IntegerVector &curr_n,
  Rcpp::IntegerVector &target_values,
  std::vector<std::size_t> &extra_n
) {
  // initialze values
  double out = 0.0;
  std::size_t curr_spp_n;
  const std::size_t n_f = pij.size();
  Rcpp::NumericVector curr_density_probs;
  std::vector<int> curr_target_values;
  std::vector<double> curr_spp_probs;
  // main calculations
  for (std::size_t i = 0; i < n_f; ++i) {
    // determine number of planning units with non-zero probabilities
    // of occupancy for species
    curr_spp_n = curr_n[i] + extra_n[i];
    // calculate probability that species target is met
    if ((curr_spp_n == 1) && (target_values[0] == 1)) {
      // manually calculate if only one selected planning unit has non-zero
      // probability of species occurrence and target value of 1
      // this is due to a bug in PoissonBinomial library
      out += pij[i][0];
    } else if (
      (curr_spp_n > 1) &&
      !(static_cast<std::size_t>(target_values[i]) > curr_spp_n)
    ) {
      // otherwise, use PoissonBinomial library to run calculations
      /// prepare probability data
      curr_spp_probs.resize(curr_spp_n);
      std::copy(
        pij[i].cbegin(),
        pij[i].cbegin() + curr_spp_n,
        curr_spp_probs.begin()
      );

      /// prepare target data
      curr_target_values.resize((curr_spp_n - target_values[i]) + 1);
      std::iota(
        curr_target_values.begin(),
        curr_target_values.end(),
        target_values[i]
      );

      /// otherwise, use the PoissonBinomial library to run the calculations
      out += Rcpp::sum(
        PoissonBinomial::dpb_na(
          Rcpp::wrap(curr_target_values),
          Rcpp::wrap(curr_spp_probs),
          false
        )
      );
    }
  }
  // return result
  return out;
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
