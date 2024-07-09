#include "rcpp_heuristic_prioritization.h"

void greedy_heuristic_prioritization(
  Eigen::MatrixXd& rij, Eigen::VectorXd& pu_costs,
  Eigen::VectorXd& pu_locked_in, Eigen::VectorXd& pu_locked_out,
  Eigen::VectorXi& target, double budget, std::vector<bool>& solution) {
  // Initialization
  // declare constants
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = rij.rows();
  const std::size_t max_target = static_cast<std::size_t>(target.maxCoeff());
  // reset solution
  std::fill(solution.begin(), solution.end(), false);
  for (std::size_t i = 0; i < n_pu; ++i)
    if ((pu_locked_in[i] > 0.5) ||
        ((pu_costs[i] < 1.0e-15) && (pu_locked_out[i] < 0.5)))
      solution[i] = true;
  // calculate solution cost
  double solution_cost = 0.0;
  for (std::size_t i = 0; i < n_pu; ++i)
    solution_cost += static_cast<double>(solution[i]) * pu_costs[i];
  // declare loop variables
  /// variables to keep track of planning units in the solution that
  /// will be considered for step-wise addition from the algorithm
  std::vector<bool> solution_rem_pu(n_pu, false);
  for (std::size_t i = 0; i < n_pu; ++i)
    if ((pu_locked_out[i] < 0.5) && (!solution[i]))
      solution_rem_pu[i] = true;
  double curr_min_feasible_pu_cost;
  double curr_ce;
  double new_ce;
  double obj;
  double curr_obj;
  double prev_obj;
  Rcpp::IntegerVector real_target = Rcpp::wrap(target);
  Rcpp::IntegerVector curr_target;
  std::size_t n_additional_units_needed;
  std::size_t curr_idx;
  std::size_t n_pu_selected = std::accumulate(
    solution.begin(), solution.end(), 0);
  std::size_t n_pu_remaining = std::accumulate(
    solution_rem_pu.begin(), solution_rem_pu.end(), 0);
  std::vector<std::size_t> extra_n(n_f);
  std::vector<std::size_t> zeros(n_f, 0);
  // variable to keep track of costs in unselected planning units
  Eigen::VectorXd costs_rem_pu = pu_costs;
  Eigen::VectorXd costs_rem_pu_sorted;
  for (std::size_t i = 0; i < n_pu; ++i)
    if (solution[i] || (pu_locked_out[i] > 0.5))
      costs_rem_pu[i] = std::numeric_limits<double>::max();
  std::size_t min_cost_pu_idx =
    std::min_element(costs_rem_pu.data(),
                     costs_rem_pu.data() + costs_rem_pu.size()) -
    costs_rem_pu.data();

  // initialize rij data for optimization
  std::vector<std::vector<double>> curr_rij(n_f);
  Rcpp::IntegerVector curr_n(n_f, 0);
  for (std::size_t i = 0; i < n_f; ++i) {
    curr_rij[i].resize(n_pu);
  }
  for (std::size_t i = 0; i < n_pu; ++i) {
    if (solution[i]) {
      for (std::size_t j = 0; j < n_f; ++j) {
        if (rij(j, i) > 1.0e-10) {
          curr_rij[j][curr_n[j]] = rij(j, i);
          ++curr_n[j];
        }
      }
    }
  }

  // calculate objective value for starting solution
  if (n_pu_selected > 0) {
    curr_target = Rcpp::pmin(real_target, curr_n);
    prev_obj = approx_expected_value_of_action(
      curr_rij, curr_n, curr_target, zeros
    );
  } else {
    prev_obj = 0.0;
    curr_target = Rcpp::pmin(real_target, 1);
  }

  // Main processing
  while(
    ((solution_cost + pu_costs[min_cost_pu_idx]) < budget) &&
    (n_pu_remaining > 0)
  ) {

    // here we cap targets to be the minimum of the "real"
    // target and the current number of planning units selected.
    // this is important so that we can account for incremental progress
    // towards meeting the targets when the solution contains fewer
    // than the number of selected planning units.
    // note we must relcalculate the previous object ive value if changing
    // the targets, so that objective values from this iteration
    // can be compared with current solution
    if (
      (n_pu_selected > 0) &&
      Rcpp::is_true(Rcpp::any(real_target != curr_n))
    ) {
      curr_target = Rcpp::pmin(real_target, curr_n);
      prev_obj = approx_expected_value_of_action(
        curr_rij, curr_n, curr_target, zeros
      );
    }

    // calculate the cost of the cheapest n-1 remaining planning units,
    // if the solution contains fewer planning units than the highest target
    if ((n_pu_selected + 1) < max_target) {
      n_additional_units_needed = max_target - n_pu_selected - 1;
      costs_rem_pu_sorted = costs_rem_pu;
      std::partial_sort(
        costs_rem_pu_sorted.data(), costs_rem_pu_sorted.data() +
        n_additional_units_needed,
        costs_rem_pu_sorted.data() + costs_rem_pu_sorted.size());
      curr_min_feasible_pu_cost = budget - std::accumulate(
        costs_rem_pu_sorted.data(),
        costs_rem_pu_sorted.data() + n_additional_units_needed, solution_cost);
    } else {
      curr_min_feasible_pu_cost = std::numeric_limits<double>::infinity();
    }

    // calculate cost-effectiveness score for each planning unit
    curr_ce = std::numeric_limits<double>::lowest();
    for (std::size_t i = 0; i < n_pu; ++i) {
      if (solution_rem_pu[i]) {
        if (
          ((pu_costs[i] <= curr_min_feasible_pu_cost) ||
            ((std::abs(pu_costs[i] - costs_rem_pu_sorted[0])) < 1.0e-15)) &&
          ((pu_costs[i] + solution_cost) <= budget)
        ) {
          // update rij data for calculating objective value
          for (std::size_t j = 0; j < n_f; ++j) {
            extra_n[j] = static_cast<std::size_t>(rij(j, i) > 1.0e-6);
            curr_rij[j][curr_n[j]] = rij(j, i);
          }
          // calculate objective value with i'th planning unit included
          obj = approx_expected_value_of_action(
            curr_rij, curr_n, curr_target, extra_n
          );

          // calculate cost effectiveness
          new_ce = (obj - prev_obj) / pu_costs[i];
          // if i'th planning unit is more cost effective then prevous ones,
          // then make it the candidate for selection in this iteration
          // of the greedy alogirthm
          if (new_ce > curr_ce) {
            curr_ce = new_ce;
            curr_idx = i;
            curr_obj = obj;
          }
        }
      }
    }

    // update loop variables
    solution_cost += pu_costs[curr_idx];
    prev_obj = curr_obj;
    costs_rem_pu[curr_idx] = std::numeric_limits<double>::max();
    solution[curr_idx] = true;
    solution_rem_pu[curr_idx] = false;
    for (std::size_t j = 0; j < n_f; ++j) {
      if (rij(j, curr_idx) > 1.0e-6) {
        curr_rij[j][curr_n[j]] = rij(j, curr_idx);
        ++curr_n[j];
      }
    }
    ++n_pu_selected;
    --n_pu_remaining;
    if (curr_idx == min_cost_pu_idx) {
      min_cost_pu_idx =
        std::min_element(costs_rem_pu.data(),
                         costs_rem_pu.data() + costs_rem_pu.size()) -
        costs_rem_pu.data();
    }
  }

  // Exports
  return;
}

// [[Rcpp::export]]
Rcpp::List rcpp_greedy_heuristic_prioritization(
  Eigen::MatrixXd rij,
  Eigen::VectorXd pu_costs,
  Eigen::VectorXd pu_locked_in,
  Eigen::VectorXd pu_locked_out,
  Eigen::VectorXi target,
  double budget) {
  // initialization
  std::vector<bool> solution(rij.cols());
  greedy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget, solution);
  // export
  return Rcpp::List::create(
    Rcpp::Named("x") = Rcpp::wrap(solution),
    Rcpp::Named("objval") = expected_value_of_action(solution, rij, target));
}
