#include "rcpp_prioritization.h"

void prioritization(
  Eigen::MatrixXd& rij, Eigen::VectorXd& pu_costs,
  Eigen::VectorXd& pu_locked_in, Eigen::VectorXd& pu_locked_out,
  Eigen::VectorXi& target, double budget, std::vector<bool>& solution) {
  // Initialization
  // declare constants
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = rij.rows();
  const std::size_t max_target_m1 =
    static_cast<std::size_t>(target.maxCoeff() - 1);
  // reset solution
  std::fill(solution.begin(), solution.end(), true);
  for (std::size_t i = 0; i < n_pu; ++i)
    if (pu_locked_out[i] > 0.5)
      solution[i] = false;
  // calculate objective value of initial solotion
  double curr_obj = expected_value_of_action(solution, rij, target);
  double solution_cost = 0.0;
  for (std::size_t i = 0; i < n_pu; ++i)
    solution_cost += static_cast<double>(solution[i]) * pu_costs[i];
  // declare loop variables
  /// variables to keep track of planning units in the solution that
  /// will be considered for step-wise removal from the algorithm
  std::vector<bool> solution_rem_pu = solution;
  for (std::size_t i = 0; i < n_pu; ++i)
    if (pu_locked_in[i] > 0.5)
      solution_rem_pu[i] = false;
  for (std::size_t i = 0; i < n_pu; ++i)
    if (pu_costs[i] < 1.0e-15)
      solution_rem_pu[i] = false;
  Eigen::VectorXd curr_alt_obj(n_pu);
  Eigen::VectorXd curr_infeasible(n_pu);
  double curr_min_feasible_pu_cost;
  double curr_ce;
  double new_ce;
  Eigen::VectorXd costs_rem_pu = pu_costs;
  Eigen::VectorXd costs_rem_pu_sorted;
  std::vector<bool> curr_solution;
  for (std::size_t i = 0; i < n_pu; ++i)
    if (!solution[i])
      costs_rem_pu[i] = std::numeric_limits<double>::max();
  std::size_t curr_idx;

  // Main processing
  while(solution_cost > budget) {
    // calculate the cost of the cheapest n-1 remaining planning units
    costs_rem_pu_sorted = costs_rem_pu;
    std::partial_sort(
      costs_rem_pu_sorted.data(), costs_rem_pu_sorted.data() + max_target_m1,
      costs_rem_pu_sorted.data() + costs_rem_pu_sorted.size());
    curr_min_feasible_pu_cost = std::accumulate(
      costs_rem_pu_sorted.data(), costs_rem_pu_sorted.data() + max_target_m1,
      0.0);

    // calculate objective value of solutions associated with dropping each
    // remaining planning unit for the solution
    curr_infeasible.setZero();
    for (std::size_t i = 0; i < n_pu; ++i) {
      // if the i'th unit is in the remaining planning units then assign Inf
      if (!solution_rem_pu[i]) continue;
      // if the cost of the i'th unit plus the n-1 planning units exceeds
      // the budget, then store index in curr_infeasible
      if ((curr_min_feasible_pu_cost + pu_costs[i]) > budget) {
        curr_infeasible[i] = 1.0;
        break;
      }
      // calculate objective value for the solution
      curr_solution = solution;
      curr_solution[i] = false;
      curr_alt_obj[i] = expected_value_of_action(curr_solution, rij, target);
    }

    // determine which planning unit to exclude in this iteration
    if (curr_infeasible.maxCoeff() > 0.5) {
      // select the planning unit which costs too much to yield a feasible
      // solution
      for (std::size_t i = 0; i < n_pu; ++i) {
        if (curr_infeasible[i] > 0.5) {
          curr_idx = i;
          break;
        }
      }
    } else {
      // select the planning unit with the lowest cost-effectiveness
      curr_ce = std::numeric_limits<double>::max();
      for (std::size_t i = 0; i < n_pu; ++i) {
        if (solution_rem_pu[i]) {
          new_ce = (curr_obj - curr_alt_obj[i]) / pu_costs[i];
          if (new_ce < curr_ce) {
            curr_ce = new_ce;
            curr_idx = i;
          }
        }
      }
    }

    // update loop variables
    curr_obj = curr_alt_obj[curr_idx];
    solution_cost -= pu_costs[curr_idx];
    costs_rem_pu[curr_idx] = std::numeric_limits<double>::max();
    solution[curr_idx] = false;
    solution_rem_pu[curr_idx] = false;
  }

  // Exports
  return;
}

// [[Rcpp::export]]
Rcpp::List rcpp_prioritization(
  Eigen::MatrixXd rij,
  Eigen::VectorXd pu_costs,
  Eigen::VectorXd pu_locked_in,
  Eigen::VectorXd pu_locked_out,
  Eigen::VectorXi target,
  double budget) {
  // initialization
  std::vector<bool> solution(rij.cols());

  // main processing
  prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget, solution);

  // export
  return Rcpp::List::create(
    Rcpp::Named("x") = Rcpp::wrap(solution),
    Rcpp::Named("objval") = expected_value_of_action(solution, rij, target));
}
