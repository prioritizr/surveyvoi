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
    static_cast<std::size_t>(target.maxCoeff());
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
  double curr_min_feasible_pu_cost;
  double curr_ce;
  double new_ce;
  double new_obj;
  double curr_alt_obj;
  double curr_cost;
  Eigen::VectorXd costs_rem_pu = pu_costs;
  Eigen::VectorXd costs_rem_pu_sorted;
  std::vector<bool> new_solution;
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

    // find most expensive planning unit in the solution
    curr_cost = std::numeric_limits<double>::lowest();
    for (std::size_t i = 0; i < n_pu; ++i) {
      if (solution_rem_pu[i]) {
        if (pu_costs[i] > curr_cost) {
          curr_cost = pu_costs[i];
          curr_idx = i;
        }
      }
    }

    // if the most expensive planning unit doesn't violate feasbility,
    // then find the least cost-effective planning unit to exclude
    if ((curr_min_feasible_pu_cost + curr_cost) <= budget) {
      curr_ce = std::numeric_limits<double>::max();
      for (std::size_t i = 0; i < n_pu; ++i) {
        if (solution_rem_pu[i]) {
          // calculate objective with i'th planning unit excluded
          new_solution = solution;
          new_solution[i] = false;
          new_obj = expected_value_of_action(new_solution, rij, target);
          // calculate cost effectiveness
          new_ce = (curr_obj - new_obj) / pu_costs[i];
          if (new_ce < curr_ce) {
            curr_alt_obj = new_obj;
            curr_ce = new_ce;
            curr_idx = i;
          }
        }
      }
    } else {
      // calculate objective with most costly planning unit excluded
      new_solution = solution;
      new_solution[curr_idx] = false;
      curr_alt_obj = expected_value_of_action(new_solution, rij, target);
    }

    // update loop variables
    curr_obj = curr_alt_obj;
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
