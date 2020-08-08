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
  double log_obj;
  double curr_cost;
  std::size_t n_additional_units_needed;
  std::size_t curr_idx;
  std::size_t n_pu_selected = std::accumulate(
    solution.begin(), solution.end(), 0);
  std::size_t n_pu_remaining = std::accumulate(
    solution_rem_pu.begin(), solution_rem_pu.end(), 0);
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
  // log data
  Eigen::VectorXd log_pu_costs = (pu_costs.array() + 1.0e-5).array().log();
  Eigen::MatrixXd log_1m_rij = rij;
  log_1m_rij = log_1m_rij.cwiseMax(1.0 - 1.0e-10);
  log_1m_matrix(log_1m_rij);
  Eigen::MatrixXd curr_log_1m_rij = log_1m_rij;
  for (std::size_t i = 0; i < n_pu; ++i) {
    if (!solution[i]) {
      curr_log_1m_rij.col(i).setZero();
    }
  }

  // Main processing
  while(((solution_cost + pu_costs[min_cost_pu_idx]) < budget) &&
  (n_pu_remaining > 0)) {
    // calculate the cost of the cheapest n-1 remaining planning units,
    // if the solution contains fewer planning units than the highest target
    if (n_pu_selected < max_target) {
      n_additional_units_needed = max_target - n_pu_selected;
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
          ((curr_min_feasible_pu_cost >= pu_costs[i]) ||
            ((std::abs(pu_costs[i] - costs_rem_pu_sorted[0])) < 1.0e-15)) &&
          ((pu_costs[i] + solution_cost) <= budget)) {
          // calculate objective with i'th planning unit excluded
          curr_log_1m_rij.col(i) = log_1m_rij.col(i);
          log_obj = log_proxy_expected_value_of_action(curr_log_1m_rij);
          curr_log_1m_rij.col(i).setZero();
          // calculate cost effectiveness
          new_ce = log_obj - log_pu_costs[i]; // note that this is logged
          if (new_ce > curr_ce) {
            curr_ce = new_ce;
            curr_idx = i;
          }
        }
      }
    }


    // update loop variables
    solution_cost += pu_costs[curr_idx];
    costs_rem_pu[curr_idx] = std::numeric_limits<double>::max();
    solution[curr_idx] = true;
    solution_rem_pu[curr_idx] = false;
    curr_log_1m_rij.col(curr_idx) = log_1m_rij.col(curr_idx);
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
