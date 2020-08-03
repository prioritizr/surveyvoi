#include "rcpp_heuristic_prioritization.h"

void stingy_heuristic_prioritization(
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
  Rcpp::IntegerVector target_values((rij.cols() - target[0]) + 1);
  std::iota(target_values.begin(), target_values.end(), target[0]);
  Eigen::MatrixXd rij_rem_pu = rij;
  for (std::size_t i = 0; i < n_pu; ++i)
    if (!solution[i])
      rij_rem_pu.col(i).setZero();

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
          rij_rem_pu.col(i).setZero();
          new_obj = approx_expected_value_of_action(rij_rem_pu, target_values);
          rij_rem_pu.col(i).array() = rij.col(i).array();
          // calculate cost effectiveness
          new_ce = (curr_obj - new_obj) / pu_costs[i];
          if (new_ce < curr_ce) {
            curr_alt_obj = new_obj;
            curr_ce = new_ce;
            curr_idx = i;
          }
        }
      }
      rij_rem_pu.col(curr_idx).setZero();
    } else {
      // calculate objective with most costly planning unit excluded
      curr_alt_obj = approx_expected_value_of_action(rij_rem_pu, target_values);
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

void greedy_heuristic_prioritization(
  Eigen::MatrixXd& rij, Eigen::VectorXd& pu_costs,
  Eigen::VectorXd& pu_locked_in, Eigen::VectorXd& pu_locked_out,
  Eigen::VectorXi& target, double budget, std::vector<bool>& solution) {
  // Initialization
  // declare constants
  const std::size_t n_pu = rij.cols();
  const std::size_t n_f = rij.rows();
  const std::size_t max_target_max =
    static_cast<std::size_t>(target.maxCoeff());
  // reset solution
  std::fill(solution.begin(), solution.end(), false);
  for (std::size_t i = 0; i < n_pu; ++i)
    if ((pu_locked_in[i] > 0.5) || (pu_costs[i] < 1.0e-15))
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
  double new_obj;
  double curr_alt_obj;
  double curr_cost;
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
  // variables for calculating objective function
  std::size_t curr_max_target =
    std::min(static_cast<double>(max_target_max),
             static_cast<double>(n_pu_selected));
  std::vector<int> target_values_raw;
  Rcpp::IntegerVector target_values((n_pu_selected - curr_max_target) + 1);
  std::iota(target_values.begin(), target_values.end(), curr_max_target);
  // rij matrices
  Eigen::MatrixXd curr_rij(rij.rows(), n_pu_selected);
  for (std::size_t i = 0, d = 0; i < n_pu; ++i) {
    if (solution[i]) {
      curr_rij.col(d) = rij.col(i);
      ++d;
    }
  }
  // calculate objective value of initial solution
  double curr_obj;
  if (n_pu_selected == 0) {
    curr_obj = 0.0;
  } else {
    curr_obj = approx_expected_value_of_action(curr_rij, target_values);
  }

  // Main processing
  while(((solution_cost + pu_costs[min_cost_pu_idx]) < budget) &&
        (n_pu_remaining > 0)) {
    // update targets
    curr_max_target =
      std::min(static_cast<double>(max_target_max),
               static_cast<double>(n_pu_selected + 1));
    target_values_raw.resize(((n_pu_selected + 1) - curr_max_target) + 1);
    std::iota(target_values_raw.begin(), target_values_raw.end(),
              curr_max_target);
    target_values = Rcpp::wrap(target_values_raw);

    // add extra column to curr_rij matrix
    curr_rij.conservativeResize(Eigen::NoChange, n_pu_selected + 1);

    // calculate the cost of the cheapest n-1 remaining planning units,
    // if the solution contains fewer planning units than the highest target
    if (n_pu_selected < max_target_max) {
      costs_rem_pu_sorted = costs_rem_pu;
      std::partial_sort(
        costs_rem_pu_sorted.data(), costs_rem_pu_sorted.data() +
        max_target_max - n_pu_selected,
        costs_rem_pu_sorted.data() + costs_rem_pu_sorted.size());
      curr_min_feasible_pu_cost = budget - std::accumulate(
        costs_rem_pu_sorted.data(),
        costs_rem_pu_sorted.data() + (max_target_max - n_pu_selected),
        solution_cost);
    } else {
      curr_min_feasible_pu_cost = std::numeric_limits<double>::infinity();
    }

    // calculate cost-effectiveness score for each planning unit
    curr_ce = std::numeric_limits<double>::lowest();
    for (std::size_t i = 0; i < n_pu; ++i) {
      if (solution_rem_pu[i]) {
        if ((curr_min_feasible_pu_cost >= pu_costs[i]) &&
            ((pu_costs[i] + solution_cost) <= budget)) {
          // calculate objective with i'th planning unit excluded
          curr_rij.col(n_pu_selected).array() = rij.col(i).array();
          new_obj = approx_expected_value_of_action(curr_rij, target_values);
          // calculate cost effectiveness
          new_ce = new_obj / pu_costs[i];
          if (new_ce > curr_ce) {
            curr_alt_obj = new_obj;
            curr_ce = new_ce;
            curr_idx = i;
          }
        }
      }
    }

    // update loop variables
    curr_obj = curr_alt_obj;
    solution_cost += pu_costs[curr_idx];
    costs_rem_pu[curr_idx] = std::numeric_limits<double>::max();
    solution[curr_idx] = true;
    solution_rem_pu[curr_idx] = false;
    curr_rij.col(n_pu_selected).array() = rij.col(curr_idx).array();;
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
Rcpp::List rcpp_stingy_heuristic_prioritization(
  Eigen::MatrixXd rij,
  Eigen::VectorXd pu_costs,
  Eigen::VectorXd pu_locked_in,
  Eigen::VectorXd pu_locked_out,
  Eigen::VectorXi target,
  double budget) {
  // initialization
  std::vector<bool> solution(rij.cols());
  stingy_heuristic_prioritization(
    rij, pu_costs, pu_locked_in, pu_locked_out, target, budget, solution);

  // export
  return Rcpp::List::create(
    Rcpp::Named("x") = Rcpp::wrap(solution),
    Rcpp::Named("objval") = expected_value_of_action(solution, rij, target));
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
