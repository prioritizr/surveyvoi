# init
load_all()
load("test.rda")

n_states_per_rep = 10000
budget <- sum(site_data$management_cost) * 0.8

# run calculations
withr::with_seed(o1$seed, {
new_evdsi <- r_approx_expected_value_of_decision_given_survey_scheme_n_states(
  rij = o1$rij, pij = o1$pij, wij = o1$wij,
  survey_features = feature_data$survey,
  survey_sensitivity = feature_data$survey_sensitivity,
  survey_specificity = feature_data$survey_specificity,
  pu_survey_solution = site_data[[o1$site_survey_scheme_column]],
  pu_survey_status = o1$site_survey_status,
  pu_survey_costs = site_data$survey_cost,
  pu_purchase_costs = site_data$management_cost,
  pu_purchase_locked_in = rep(FALSE, nrow(site_data)),
  pu_env_data = o1$ejx,
  xgb_parameters = o1$xgb_parameters,
  xgb_nrounds = o1$xgb_nrounds,
  xgb_train_folds = lapply(o1$xgb_folds, `[[`, "train"),
  xgb_test_folds = lapply(o1$xgb_folds, `[[`, "test"),
  obj_fun_preweight = feature_data$preweight,
  obj_fun_postweight  = feature_data$postweight,
  obj_fun_target = feature_data$target,
  n_approx_obj_fun_points = o1$n_approx_obj_fun_points,
  total_budget = budget,
  optim_gap = o1$optimality_gap,
  n_approx_replicates = 1, n_approx_states_per_replicate = n_states_per_rep)
})

withr::with_seed(o1$seed, {
new_edpi <- r_approx_expected_value_of_decision_given_perfect_info_n_states(
  prior_data = o1$pij, pu_costs = site_data$management_cost,
  pu_locked_in = rep(FALSE, nrow(site_data)),
  preweight = feature_data$preweight,
  postweight = feature_data$postweight,
  target = feature_data$target,
  n_approx_obj_fun_points = o1$n_approx_obj_fun_points,
  budget = budget, gap = o1$optimality_gap,
  n_replicates = 1, n_states_per_replicate = n_states_per_rep)
})

withr::with_seed(o1$seed, {
new_edci <- r_approx_expected_value_of_decision_given_current_info_n_states(
  prior_data = o1$pij, pu_costs = site_data$management_cost,
  pu_locked_in = rep(FALSE, nrow(site_data)),
  preweight = feature_data$preweight,
  postweight = feature_data$postweight,
  target = feature_data$target,
  n_approx_obj_fun_points = o1$n_approx_obj_fun_points,
  budget = budget, gap = o1$optimality_gap,
  n_replicates = 1, n_states_per_replicate = n_states_per_rep)
})

# rcpp answer
print(new_evdsi)
print(new_edpi)
print(new_edci)
