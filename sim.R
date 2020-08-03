# simulate site data
n_sites <- 100
n_features <- 5
site_data <- simulate_site_data(
  n_sites = n_sites, n_features = site_data,
  proportion_of_sites_missing_data =  0.5,
  n_env_vars = 2, output_probabilities = FALSE)
feature_data <- simulate_feature_data(
  n_features = n_features,n_features proportion_of_survey_features = 1)
for (i in seq_len(nrow(feature_data))) {
  site_data[[paste0("p", i)]] <- runif(nrow(site_data))
}
xgb_parameters <- list(list(
  max_depth =  2,
  eta = 0.2,
  lambda = 0.001,
  subsample = 0.6,
  colsample_bytree = 0.6,
  scale_pos_weight = 1,
  nrounds = 5,
  objective = "binary:logistic"))[rep(1, nrow(feature_data))]

# define budget
total_budget <- 0.6 * sum(site_data$management_cost)
site_data$survey_cost <- 0
site_data$candidate <- !is.na(site_data$f1)
site_data$scheme <- site_data$candidate & (cumsum(site_data$candidate) <= 2)

evd_ci <- evdci(
  site_data = site_data,
  feature_data = feature_data,
  site_occupancy_columns = paste0("f", seq_len(nrow(feature_data))),
  site_probability_columns = paste0("p", seq_len(nrow(feature_data))),
  site_management_cost_column = "management_cost",
  feature_survey_sensitivity_column = "survey_sensitivity",
  feature_survey_specificity_column = "survey_specificity",
  feature_model_sensitivity_column = "model_sensitivity",
  feature_model_specificity_column = "model_specificity",
  feature_target_column = "target",
  total_budget = total_budget)

evd_si <- approx_evdsi(
  site_data = site_data,
  feature_data = feature_data,
  site_occupancy_columns = paste0("f", seq_len(nrow(feature_data))),
  site_probability_columns = paste0("p", seq_len(nrow(feature_data))),
  site_env_vars_columns = c("e1", "e2"),
  site_survey_scheme_column = "scheme",
  site_management_cost_column = "management_cost",
  site_survey_cost_column = "survey_cost",
  feature_survey_column = "survey",
  feature_survey_sensitivity_column = "survey_sensitivity",
  feature_survey_specificity_column = "survey_specificity",
  feature_model_sensitivity_column = "model_sensitivity",
  feature_model_specificity_column = "model_specificity",
  feature_target_column = "target",
  total_budget = total_budget,
  xgb_parameters = xgb_parameters,
  n_approx_outcomes_per_replicate = 1000,
  n_approx_replicates = 1,
  xgb_n_folds = rep(2, nrow(feature_data)))

# print results
print(evd_ci)
print(evd_si)
