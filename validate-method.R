# load data
load_all()
load("../surveyvoi-manuscript/data/final/results.rda")

# set budget
total_budgets <-
  round(voi_parameters$total_budget_proportion *
        sum(site_data$management_cost + site_data$survey_cost))
b <- total_budgets[1]
n_approx_states <- 100000

# current
v1 <-
  site_data %>%
  dplyr::mutate(
    survey_cost = survey_cost * general_parameters$scale_factor,
    management_cost = management_cost * general_parameters$scale_factor) %>%
  surveyvoi::approx_evdci(
    feature_data = species_data,
    site_occupancy_columns = species_data$code,
    site_probability_columns = paste0(species_data$code, "_model"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = b * general_parameters$scale_factor,
    site_management_locked_in_column = "protected",
    site_management_locked_out_column = "management_locked_out",
    prior_matrix = prior_matrix,
    optimality_gap = voi_parameters$gap,
    n_approx_states  = n_approx_states)

# current minus cost of additional surveys with no information gained
v2 <-
  site_data %>%
  dplyr::mutate(
    survey_cost = survey_cost * general_parameters$scale_factor,
    management_cost = management_cost * general_parameters$scale_factor) %>%
  surveyvoi::approx_evdci(
    feature_data = species_data,
    site_occupancy_columns = species_data$code,
    site_probability_columns = paste0(species_data$code, "_model"),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_model_specificity_column = "model_specificity",
    feature_preweight_column = "preweight",
    feature_postweight_column = "postweight",
    feature_target_column = "target",
    total_budget = (b - sum(site_data$survey_cost)) * general_parameters$scale_factor,
    site_management_locked_in_column = "protected",
    site_management_locked_out_column = "management_locked_out",
    prior_matrix = prior_matrix,
    optimality_gap = voi_parameters$gap,
    n_approx_states  = n_approx_states)

# sample information

v3 <-
site_data %>%
dplyr::mutate(
    survey_cost = survey_cost * general_parameters$scale_factor,
    management_cost =
      management_cost * general_parameters$scale_factor) %>%
surveyvoi::approx_evdsi(
  feature_data = species_data,
  site_survey_scheme_column = "geo_scheme_388",
  site_occupancy_columns = species_data$code,
  site_probability_columns = paste0(species_data$code, "_model"),
  site_weight_columns = paste0(species_data$code, "_n_samples"),
  site_env_vars_columns = env_names,
  site_survey_cost_column = "survey_cost",
  site_management_cost_column = "management_cost",
  feature_survey_column = "survey",
  feature_survey_sensitivity_column = "survey_sensitivity",
  feature_survey_specificity_column = "survey_specificity",
  feature_model_sensitivity_column = "model_sensitivity",
  feature_model_specificity_column = "model_specificity",
  feature_preweight_column = "preweight",
  feature_postweight_column = "postweight",
  feature_target_column = "target",
  total_budget = b * general_parameters$scale_factor,
  xgb_parameters = occ_models$parameters,
  site_management_locked_in_column = "protected",
  site_management_locked_out_column = "management_locked_out",
  prior_matrix = prior_matrix,
  optimality_gap = voi_parameters$gap,
  n_approx_replicates = voi_parameters$n_approx_replicates,
  n_approx_outcomes_per_replicate = 10,
  method_approx_outcomes = voi_parameters$method_approx_outcomes,
  seed = voi_parameters$seed,
  n_approx_states  = n_approx_states)

# print results
print(v1) # should be highest value
print(v2) # should be lowest value
print(v3) # should be somewhere in between v1 and v2
