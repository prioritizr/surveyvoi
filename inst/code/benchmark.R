# define parameters
n_sites <- 15
n_features <- 1
n_benchmark_replicates = 100

# set RNG seeds for reproducibility
library(RandomFields)
set.seed(505)
RFoptions(seed = 505)

# data
site_data <- simulate_site_data(
  n_sites = n_sites, n_features = n_features,
  proportion_of_sites_missing_data = 0.5)
total_budget <- sum(site_data$management_cost) * 0.6
feature_data <- simulate_feature_data(
  n_features = n_features, proportion_of_survey_features = 1)
n_approx_states <- n_states(n_sites, n_features)

# set microbenchmark function
bench_function <- function(..) {
  set.seed(500)
  approx_expected_value_of_decision_given_current_information(
    site_data = site_data,
    feature_data = feature_data,
    site_occupancy_columns = feature_data$name,
    site_probability_columns = gsub("f", "p", feature_data$name),
    site_management_cost_column = "management_cost",
    feature_survey_sensitivity_column = "survey_sensitivity",
    feature_survey_specificity_column = "survey_specificity",
    feature_model_sensitivity_column = "model_sensitivity",
    feature_alpha_column = "alpha",
    feature_gamma_column = "gamma",
    total_budget = total_budget,
    n_approx_replicates = 25,
    n_approx_states_per_replicate = n_approx_states,
    seed = 100)
  invisible(TRUE)
}

# run benchmark
result <- microbenchmark::microbenchmark(bench_function(),
  times = n_benchmark_replicates)

# print result
print(result)
