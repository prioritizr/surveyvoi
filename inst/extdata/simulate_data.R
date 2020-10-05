# Initialization
## load packages
library(sf)
library(RandomFields)

## load functions
source("R/simulate_feature_data.R")
source("R/simulate_site_data.R")

## set seed for reproducibility
seed <- 600

## set simulation parameters
number_features <- 3
number_sites <- 15
proportion_of_survey_features <- 1
proportion_of_sites_missing_data <- 0.75
n_env_vars <- 2
survey_cost_intensity <- 10
survey_cost_radius <- 0.5
management_cost_intensity <- 80
management_cost_radius <- 0.5
max_number_surveys_per_site <- 5
output_probabilities <- TRUE

# Simulate data
## features
set.seed(seed)
RFoptions(seed = seed)
sim_features <- simulate_feature_data(
  number_features, proportion_of_survey_features)

## sites
set.seed(seed)
RFoptions(seed = seed)
sim_sites <- simulate_site_data(
  n_sites = number_sites,
  n_features = number_features,
  proportion_of_sites_missing_data = proportion_of_sites_missing_data,
  n_env_vars = n_env_vars,
  survey_cost_intensity = survey_cost_intensity,
  survey_cost_radius = survey_cost_radius,
  management_cost_intensity = management_cost_intensity,
  management_cost_radius = management_cost_radius,
  max_number_surveys_per_site = max_number_surveys_per_site,
  output_probabilities = TRUE)

# Exports
save(sim_sites, file = "data/sim_sites.rda", compress = "xz")
save(sim_features, file = "data/sim_features.rda", compress = "xz")
