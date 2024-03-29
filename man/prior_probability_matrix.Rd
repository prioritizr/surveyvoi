% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_probability_matrix.R
\name{prior_probability_matrix}
\alias{prior_probability_matrix}
\title{Prior probability matrix}
\usage{
prior_probability_matrix(
  site_data,
  feature_data,
  site_detection_columns,
  site_n_surveys_columns,
  site_probability_columns,
  feature_survey_sensitivity_column,
  feature_survey_specificity_column,
  feature_model_sensitivity_column,
  feature_model_specificity_column
)
}
\arguments{
\item{site_data}{\code{\link[sf:sf]{sf::sf()}} object with site data.}

\item{feature_data}{\code{\link[base:data.frame]{base::data.frame()}} object with feature data.}

\item{site_detection_columns}{\code{character} names of \code{numeric}
columns in the argument to \code{site_data} that contain the proportion of
surveys conducted within each site that detected each feature.
Each column should correspond to a different feature, and contain
a proportion value (between zero and one). If a site has
not previously been surveyed, a value of zero should be used.}

\item{site_n_surveys_columns}{\code{character} names of \code{numeric}
columns in the argument to \code{site_data} that contain the total
number of surveys conducted for each each feature within each site.
Each column should correspond to a different feature, and contain
a non-negative integer number (e.g. 0, 1, 2, 3). If a site has
not previously been surveyed, a value of zero should be used.}

\item{site_probability_columns}{\code{character} names of \code{numeric}
columns in the argument to \code{site_data} that contain modeled
probabilities of occupancy for each feature in each site.
Each column should correspond to a different feature, and contain
probability data (values between zero and one). No missing (\code{NA})
values are permitted in these columns.}

\item{feature_survey_sensitivity_column}{\code{character} name of the
column in the argument to \code{feature_data} that contains
probability of future surveys correctly detecting a presence of each
feature in a given site (i.e. the sensitivity of the survey methodology).
This column should have \code{numeric} values that are between zero and
one. No missing (\code{NA}) values are permitted in this column.}

\item{feature_survey_specificity_column}{\code{character} name of the
column in the argument to \code{feature_data} that contains
probability of future surveys correctly detecting an absence of each
feature in a given site (i.e. the specificity of the survey methodology).
This column should have \code{numeric} values that are between zero and
one. No missing (\code{NA}) values are permitted in this column.}

\item{feature_model_sensitivity_column}{\code{character} name of the
column in the argument to \code{feature_data} that contains
probability of the initial models correctly predicting a presence of each
feature in a given site (i.e. the sensitivity of the models).
This column should have \code{numeric} values that are between zero and
one. No missing (\code{NA}) values are permitted in this column.
This should ideally be calculated using
\code{\link[=fit_xgb_occupancy_models]{fit_xgb_occupancy_models()}} or
\code{\link[=fit_hglm_occupancy_models]{fit_hglm_occupancy_models()}}.}

\item{feature_model_specificity_column}{\code{character} name of the
column in the argument to \code{feature_data} that contains
probability of the initial models correctly predicting an absence of each
feature in a given site (i.e. the specificity of the models).
This column should have \code{numeric} values that are between zero and
one. No missing (\code{NA}) values are permitted in this column.
This should ideally be calculated using
\code{\link[=fit_xgb_occupancy_models]{fit_xgb_occupancy_models()}} or
\code{\link[=fit_hglm_occupancy_models]{fit_hglm_occupancy_models()}}.}
}
\value{
A \code{matrix} object containing the prior probabilities of each
feature occupying each site. Each row corresponds to a different
feature and each column corresponds to a different site.
}
\description{
Create prior probability matrix for the value of information analysis.
}
\examples{
# set seeds for reproducibility
set.seed(123)

# load example site data
data(sim_sites)
print(sim_sites)

# load example feature data
data(sim_features)
print(sim_features)

# calculate prior probability matrix
prior_matrix <- prior_probability_matrix(
  sim_sites, sim_features,
  c("f1", "f2", "f3"), c("n1", "n2", "n3"), c("p1", "p2", "p3"),
  "survey_sensitivity", "survey_specificity",
  "model_sensitivity", "model_specificity")

# preview prior probability matrix
print(prior_matrix)
}
