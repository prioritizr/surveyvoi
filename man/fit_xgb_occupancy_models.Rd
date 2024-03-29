% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_xgb_occupancy_models.R
\name{fit_xgb_occupancy_models}
\alias{fit_xgb_occupancy_models}
\title{Fit boosted regression tree models to predict occupancy}
\usage{
fit_xgb_occupancy_models(
  site_data,
  feature_data,
  site_detection_columns,
  site_n_surveys_columns,
  site_env_vars_columns,
  feature_survey_sensitivity_column,
  feature_survey_specificity_column,
  xgb_tuning_parameters,
  xgb_early_stopping_rounds = rep(20, length(site_detection_columns)),
  xgb_n_rounds = rep(100, length(site_detection_columns)),
  n_folds = rep(5, length(site_detection_columns)),
  n_threads = 1,
  seed = 500,
  verbose = FALSE
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

\item{site_env_vars_columns}{\code{character} names of columns in the
argument to \code{site_data} that contain environmental information
for fitting updated occupancy models based on possible survey outcomes.
Each column should correspond to a different environmental variable,
and contain \code{numeric}, \code{factor}, or \code{character} data.
No missing (\code{NA}) values are permitted in these columns.}

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

\item{xgb_tuning_parameters}{\code{list} object containing the candidate
parameter values for fitting models. Valid parameters include:
\code{"max_depth"}, \code{"eta"}, \code{"lambda"},
\code{"min_child_weight"}, \code{"subsample"}, \code{"colsample_by_tree"},
\code{"objective"}. See documentation for the \code{params} argument in
\code{\link[xgboost:xgb.train]{xgboost::xgb.train()}} for more information.}

\item{xgb_early_stopping_rounds}{\code{numeric} model rounds for parameter
tuning. See \code{\link[xgboost:xgb.train]{xgboost::xgboost()}} for more information.
Defaults to 10 for each feature.}

\item{xgb_n_rounds}{\code{numeric} model rounds for model fitting
See \code{\link[xgboost:xgb.train]{xgboost::xgboost()}} for more information.
Defaults to 100 for each feature.}

\item{n_folds}{\code{numeric} number of folds to split the training
data into when fitting models for each feature.
Defaults to 5 for each feature.}

\item{n_threads}{\code{integer} number of threads to use for parameter
tuning. Defaults to 1.}

\item{seed}{\code{integer} initial random number generator state for model
fitting. Defaults to 500.}

\item{verbose}{\code{logical} indicating if information should be
printed during computations. Defaults to \code{FALSE}.}
}
\value{
A \code{list} object containing:
\describe{

\item{parameters}{\code{list} of \code{list} objects containing the best
tuning parameters for each feature.}

\item{predictions}{\code{\link[tibble:tibble]{tibble::tibble()}} object containing
predictions for each feature.}

\item{performance}{\code{\link[tibble:tibble]{tibble::tibble()}} object containing the
performance of the best models for each feature. It contains the following
columns:

\describe{
\item{feature}{name of the feature.}
\item{train_tss_mean}{
mean TSS statistic for models calculated using training data in
cross-validation.}
\item{train_tss_std}{
standard deviation in TSS statistics for models calculated using training
data in cross-validation.}
\item{train_sensitivity_mean}{
mean sensitivity statistic for models calculated using training data in
cross-validation.}
\item{train_sensitivity_std}{
standard deviation in sensitivity statistics for models calculated using
training data in cross-validation.}
\item{train_specificity_mean}{
mean specificity statistic for models calculated using training data in
cross-validation.}
\item{train_specificity_std}{
standard deviation in specificity statistics for models calculated using
training data in cross-validation.}
\item{test_tss_mean}{
mean TSS statistic for models calculated using test data in
cross-validation.}
\item{test_tss_std}{
standard deviation in TSS statistics for models calculated using test
data in cross-validation.}
\item{test_sensitivity_mean}{
mean sensitivity statistic for models calculated using test data in
cross-validation.}
\item{test_sensitivity_std}{
standard deviation in sensitivity statistics for models calculated using
test data in cross-validation.}
\item{test_specificity_mean}{
mean specificity statistic for models calculated using test data in
cross-validation.}
\item{test_specificity_std}{
standard deviation in specificity statistics for models calculated using
test data in cross-validation.}
}}

}
}
\description{
Estimate probability of occupancy for a set of features in a set of
planning units. Models are fitted using gradient boosted trees (via
\code{\link[xgboost:xgb.train]{xgboost::xgb.train()}}).
}
\details{
This function (i) prepares the data for model fitting, (ii) calibrates
the tuning parameters for model fitting (see \code{\link[xgboost:xgb.train]{xgboost::xgb.train()}}
for details on tuning parameters), (iii) generate predictions using
the best found tuning parameters, and (iv) assess the performance of the
best supported models. These analyses are performed separately for each
feature. For a given feature:

\enumerate{

\item The data are prepared for model fitting by partitioning the data using
k-fold cross-validation (set via argument to \code{n_folds}). The
training and evaluation folds are constructed
in such a manner as to ensure that each training and evaluation
fold contains at least one presence and one absence observation.

\item A grid search method is used to tune the model parameters. The
candidate values for each parameter (specified via \code{parameters}) are
used to generate a full set of parameter combinations, and these
parameter combinations are subsequently used for tuning the models.
To account for unbalanced datasets, the
\code{scale_pos_weight} \code{\link[xgboost:xgb.train]{xgboost::xgboost()}} parameter
is calculated as the mean value across each of the training folds
(i.e. number of absence divided by number of presences per feature).
For a given parameter combination, models are fit using k-fold cross-
validation (via \code{\link[xgboost:xgb.cv]{xgboost::xgb.cv()}}) -- using the previously
mentioned training and evaluation folds -- and the True Skill Statistic
(TSS) calculated using the data held out from each fold is
used to quantify the performance (i.e. \code{"test_tss_mean"} column in
output). These models are also fitted using the
\code{early_stopping_rounds} parameter to reduce time-spent
tuning models. If relevant, they are also fitted using the supplied weights
(per by the argument to \code{site_weights_data}). After exploring the
full set of parameter combinations, the best parameter combination is
identified, and the associated parameter values and models are stored for
later use.

\item The cross-validation models associated with the best parameter
combination are used to generate predict the average probability that the
feature occupies each site. These predictions include sites that have
been surveyed before, and also sites that have not been surveyed before.

\item The performance of the cross-validation models is evaluated.
Specifically, the TSS, sensitivity, and specificity statistics are
calculated (if relevant, weighted by the argument to
\code{site_weights_data}). These performance values are calculated using
the models' training and evaluation folds.

}
}
\examples{
\dontrun{
# set seeds for reproducibility
set.seed(123)

# simulate data for 30 sites, 2 features, and 3 environmental variables
site_data <- simulate_site_data(
  n_sites = 30, n_features = 2, n_env_vars = 3, prop = 0.1)
feature_data <- simulate_feature_data(n_features = 2, prop = 1)

# create list of possible tuning parameters for modeling
parameters <- list(eta = seq(0.1, 0.5, length.out = 3),
                   lambda = 10 ^ seq(-1.0, 0.0, length.out = 3),
                   objective = "binary:logistic")

# fit models
# note that we use 10 random search iterations here so that the example
# finishes quickly, you would probably want something like 1000+
results <- fit_xgb_occupancy_models(
   site_data, feature_data,
   c("f1", "f2"), c("n1", "n2"), c("e1", "e2", "e3"),
   "survey_sensitivity", "survey_specificity",
   n_folds = rep(5, 2), xgb_early_stopping_rounds = rep(100, 2),
   xgb_tuning_parameters = parameters, n_threads = 1)

# print best found model parameters
print(results$parameters)

# print model predictions
print(results$predictions)

# print model performance
print(results$performance, width = Inf)
}
}
