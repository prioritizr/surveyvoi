% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_hglm_occupancy_models.R
\name{fit_hglm_occupancy_models}
\alias{fit_hglm_occupancy_models}
\title{Fit hierarchical generalized linear models to predict occupancy}
\usage{
fit_hglm_occupancy_models(
  site_data,
  feature_data,
  site_detection_columns,
  site_n_surveys_columns,
  site_env_vars_columns,
  feature_survey_sensitivity_column,
  feature_survey_specificity_column,
  jags_n_samples = rep(10000, length(site_detection_columns)),
  jags_n_burnin = rep(1000, length(site_detection_columns)),
  jags_n_thin = rep(100, length(site_detection_columns)),
  jags_n_adapt = rep(1000, length(site_detection_columns)),
  jags_n_chains = rep(4, length(site_detection_columns)),
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

\item{jags_n_samples}{\code{integer} number of sample to generate per chain
for MCMC analyses.
See documentation for the \code{sample} parameter
in \code{\link[runjags:run.jags]{runjags::run.jags()}} for more information).
Defaults to 10,000 for each feature.}

\item{jags_n_burnin}{\code{integer} number of warm up iterations for MCMC
analyses.
See documentation for the \code{burnin} parameter
in \code{\link[runjags:run.jags]{runjags::run.jags()}} for more information).
Defaults to 10,000 for each feature.}

\item{jags_n_thin}{\code{integer} number of thinning iterations for MCMC
analyses.
See documentation for the \code{thin} parameter
in \code{\link[runjags:run.jags]{runjags::run.jags()}} for more information).
Defaults to 100 for each feature.}

\item{jags_n_adapt}{\code{integer} number of adapting iterations for MCMC
analyses.
See documentation for the \code{adapt} parameter
in \code{\link[runjags:run.jags]{runjags::run.jags()}} for more information).
Defaults to 1,000 for each feature.}

\item{jags_n_chains}{\code{integer} total number of chains for MCMC analyses.
See documentation for the \code{n.chains} parameter
in \code{\link[runjags:run.jags]{runjags::run.jags()}} for more information).
Defaults to 4 for each fold for each feature.}

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

\item{models}{\code{list} of \code{list} objects containing the models.}

\item{predictions}{\code{\link[tibble:tibble]{tibble::tibble()}} object containing
predictions for each feature.}

\item{performance}{\code{\link[tibble:tibble]{tibble::tibble()}} object containing the
performance of the best models for each feature. It contains the following
columns:

\describe{
\item{feature}{name of the feature.}
\item{max_mpsrf}{maximum multi-variate potential scale reduction factor
(MPSRF) value for the models. A MPSRF value less than 1.05 means that all
coefficients in a given model have converged, and so a value less than
1.05 in this column means that all the models fit for a given feature
have successfully converged.}
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
planning units. Models are fitted as hierarchical generalized linear models
that account for for imperfect detection (following Royle & Link 2006)
using JAGS (via \code{\link[runjags:run.jags]{runjags::run.jags()}}). To limit over-fitting,
covariate coefficients are sampled using a Laplace prior distribution
(equivalent to L1 regularization used in machine learning contexts)
(Park & Casella 2008).
}
\details{
This function (i) prepares the data for model fitting,
(ii) fits the models, and (iii) assesses the performance of the
models. These analyses are performed separately for each
feature. For a given feature:

\enumerate{

\item The data are prepared for model fitting by partitioning the data using
k-fold cross-validation (set via argument to \code{n_folds}). The
training and evaluation folds are constructed
in such a manner as to ensure that each training and evaluation
fold contains at least one presence and one absence observation.

\item A model for fit separately for each fold (see
\code{inst/jags/model.jags} for model code). To assess convergence,
the multi-variate potential scale reduction factor
(MPSRF) statistic is calculated for each model.

\item The performance of the cross-validation models is evaluated.
Specifically, the TSS, sensitivity, and specificity statistics are
calculated (if relevant, weighted by the argument to
\code{site_weights_data}). These performance values are calculated using
the models' training and evaluation folds. To assess convergence,
the maximum MPSRF statistic for the models fit for each feature
is calculated.

}
}
\section{Dependencies}{

This function requires the
\href{https://mcmc-jags.sourceforge.io/}{JAGS software} to be installed.
For information on installing the JAGS software, please consult
the documentation for the \pkg{rjags} package.
}

\examples{
\dontrun{
# set seeds for reproducibility
set.seed(123)

# simulate data for 200 sites, 2 features, and 3 environmental variables
site_data <- simulate_site_data(n_sites = 30, n_features = 2, prop = 0.1)
feature_data <- simulate_feature_data(n_features = 2, prop = 1)

# print JAGS model code
cat(readLines(system.file("jags", "model.jags", package = "surveyvoi")),
    sep = "\n")

# fit models
# note that we use a small number of MCMC iterations so that the example
# finishes quickly, you probably want to use the defaults for real work
results <- fit_hglm_occupancy_models(
   site_data, feature_data,
   c("f1", "f2"), c("n1", "n2"), c("e1", "e2", "e3"),
   "survey_sensitivity", "survey_specificity",
   n_folds = rep(5, 2),
   jags_n_samples = rep(250, 2), jags_n_burnin = rep(250, 2),
   jags_n_thin = rep(1, 2), jags_n_adapt = rep(100, 2),
   n_threads = 1)

# print model predictions
print(results$predictions)

# print model performance
print(results$performance, width = Inf)
}
}
\references{
Park T & Casella G (2008) The Bayesian lasso.
\emph{Journal of the American Statistical Association}, 103: 681--686.

Royle JA & Link WA (2006) Generalized site occupancy models allowing for
false positive and false negative errors. \emph{Ecology}, 87: 835--841.
}
