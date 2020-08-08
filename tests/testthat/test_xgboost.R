context("xgboost")

r_xgboost <- function(y, x_train, predict_x, xgb_parameters, xgb_nrounds) {
  set.seed(as.numeric(xgb_parameters$seed))
  model <- xgboost::xgboost(
    data = x_train, label = y, nrounds = xgb_nrounds,
    eta = as.numeric(xgb_parameters$eta),
    objective = xgb_parameters$objective,
    verbose = FALSE)
  set.seed(as.numeric(xgb_parameters$seed))
  predict(model, predict_x)
}

test_that("works", {
  ## set parameters
  set.seed(100)
  n_obs_train <- 10000
  n_obs_predict <- 100
  n_vars <- 8
  xgb_parameters <- list(eta = "1", seed = "0", objective = "binary:logistic")
  xgb_nrounds <- 10
  ## simulate matrix of predictor variables for model fitting
  ## including an intercept term
  x_train <- matrix(rnorm(n_vars * n_obs_train), nrow = n_obs_train)
  x_train[, 1] <- 1
  ## simulate matrix of predictor variables for model prediction
  ## including an intercept term
  x_predict <- matrix(rnorm(n_vars * n_obs_predict), nrow = n_obs_predict)
  x_predict[, 1] <- 1
  ## simulate coefficients
  beta <- rnorm(n_vars)
  ## create training labels
  y <- (x_train %*% beta)
  y <- c(exp(y) / (1 + exp(y)))
  ## create predictions
  yhat <- (x_predict %*% beta)
  yhat <- c(exp(yhat) / (1 + exp(yhat)))
  y[] <- rbinom(length(y), 1, y[])
  ## generate predictions
  r1 <- rcpp_xgboost(y, x_train, x_predict, xgb_parameters, xgb_nrounds)
  r2 <- r_xgboost(y, x_train, x_predict, xgb_parameters, xgb_nrounds)
  ## tests
  expect_equal(r1, r2)
})

test_that("reproducible", {
  ## set parameters
  set.seed(100)
  n_obs_train <- 100
  n_obs_predict <- 3
  n_vars <- 8
  n_reps <- 100
  xgb_parameters <- list(eta = "1", seed = "0", colsample_bytree = "0.2",
                         objective = "binary:logistic", nthreads = "1")
  xgb_nrounds <- 90
  ## simulate matrix of predictor variables for model fitting
  ## including an intercept term
  x_train <- matrix(rnorm(n_vars * n_obs_train), nrow = n_obs_train)
  x_train[, 1] <- 1
  ## simulate matrix of predictor variables for model prediction
  ## including an intercept term
  x_predict <- matrix(rnorm(n_vars * n_obs_predict), nrow = n_obs_predict)
  x_predict[, 1] <- 1
  ## simulate coefficients
  beta <- rnorm(n_vars)
  ## create training labels
  y <- (x_train %*% beta)
  y <- c(exp(y) / (1 + exp(y)))
  ## create predictions
  yhat <- (x_predict %*% beta)
  yhat <- c(exp(yhat) / (1 + exp(yhat)))
  y[] <- rbinom(length(y), 1, y[])
  ## generate predictions
  r <- sapply(seq_len(n_reps), function(i) {
    set.seed(123)
    rcpp_xgboost(y, x_train, x_predict, xgb_parameters, xgb_nrounds)
  })
  ## tests
  for (i in seq_len(n_reps))
    expect_equal(r[, 1], r[, i])
})
