# load packages
load_all()
library(xgboost)
load("../surveyvoi-manuscript/test.rda")
spp_col <- "acerspi"

# load data
train_ind <- !is.na(site_data[[spp_col]])
x <- as.matrix(sf::st_drop_geometry(site_data)[train_ind, env_names])
y <- site_data[[spp_col]][train_ind]
spw <- 10

# create folds
dtrain <- xgb.DMatrix(x, label = y)
folds <- create_folds(y, 5)

# tune model
withr::with_seed(500, {
  cv <- xgb.cv(data = dtrain, nrounds = 1000, early_stopping_round = 10,
               nthread = 1, folds = folds$test, train_folds = folds$train,
               metrics = "auc", max_depth = 20, eta = 0.1, verbose = TRUE,
              # subsample = 0.2, colsample_bytree = 0.2,
              seed = 500,
               objective = "binary:logistic")
})

# refit model
models <- lapply(seq_along(folds[[2]]), function(i) {
  withr::with_seed(500, {
    set.seed(500)
  xgb.train(data = xgboost::slice(dtrain, folds$train[[i]]),
            nrounds = cv$best_iteration, nthread = 1, verbose = TRUE,
            metrics = "auc", max_depth = 20, eta = 0.1,
           # subsample = 0.2, colsample_bytree = 0.2,
            objective = "binary:logistic")
  })
})

# generate predictions
test_pred <- lapply(seq_along(models), function(i) {
  withr::with_seed(500, {
    predict(models[[i]], xgboost::slice(dtrain, folds$test[[i]]),
            ntreelimit = cv$best_iteration)
  })
})

# calculate AUC
test_auc <- sapply(seq_along(models), function(i) {
  Metrics::auc(y[folds$test[[i]]], test_pred[[i]])
})

# print results
print(c(mean(test_auc), sd(test_auc)))
print(c(cv$evaluation_log$test_auc_mean[cv$best_iteration],
        cv$evaluation_log$test_auc_std[cv$best_iteration]))
