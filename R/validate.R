validate_site_occupancy_data <- function(site_data, site_occupancy_columns) {
  assertthat::assert_that(
    all(sapply(site_occupancy_columns,
               function(x) is.numeric(site_data[[x]]))),
    msg = "site_data values in site_occupancy_columns must be numeric")
  assertthat::assert_that(
    all(sapply(site_occupancy_columns,
               function(x) max(site_data[[x]], na.rm = TRUE) <= 1)),
    msg = "site_data values in site_occupancy_columns must be <= 1")
  assertthat::assert_that(
    all(sapply(site_occupancy_columns,
               function(x) min(site_data[[x]], na.rm = TRUE) >= 0)),
    msg = "site_data values in site_occupancy_columns must be >= 0")
  assertthat::assert_that(
    all(sapply(site_occupancy_columns,
               function(x) any(site_data[[x]] == 0, na.rm = TRUE))),
    msg = paste("site_data values in site_occupancy_columns require at",
                "least one absence per feature"))
  assertthat::assert_that(
    all(sapply(site_occupancy_columns,
               function(x) any(site_data[[x]] == 1, na.rm = TRUE))),
    msg = paste("site_data values in site_occupancy_columns require at",
                "least one presence per feature"))
  invisible(TRUE)
}

validate_site_prior_data <- function(site_data, site_probability_columns) {
  assertthat::assert_that(
    all(sapply(site_probability_columns,
               function(x) is.numeric(site_data[[x]]))),
    msg = "site_data values in site_probability_columns must be numeric")
  assertthat::assert_that(
    all(sapply(site_probability_columns,
               function(x) all(!is.na(site_data[[x]])))),
    msg = "site_data values in site_probability_columns must not be NA")
  assertthat::assert_that(
    all(sapply(site_probability_columns,
               function(x)
                 all(site_data[[x]] >= 0) && all(site_data[[x]] <= 1))),
    msg = paste("site_data values in site_probability_columns must be between",
                "0 and 1"))
  invisible(TRUE)
}

validate_prior_data <- function(prior_matrix, n_sites, n_features) {
  assertthat::assert_that(
    is.matrix(prior_matrix), is.numeric(prior_matrix),
    all(is.finite(c(prior_matrix)),
    all(prior_matrix > 0), all(prior_matrix < 1)),
    identical(ncol(prior_matrix), n_sites),
    identical(nrow(prior_matrix), n_features))
}

validate_site_weight_data <- function(site_data, site_occupancy_columns,
  site_weight_columns) {
  assertthat::assert_that(
    is.character(site_weight_columns),
    identical(length(site_weight_columns), length(site_occupancy_columns)),
    all(assertthat::has_name(site_data, site_weight_columns)),
    assertthat::noNA(site_weight_columns))
  assertthat::assert_that(
    all(sapply(site_weight_columns,
               function(x) is.numeric(site_data[[x]]))),
    msg = "site_data values in site_weight_columns must be numeric")
  assertthat::assert_that(
    all(sapply(site_weight_columns,
               function(x) all(is.finite(site_data[[x]])))),
    msg = "site_data values in site_weight_columns must not be NA")
}

validate_target_data <- function(feature_data, feature_target_column) {
  assertthat::assert_that(
    inherits(feature_data, "data.frame"),
    assertthat::is.string(feature_target_column),
    assertthat::has_name(feature_data, feature_target_column))
  assertthat::assert_that(
    all(vapply(feature_data[[feature_target_column]], assertthat::is.count,
           logical(1))),
    msg = paste("feature target values must be count values",
                "(i.e. integer values >= 1"))
  assertthat::assert_that(
    dplyr::n_distinct(feature_data[[feature_target_column]]) == 1,
    msg = paste("all features must have exactly the same target value"))
  invisible(TRUE)
}

validate_xgboost_tuning_parameters <- function(x) {
  param_names <- c("max_depth", "eta", "lambda", "subsample",
                   "colsample_bytree", "objective", "tree_method")
  assertthat::assert_that(
    all(names(x) %in% param_names),
    msg = paste("argument to xgb_tuning_parameters has unrecognised elements:",
               paste(setdiff(names(x), param_names),
                     collapse = ", ")))
  if ("tree_method" %in% names(x))
    assertthat::assert_that(
      all(x$tree_method %in% c("auto", "hist", "exact", "approx")),
      msg = "argument to xgb_tuning_parameters has invalid tree_method value")
  invisible(TRUE)
}
