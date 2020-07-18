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

validate_xgboost_parameters <- function(x) {
  param_names <- c("scale_pos_weight", "max_depth", "eta", "nrounds",
                   "lambda", "subsample",  "colsample_bytree", "objective",
                   "tree_method")
  lapply(x, function(z) {
    assertthat::assert_that(
      assertthat::is.string(z$objective),
      msg = "a feature is missing the objective parameter in xgb_parameters")
    assertthat::assert_that(
      assertthat::is.number(z$scale_pos_weight),
      msg = paste("a feature is missing the scale_pos_weight parameter",
                  "in xgb_parameters"))
    assertthat::assert_that(
      assertthat::is.number(z$nrounds),
      msg = paste("a feature is missing the nrounds parameter",
                  "in xgb_parameters"))
    if (!is.null(z$tree_method)) {
    assertthat::assert_that(
      assertthat::is.string(z$tree_method),
      z$tree_method %in% c("auto", "exact", "hist", "approx"),
      msg = paste("invalid tree_method paramter in xgb_parameters"))
    }
    extra_names <- names(z)[!names(z) %in% param_names]
    assertthat::assert_that(
      length(extra_names) == 0,
      msg = paste("argument to xgb_parameters has unrecognized parameters:",
                  paste(extra_names, collapse = ",")))
    invisible(TRUE)
  })
  invisible(TRUE)
}
