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

validate_xgb_parameters <- function(x, n_folds) {
  vapply(seq_along(xgb_parameters), FUN.VALUE = logival(1), function(i) {
    out <- xgb_parameters[[i]]$nrounds
    if (is.null(out))
      stop(paste0("argument to xgb_parameters[[", i,
                  "]] is missing nrounds element"))
    if (!identical(length(out), n_folds[[i]]))
      stop(paste0("argument to xgb_parameters[[", i,
                  "]]$nrounds does not have a value for each training fold"))
    invisible(TRUE)
  })
  vapply(seq_along(xgb_parameters), FUN.VALUE = logical(1), function(i) {
    out <- xgb_parameters[[i]]$nrounds
    if (is.null(out))
      stop(paste0("argument to xgb_parameters[[", i,
                  "]] is missing scale_pos_weight element"))
    if (!identical(length(out), n_folds[[i]]))
      stop(paste0("argument to xgb_parameters[[", i, "]]$scale_pos_weight",
                  "does not have a value for each training fold"))
    invisible(TRUE)
  })
  vapply(xgb_parameters, FUN.VALUE = logical(1), function(x) {
    out <- x[names(x) != "nrounds"]
    xgb_param_names <- c("max_depth", "eta", "lambda",
                         "subsample", "colsample_bytree", "objective")
    extra_names <- names(out)[!names(out) %in% xgb_param_names]
    assertthat::assert_that(
      length(extra_names) == 0,
      msg = paste0("argument to xgb_parameters has unrecognized parameters: ,",
                   paste(extra_names, collapse = ",")))
    invisible(TRUE)
  })
  invisible(TRUE)
}
