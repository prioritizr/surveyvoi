validate_site_detection_data <- function(site_data, column_names,
                                         check_zeros = TRUE) {
  ## check that data are numeric
  is_valid <- sapply(column_names, function(x) is.numeric(site_data[[x]]))
  assertthat::assert_that(
    all(is_valid),
    msg = paste0(paste_list(paste0("\"", column_names[!is_valid], "\"")),
        " columns in site_data are not numeric"))
  ## check that data have no NA values
  is_valid <- sapply(column_names, function(x) assertthat::noNA(site_data[[x]]))
  assertthat::assert_that(
    all(is_valid),
    msg = paste0(paste_list(paste0("\"", column_names[!is_valid], "\"")),
         " columns in site_data have NA values"))
  ## check that data have valid proportion values
  is_valid <- sapply(column_names, function(x) {
    all(site_data[[x]] >= 0 & site_data[[x]] <= 1)
  })
  assertthat::assert_that(
    all(is_valid),
    msg = paste0(paste_list(paste0("\"", column_names[!is_valid], "\"")),
         " columns in site_data must have values >= 0 and <= 1"))
  ## check that each species has been detected at least once
  if (isTRUE(check_zeros)) {
    is_valid <- sapply(column_names, function(x) max(site_data[[x]]) > 0)
    assertthat::assert_that(
      all(is_valid),
      msg = paste0(paste_list(paste0("\"", column_names[!is_valid], "\"")),
         " columns in site_data need a maximum value > 0",
         " (i.e. it needs to be detected at least once)"))
  }
  invisible(TRUE)
}

validate_site_n_surveys_data <- function(site_data, column_names,
                                         check_zeros = TRUE) {
  ## check that data are integer
  is_valid <- sapply(column_names, function(x) {
    max(abs(round(site_data[[x]]) - site_data[[x]])) < 1e-10
  })
  assertthat::assert_that(
    all(is_valid),
    msg = paste0(paste_list(paste0("\"", column_names[!is_valid], "\"")),
        " columns in site_data are not whole numbers"))
  ## check that data have no NA values
  is_valid <- sapply(column_names, function(x) assertthat::noNA(site_data[[x]]))
  assertthat::assert_that(
    all(is_valid),
    msg = paste0(paste_list(paste0("\"", column_names[!is_valid], "\"")),
         " columns in site_data have NA values"))
  ## check that data have values >= 0
  is_valid <- sapply(column_names, function(x) all(site_data[[x]] >= 0))
  assertthat::assert_that(
    all(is_valid),
    msg = paste0(paste_list(paste0("\"", column_names[!is_valid], "\"")),
         " columns in site_data have values < 0"))
  ## check that each species has been detected at least once
  if (isTRUE(check_zeros)) {
    is_valid <- sapply(column_names, function(x) max(site_data[[x]]) > 0)
    assertthat::assert_that(
      all(is_valid),
      msg = paste0(paste_list(paste0("\"", column_names[!is_valid], "\"")),
         " columns in site_data need a maximum value > 0",
         " (i.e. it needs to be surveyed at least in one site)"))
  }
  invisible(TRUE)
}

validate_site_probability_data <- function(
  site_data, site_probability_columns) {
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
    is.matrix(prior_matrix), is.numeric(c(prior_matrix)),
    all(is.finite(c(prior_matrix)),
    all(prior_matrix > 0), all(prior_matrix < 1)),
    identical(ncol(prior_matrix), n_sites),
    identical(nrow(prior_matrix), n_features))
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

paste_list <- function(x) {
  if (length(x) == 1)
    return(x)
  if (length(x) == 2)
    return(paste(x[1], "and", x[2]))
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}
