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
    all(prior_matrix >= 0), all(prior_matrix <= 1)),
    identical(ncol(prior_matrix), n_sites),
    identical(nrow(prior_matrix), n_features))
}