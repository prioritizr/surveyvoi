#' @include internal.R
NULL

#' Plot conservation value
#'
#' Create a plot showing the relative conservation value for different
#' conservation features.
#'
#' @inheritParams approx_evdsi
#'
#' @param plot \code{logical} should the plot be returned? If \code{FALSE}
#'   the underlying data (i.e. a \code{tibble}{tibble} object) is returned
#'   instead of the plot. Defaults to \code{TRUE}.
#'
#' @details This function creates a plot to visualize the conservation value
#'   for each feature. The x-axis shows the number of protected sites that
#'   are occupied by a given feature, and the y-axis represents the amount of
#'   conservation value with that number of protected sites. Lines correspond
#'   to features. Each feature will have a line that is comprised of a
#'   solid component and a dashed component. The solid component of the line
#'   corresponds to the range of values based on recorded
#'   presences and absences in the site data. The dashed component of the
#'   line corresponds to the additional values that could
#'   be obtained if additional surveys were conducted and the feature was
#'   detected in all sites that are currently unsurveyed. Features that are
#'   associated with greater values (greater y-axis
#'   values) at lower numbers of protected sites (smaller x-axis values) will,
#'   broadly speaking, be considered more important for protection.
#'
#' @return \code{\link[ggplot2]{ggplot}} plot object.
#'
#' @examples
#' # set seeds for reproducibility
#' library(RandomFields)
#' set.seed(123)
#' RFoptions(seed = 123)
#'
#' # simulate data
#' site_data <- simulate_site_data(n_sites = 5, n_features = 4, prop = 0.5)
#' feature_data <- simulate_feature_data(n_features = 4, prop = 0.5)
#'
#' # preview simulated data
#' print(site_data)
#' print(feature_data)
#'
#' # plot conservation value for each feature
#' plot_conservation_value(
#'  site_data = site_data,
#'  feature_data = feature_data,
#'  site_occupancy_columns = paste0("f", seq_len(4)),
#'  feature_preweight_column = "preweight",
#'  feature_postweight_column  = "postweight",
#'  feature_target_column = "target")
#'
#' @export
plot_conservation_value <- function(
  site_data,
  feature_data,
  site_occupancy_columns,
  feature_preweight_column,
  feature_postweight_column,
  feature_target_column,
  plot = TRUE) {
  # assert arguments are valid
  assertthat::assert_that(
  ## site_data
  inherits(site_data, "sf"), ncol(site_data) > 0,
  nrow(site_data) > 0,
  ## feature_data
  inherits(feature_data, "data.frame"), ncol(feature_data) > 0,
  nrow(feature_data) > 0,
  ## site_occupancy_columns
  is.character(site_occupancy_columns),
  identical(nrow(feature_data), length(site_occupancy_columns)),
  assertthat::noNA(site_occupancy_columns),
  all(assertthat::has_name(site_data, site_occupancy_columns)),
  ## feature_preweight_column
  assertthat::is.string(feature_preweight_column),
  all(assertthat::has_name(feature_data, feature_preweight_column)),
  is.numeric(feature_data[[feature_preweight_column]]),
  assertthat::noNA(feature_data[[feature_preweight_column]]),
  all(feature_data[[feature_preweight_column]] >= 0),
  ## feature_postweight_column
  assertthat::is.string(feature_postweight_column),
  all(assertthat::has_name(feature_data, feature_postweight_column)),
  is.numeric(feature_data[[feature_postweight_column]]),
  assertthat::noNA(feature_data[[feature_postweight_column]]),
  all(feature_data[[feature_postweight_column]] >= 0),
  ## feature_target_column
  assertthat::is.string(feature_target_column),
  all(assertthat::has_name(feature_data, feature_target_column)),
  is.numeric(feature_data[[feature_target_column]]),
  assertthat::noNA(feature_data[[feature_target_column]]),
  all(feature_data[[feature_target_column]] >= 0),
  ## plot
  assertthat::is.flag(plot), assertthat::noNA(plot))
  ## validate rij values
  validate_site_occupancy_data(site_data, site_occupancy_columns)

  # drop spatial information
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)

  # preliminary processing
  rij <- t(as.matrix(site_data[, site_occupancy_columns]))
  n_survey_held <- rowSums(rij, na.rm = TRUE)
  rij2 <- rij
  rij2[is.na(rij2)] <- 1
  n_unsurveyed_held <- rowSums(rij2, na.rm = TRUE)

  # create data to plot
  d <- plyr::ldply(seq_len(nrow(feature_data)), function(i) {
    survey_held <- seq(0, n_survey_held[i])
    survey_value <- vapply(
      survey_held, rcpp_conservation_value_amount, numeric(1),
      preweight = feature_data[[feature_preweight_column]][[i]],
      postweight = feature_data[[feature_postweight_column]][[i]],
      target = feature_data[[feature_target_column]][[i]],
      total = nrow(site_data))
    unsurveyed_held <- seq(n_survey_held[i], n_unsurveyed_held[i])
    unsurveyed_value <- vapply(
      unsurveyed_held, rcpp_conservation_value_amount, numeric(1),
      preweight = feature_data[[feature_preweight_column]][[i]],
      postweight = feature_data[[feature_postweight_column]][[i]],
      target = feature_data[[feature_target_column]][[i]],
      total = nrow(site_data))
  tibble::tibble(feature = as.character(i),
                 held = c(survey_held, unsurveyed_held),
                 value = c(survey_value, unsurveyed_value),
                 type = c(rep("observed", length(survey_value)),
                          rep("potential", length(unsurveyed_value))))
  }) %>% tibble::as_tibble()

  # if not plot then return result
  if (!isTRUE(plot)) return(d)

  # create plot
  p <-
    ggplot2::ggplot(
      data = d,
      mapping = ggplot2::aes(x = held, y = value, color = feature,
                        linetype = type)) +
    ggplot2::geom_line() +
    ggplot2::scale_linetype_manual(
      values = c("observed" = "solid", "potential" = "dashed")) +
    ggplot2::labs(
      color = "Feature",
      linetype = "Data type",
      x = "Number of occupied sites held in prioritization",
      y = "Conservation value")

  # return plot
  p
}
