#' Create K-folds using site data.
#'
#' Create k-folds given survey data in multiple sites.
#'
#' @param prop_detected `integer` proportion of surveys for each site
#'   within which the species was detected. Each
#'   element corresponds to a different site, and values indicate the
#'   proportion of times a species was detected within a given site.
#'   If a site does not have any detections, then a value of zero should be
#'   used (not `NA`).
#'
#' @param n_total `integer` number of total surveys conducted within
#'   each site.
#'   Each element corresponds to a different site, and values indicate the
#'   number of surveys conducted within the given site.
#'   If a site does not have any non-detections, then a value of zero should be
#'   used (not `NA`).
#'
#' @param n `numeric` number of folds.
#'
#' @param index `integer` indices associated with each site.
#'   Defaults to a sequence ranging from 1 to the cardinality of the
#'   argument to `x` (i.e. `seq_along(x)`).
#'
#' @param seed `numeric` random number generated seed for generating
#'   folds. Defaults to 500.
#'
#' @details
#'  The sites will be stratified into folds will be stratified to ensure that
#'  each fold contains least one detection and one non-detection in the
#'  training and test datasets for subsequent model fitting. Note that
#'  sites with have zero detections and zero non-detections are
#'  randomly allocated to folds.
#'
#' @return A `list` of `list` objects containing the
#'  indices excluded from each fold.
#'
#' @noRd
create_site_folds <- function(
  prop_detected, n_total, n, index = seq_along(prop_detected), seed = 500) {
  # assert arguments are valid
  assertthat::assert_that(
    is.numeric(prop_detected), length(prop_detected) > 0,
    assertthat::noNA(prop_detected),
    all(prop_detected >= 0), all(prop_detected <= 1),
    is.numeric(n_total), length(n_total) > 0,
    all(n_total >= 0), assertthat::noNA(n_total),
    identical(length(n_total), length(n_total)),
    assertthat::is.count(n),
    assertthat::noNA(n),
    is.numeric(index), assertthat::noNA(index),
    identical(length(prop_detected), length(index)),
    assertthat::is.count(seed))
  assertthat::assert_that(
    max(abs(round(n_total) - n_total)) < 1e-10,
    msg = "argument to n_total does not contain whole numbers")
  assertthat::assert_that(sum(round(prop_detected * n_total) > 0) >= n,
    msg = "not enough presences to create the specified number of folds")
  assertthat::assert_that(sum(round((1 - prop_detected) * n_total) > 0) >= n,
    msg = "not enough absence to create the specified number of folds")

  # initialization
  n_det <- round(prop_detected * n_total)
  n_nondet <- round((1 - prop_detected) * n_total)
  site_data <- tibble::tibble(
    idx = index, n_det = n_det, n_nondet = n_nondet, n_total = n_total)
  obs_y <- c(rep(rep(1, length(n_det)), n_det),
             rep(rep(0, length(n_nondet)), n_nondet))
  obs_index <- c(rep(index, n_det), rep(index, n_nondet))
  obs_data <- tibble::tibble(y = obs_y, y2 = obs_y, idx = obs_index,
                             idf = factor(as.character(obs_index)))
  obs_data$y2[obs_data$y < 0.5] <- -1

  # organize site data with observations into folds
  withr::with_seed(seed, {
    # create folds
    obs_data2 <- groupdata2::partition(
      obs_data, p = rep(1 / n, n - 1), num_col = "y2", id_col = "idf",
      list_out = FALSE)
  })

  # find valid fold
  fold_columns <- setdiff(names(obs_data2), names(obs_data))
  found_valid <- FALSE
  for (f in fold_columns) {
    ## format fold data
    obs_data2$fold <- as.integer(as.character(obs_data2[[f]]))
    ## calculate statistics to determine if folding scheme is valid
    n_det_per_fold <-
      stats::aggregate(obs_data2$y, by = list(obs_data2$fold), sum)$x
    n_nondet_per_fold <-
      stats::aggregate(1 - obs_data2$y, by = list(obs_data2$fold), sum)$x
    ## if folding scheme is valid, then keep it
    if (all(n_det_per_fold > 0) && all(n_nondet_per_fold > 0)) {
      found_valid <- TRUE
      break()
    }
  }

  # throw error if no valid folding scheme was found, then throw error
  assertthat::assert_that(found_valid,
    msg = paste("could not find any valid folding schemes that have at least",
                "one detection and non-detection per fold, try again with a",
                "different seed."))

  # determine which fold each site belongs to
  site_data$fold <- vapply(site_data$idx, FUN.VALUE = integer(1), function(x) {
    if (x %in% obs_data2$idx) {
      out <- as.integer(obs_data2$fold[obs_data2$idx == x][[1]])
    } else {
      out <- NA_integer_
    }
  })

  # randomly allocate any sites that are missing fold values
  # (because they have no previous detections or non-detections)
  na_pos <- is.na(site_data$fold)
  if (any(na_pos)) {
    withr::with_seed(seed, {
      site_data$fold[na_pos] <-
        sample(seq_len(n), sum(na_pos), replace = sum(na_pos) > n)
    })
  }

  # extract indices for folds
  site_data2 <- split(site_data, site_data$fold)
  train <- lapply(site_data2, function(i) setdiff(index, i$idx))
  test <- lapply(site_data2, function(i) i$idx)

  # return result
  list(train = train, test = test)
}

#' Is JAGS installed?
#'
#' Check if JAGS is installed.
#'
#' @return A `logical` indicating if JAGS is installed.
#'
#' @noRd
is_jags_installed <- function() {
  x <- try(runjags::findjags(), silent = TRUE)
  if (inherits(x, "try-error")) return(FALSE)
  if (!assertthat::is.string(x)) return(FALSE)
  if (!file.exists(x)) return(FALSE)
  TRUE
}
