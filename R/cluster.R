#' @include internal.R
NULL

#' Start parallel processing cluster
#'
#' Start a new cluster for parallel processing.
#'
#' @param n `integer` number of cores.
#'
#' @param names `character` name of objects to send to cluster if
#'  using `PSOCK` cluster.
#'
#' @param type `character` type of cluster. Available options are
#'  `"FORK"` or `"PSOCK"`. Defaults to `NULL` which will correspond to
#'   `"FORK"` on Unix systems, and `"PSOCK"` on other systems.

#' @param envir `environment` with objects. Defaults to parent environment
#'  (`parent.env()`).
#'
#' @return `cluster` object.
#'
#' @noRd
start_cluster <- function(n, names, type = NULL) {
  # assert arguments are valid
  assertthat::assert_that(
    assertthat::is.count(n),
    assertthat::noNA(n),
    is.character(names),
    assertthat::noNA(names),
    inherits(type, c("character", "NULL")))
  # determine cluster type
  if (is.null(type)) {
    type <- ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK")
  }
  if (nchar(Sys.getenv("TESTTHAT")) > 0) {
    # this function only works with PSOCK clusters when tested in testthat v3+
    type <- "PSOCK"
  }
  # validate cluster type
  assertthat::assert_that(
    assertthat::is.string(type),
    assertthat::noNA(type),
    isTRUE(type %in% c("FORK", "PSOCK")))
  # get environment with data
  envir <- parent.frame()
  # create cluster
  withr::with_environment(
    envir, {
    cl <- parallel::makeCluster(n, type)
    doParallel::registerDoParallel(cl)
  })
  # set up PSOCK cluster if needed
  if (identical(type, "PSOCK")) {
    ## export data to cluster
    parallel::clusterExport(cl, varlist = names, envir = envir)
  }
  # return cluster
  cl
}

#' Stop parallel processing cluster
#'
#' Stop a cluster for parallel processing.
#'
#' @param x `cluster` object.
#'
#' @return `cluster` object.
#'
#' @noRd
stop_cluster <- function(x) {
  # assert arguments are valid
  assertthat::assert_that(inherits(x, "cluster"))
  # stop cluster
  doParallel::stopImplicitCluster()
  parallel::stopCluster(x)
}
