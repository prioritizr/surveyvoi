skip_if_parallel_processing_not_available <- function() {
  skip_if(
    (!requireNamespace("surveyvoi", quietly = TRUE) &&
      !identical(.Platform$OS.type, "unix")) ||
    (!"surveyvoi" %in% rownames(installed.packages())),
    message = "parallel processing not available"
  )
}

skip_on_fedora <- function() {
  os_name <- utils::sessionInfo()$running
  is_fedora <- TRUE
  if (
    is.character(os_name) &&
    identical(length(os_name), 1L) &&
    all(!is.na(os_name))
  ) {
    is_fedora <- any(grepl("fedora", tolower(os_name), fixed = TRUE))
  }
  skip_if(isTRUE(is_fedora), message = "on Fedora")
}
