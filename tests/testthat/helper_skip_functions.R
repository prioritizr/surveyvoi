skip_if_parallel_processing_not_available <- function() {
  skip_if(
    (!requireNamespace("surveyvoi", quietly = TRUE) &&
      !identical(.Platform$OS.type, "unix")) ||
    (!"surveyvoi" %in% rownames(installed.packages())),
    message = "parallel processing not available"
  )
}
