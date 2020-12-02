#' @include internal.R
NULL

#' @import sf
#' @import Matrix
#' @import nloptr
#' @importFrom Rcpp sourceCpp
NULL

#' @useDynLib surveyvoi, .registration = TRUE
NULL

#' surveyvoi: Survey Value of Information
#'
#' The \pkg{surveyvoi} package is a decision support tool for prioritizing
#' sites for ecological surveys based on their potential to improve
#' conservation plans for managing biodiversity (e.g. plans for establishing
#' protected areas). Given a set of sites that could potentially be acquired
#' for conservation management -- wherein some sites have previously been
#' surveyed and other sites have not -- this package provides functionality to
#' identify which sites should be surveyed because doing so would likely lead
#' to superior conservation management plans. It uses value of information
#' analyses to measure return on investment for competing survey schemes.
#' Furthermore, by directly maximizing the expected value of sample
#' information, optimal survey schemes can also be identified.
#'
#' @name surveyvoi
#' @docType package
NULL
