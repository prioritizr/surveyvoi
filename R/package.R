#' @include internal.R
NULL

#' @import sf
#' @import Matrix
NULL

#' @useDynLib surveyvoi, .registration = TRUE
NULL

#' \pkg{surveyvoi}
#'
#' The \pkg{surveyvoi} package is decision support tool that will prioritize
#' sites for surveys based on their potential to improve conservation plans for
#' managing biodiversity (e.g. plans for establishing protected areas). Given a
#' set of sites that could potentially be acquired for conservation management
#' -- wherein some sites have previously been surveyed and other sites have not
#' -- this package aims to identify which sites should be surveyed because
#' doing so could lead to vastly superior conservation management plans.
#' Methods are provided to calculate the expected value of current information,
#' and the expected value of sample information given a given survey scheme.
#' Furthermore, by directly maximizing the value of sample information, optimal
#' survey schemes can also be identified.
#' This package requires the 'Gurobi' software suite (https://www.gurobi.com).
#'
#' @name surveyvoi
#' @docType package
NULL
