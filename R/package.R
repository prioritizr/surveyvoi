#' @include internal.R
NULL

#' @import sf
#' @import Matrix
#' @importFrom Rcpp sourceCpp
#' @import nloptr
#' @useDynLib surveyvoi, .registration = TRUE
NULL

#' surveyvoi: Survey Value of Information
#'
#' Decision support tool for prioritizing sites for ecological
#' surveys based on their potential to improve plans for conserving
#' biodiversity (e.g. plans for establishing protected areas). Given a set of
#' sites that could potentially be acquired for conservation management --
#' wherein some sites have previously been surveyed and other sites have not
#' -- it can be used to generate and evaluate plans for
#' additional surveys. Specifically, plans for ecological surveys can be
#' generated using various conventional approaches (e.g. maximizing expected
#' species richness, geographic coverage, diversity of sampled environmental
#' conditions) and directly maximizing return on investment via optimization
#' algorithms. After generating such survey plans, they can be evaluated using
#' value of information analysis. Please note that several functions depend on
#' the 'Gurobi' optimization software (available from <https://www.gurobi.com>)
#' and the \pkg{gurobi} R package (installation instructions available for
#' [Linux](https://www.gurobi.com/documentation/9.1/quickstart_linux/r_ins_the_r_package.html), [Windows](https://www.gurobi.com/documentation/9.1/quickstart_windows/r_ins_the_r_package.html), and [Mac OS](https://www.gurobi.com/documentation/9.1/quickstart_mac/r_ins_the_r_package.html)).
#' See the package vignette for a tutorial
#' (accessible using the code `vignettes("surveyvoi")`).
#'
#' @name surveyvoi
#' @docType package
NULL
