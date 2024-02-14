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
#' @description
#' Decision support tool for prioritizing sites for ecological
#' surveys based on their potential to improve plans for conserving
#' biodiversity (e.g. plans for establishing protected areas). Given a set of
#' sites that could potentially be acquired for conservation management --
#' wherein some sites have previously been surveyed and other sites have not
#' -- it can be used to generate and evaluate plans for
#' additional surveys. Specifically, plans for ecological surveys can be
#' generated using various conventional approaches (e.g. maximizing expected
#' species richness, geographic coverage, diversity of sampled environmental
#' conditions) and by maximizing value of information. After generating
#' plans for surveys, they can also be evaluated using
#' value of information analysis.
#'
#' @details
#' Please note that several functions depend on
#' the 'Gurobi' optimization software (available from <https://www.gurobi.com>)
#' and the \pkg{gurobi} R package (installation instructions
#' [available online for Linux, Windows, and Mac OS](https://support.gurobi.com/hc/en-us/articles/4534161999889-How-do-I-install-Gurobi-Optimizer)).
#' Additionally, the JAGS software
#' (available from <https://mcmc-jags.sourceforge.io/>) is required to fit
#' hierarchical generalized linear models.
#'
#' @seealso
#' The package vignette provides a tutorial
#' (accessible using the code `vignettes("surveyvoi")`).
#'
#' @name surveyvoi
#'
#' @aliases surveyvoi-package
#'
#' @author
#' Package authors:
#' * Jeffrey O. Hanson \email{jeffrey.hanson@uqconnect.edu.au} [ORCID](https://orcid.org/0000-0002-4716-6134)
#' * Iadine Chadès \email{iadine.chades@csiro.au} [ORCID](https://orcid.org/0000-0002-7442-2850)
#' * Emma J. Hudgins \email{emma.hudgins@mail.mcgill.ca} [ORCID](https://orcid.org/0000-0002-8402-5111)
#' * Joseph R. Bennett \email{joseph.bennett@carleton.ca} [ORCID](https://orcid.org/0000-0002-3901-9513)
"_PACKAGE"

# ensure package checks pass
#' @importFrom utils zip
NULL
