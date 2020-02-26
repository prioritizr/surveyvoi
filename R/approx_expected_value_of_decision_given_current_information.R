#' Approximate expected value of the decision given current information
#'
#' Calculate the \emph{expected value of the conservation management decision
#' given current information}. This metric describes the value of the management
#' decision that is expected when the decision maker is limited to
#' existing biodiversity data (i.e. survey data and environmental niche models).
#'
#' @param site_data \code{\link[sf]{sf}} object with site data.
#'
#' @param feature_data \code{\link[base]{data.frame}} object with feature data.
#'
#' @param site_occupancy_columns \code{character} names of \code{numeric}
#'   columns in the
#'   argument to \code{site_data} that contain presence/absence data.
#'   Each column should correspond to a different feature, and contain
#'   binary presence/absence data (zeros or ones) indicating if the
#'   feature was detected in a previous survey or not. If a site has not
#'   been surveyed before, then missing (\code{NA}) values should be used.
#'
#' @param site_probability_columns \code{character} names of \code{numeric}
#'   columns in the argument to \code{site_data} that contain modelled
#'   probabilities of occupancy for each feature in each site.
#'   Each column should correspond to a different feature, and contain
#'   probability data (values between zero and one). No missing (\code{NA})
#'   values are permitted in these columns.
#'
#' @param site_management_cost_column \code{character} name of column in the
#'   argument to \code{site_data} that contains costs for managing each
#'   site for conservation. This column should have \code{numeric} values that
#'   are equal to or greater than zero. No missing (\code{NA}) values are
#'   permitted in this column.
#'
#' @param feature_survey_sensitivity_column \code{character} name of the
#'   column in the argument to \code{feature_data} that contains
#'   probability of future surveys correctly detecting a presence of each
#'   feature in a given site (i.e. the sensitivity of the survey methodology).
#'   This column should have \code{numeric} values that are between zero and
#'   one. No missing (\code{NA}) values are permitted in this column.
#'
#' @param feature_survey_specificity_column \code{character} name of the
#'   column in the argument to \code{feature_data} that contains
#'   probability of future surveys correctly detecting an absence of each
#'   feature in a given site (i.e. the specificity of the survey methodology).
#'   This column should have \code{numeric} values that are between zero and
#'   one. No missing (\code{NA}) values are permitted in this column.
#'
#' @param feature_model_sensitivity_column \code{character} name of the
#'   column in the argument to \code{feature_data} that contains
#'   probability of the initial models correctly predicting a presence of each
#'   feature in a given site (i.e. the sensitivity of the models).
#'   This column should have \code{numeric} values that are between zero and
#'   one. No missing (\code{NA}) values are permitted in this column.
#'   This should ideally be calculated using \code{\link{fit_occupancy_models}}.
#'
#' @param feature_alpha_column \code{character} name of the column in the
#'   argument to \code{feature_data} that contains the \eqn{\alpha}
#'   values used to parametrize
#'   the conservation benefit of managing of each feature.
#'   This column should have \code{numeric} values that
#'   are equal to or greater than zero. No missing (\code{NA}) values are
#'   permitted in this column.
#'
#' @param feature_gamma_column \code{character} name of the column in the
#'   argument to \code{feature_data} that contains the \eqn{\gamma}
#'   values used to parametrize the conservation benefit of managing of each
#'   feature.
#'   This column should have \code{numeric} values that
#'   are equal to or greater than zero. No missing (\code{NA}) values are
#'   permitted in this column.
#'
#' @param total_budget \code{numeric} maximum expenditure permitted
#'   for conducting surveys and managing sites for conservation.
#'
#' @param site_management_locked_in_column \code{character} name of the column
#'   in the argument to \code{site_data} that contains \code{logical}
#'   (\code{TRUE} / \code{FALSE}) values indicating which sites should
#'   be locked in for (\code{TRUE}) being managed for conservation or
#'   (\code{FALSE}) not. No missing (\code{NA}) values are permitted in this
#'   column. This is useful if some sites have already been earmarked for
#'   conservation, or if some sites are already being managed for conservation.
#'   Defaults to \code{NULL} such that no sites are locked in.
#'
#' @param prior_matrix \code{numeric} \code{matrix} containing
#'  the prior probability of each feature occupying each site.
#'  Rows correspond to features, and columns correspond to sites.
#'  Defaults to \code{NULL} such that prior data is calculated automatically
#'  using \code{\link{prior_probability_matrix}}.
#'
#' @param n_approx_obj_fun_points \code{integer} number of points to use
#'  for approximating the piecewise-linear components of the objective
#'  function. Greater values result in more precise calculations, but
#'  also incur more computational costs. Defaults to 1000.
#'
#' @param n_approx_replicates \code{integer} number of replicates to use for
#'   approximating the expected value of a given management action.
#'   Valid arguments include an integer value greater than zero (e.g.
#'   100), or \code{NULL}. If \code{NULL} is supplied as an argument, then
#'   the expected value of a given management action is calculated precisely
#'   and not using approximation methods. Defaults to \code{NULL}.
#'
#' @param n_approx_states_per_replicate \code{integer} number of states to use
#'   per replicate for approximating the expected value of a given management
#'   action. Valid arguments include an integer value greater than zero (e.g.
#'   10000), or \code{NULL}. If \code{NULL} is supplied as an argument, then
#'   the expected value of a given management action is calculated precisely
#'   and not using approximation methods. Defaults to \code{NULL}.
#'
#' @param optimality_gap \code{numeric} relative optimality gap for generating
#'   conservation prioritizations. A value of zero indicates that
#'   prioritizations must be solved to optimality. A value of 0.1 indicates
#'   prioritizations must be within 10\% of optimality. Defaults to 0.
#'
#' @param seed \code{integer} random seed used for generating states. Defaults
#'   to 500.
#'
#' @details
#' Let \eqn{I} denote the set of feature (indexed by
#' \eqn{i}) and \eqn{J} the set of planning units (indexed by \eqn{j}).
#' To describe the results of previous surveys, let \eqn{D_{ij}} indicate the
#' detection/non-detection of feature \eqn{i \in I} in planning unit
#' \eqn{j \in J} during previous surveys using zeros and ones
#' (specified via \code{site_occupancy_columns}).
#' Also let \eqn{U_j} indicate if planning unit \eqn{j \in J} has previously
#' been surveyed or not.
#' To describe the accuracy of the surveys, let \eqn{S_i} denote the
#' sensitivity and \eqn{N_i} denote the specificity of the surveying
#' methodology for features \eqn{i \in I} (specified via
#' \code{feature_survey_sensitivity_column} and
#' \code{feature_survey_specificity_column} respectively).
#' Next, let \eqn{H_{ij}} denote the modelled probabilities of
#' the features \eqn{i \in I} occupying planning units \eqn{j \in J}
#' (specified via \code{site_probability_columns}).
#' To describe the accuracy of the environmental niche models, let \eqn{{S'}_i}
#' denote the sensitivity of the models for features \eqn{i \in I} (specified
#' via \code{feature_model_sensitivity_column}).
#' We can calculate the prior probability of each feature \eqn{i \in I}
#' occupying each site \eqn{j \in J} following
#' (or manually specified via \code{prior_matrix}):
#'
#' \deqn{
#' Q_{ij} = \\
#' S_i, \text{ if } D_j = 1, H_{ij} = 1 \space (i \text{ detected in } j) \\
#' 1 - N_i, \text{ else if } D_j = 1, H_{ij} = 0 \space (i \text{ not detected in } j) \\
#' {S'}_i {H'}_{ij}, \text{ else if } D_j = 0 \space (j \text{ not surveyed)} \\
#' }
#'
#' Since we do not know which features truly occur in which sites, let \eqn{S}
#' denote the full range of possible states that affect the evaluation of the
#' management action (indexed by \eqn{s}). Here, \eqn{S} represents all possible
#' combinations of different species occurring in different sites.
#' The total number of states \eqn{s \in S} is incredibly large for even a small
#' conservation planning problem. For example, a conservation planning problem
#' involving ten sites and four features (species) has approximately
#' \eqn{1.09 \times 10^{12}} states. In other words, one trillion ninety-nine
#' billion five hundred eleven million six hundred twenty-seven thousand seven
#' hundred seventy-five states. Unfortunately, this means that it is not
#' computationally feasible to calculate exact values of information for
#' conservation problems involving many sites. Let \eqn{G_{ijs}} indicate which
#' species \eqn{i \in I} occur in which planning units \eqn{j \in J} given
#' state \eqn{s \in S}. We calculate the prior probability of each state
#' \eqn{s \in S} being the true state following:
#'
#' \deqn{P_s = \\
#' Q_{ij}, \text{ if } G_{ijs} = 1 \\
#' 1 - Q_{ij}, \text{ else } \\
#' }
#'
#' Under the \emph{approximation method}, we define \eqn{S'} which contains a
#' subset of sites \eqn{s \in S}. The prior probability of each state
#' \eqn{s \in S'} being the true state is calculated following:
#'
#' \deqn{
#' {P'}_s = \frac{P_s}{\sum_{k \in S'} P_k}
#' }
#'
#' Here, the management action is to purchase a set of sites for
#' protected area establishment. In other words, the management action is to
#' implement a spatial prioritisation for protected area establishment. The
#' management action is subject to a budget that limits the cost of purchasing
#' sites. Let \eqn{A_j} denote the cost of acquiring each site
#' \eqn{j \in J} (specified via \code{site_management_cost_column}),
#' and \eqn{b} denote the total budget available for
#' conservation (specified via \code{total_budget}). Furthermore,
#' let \eqn{T_j} denote which sites \eqn{j \in j} are locked into the
#' prioritization (e.g. sites in existing protected area system).
#' Given these data, let \eqn{Z} denote the
#' set of all possible management actions (indexed by \eqn{z}). To describe
#' each of the possible management actions \eqn{z \in Z}, let \eqn{X_{jz}}
#' indicate which sites \eqn{j \in J} are purchased for protection (using zeros
#' and ones) under a given management action \eqn{z \in Z}. To evaluate the
#' conservation value of a given management action, we will adopt equations
#' used in the objective function of the \emph{Zonation} decision support tool
#' (Moilanen 2007). Let \eqn{\alpha_i} and \eqn{\gamma_i} denote scaling terms
#' for species \eqn{i \in I} to parametrise trade-offs in the protection of
#' different species. These scaling terms are crucial to ensure that competing
#' management decisions are evaluated in a manner conforming to the principle
#' of complementarity (Vane-Wright \emph{et al.} 1991).
#' The value of a management action (\eqn{z}) for a given state (\eqn{s}) is
#' calculated following, wherein \eqn{Y_j} indicates
#' if each site \eqn{j \in J} is selected in the (\eqn{z}) management action or
#' not (using zeros and ones):
#'
#' \deqn{
#' V(z, s) = \sum_{i \in I} \alpha_i \left( \sum_{j \in J} Y_{j} G_{ijs} \right) ^{\gamma_i} \\
#' }
#'
#' We can now calculate the
#' \emph{approximate expected value of the management decision given
#' current information} (\eqn{\text{EV'}_{\text{current}}}). To achieve this, we
#' assume that the decision maker will act optimally. Thus the
#' \emph{approximate expected value of the management decision given current
#' information} is the \emph{approximate expected value of the best management
#' action given current information}. Let \eqn{z'} denote an optimal management
#' action \eqn{z \in Z} given current information. This can be expressed as
#' follows:
#'
#' \deqn{
#' \text{EV'}_{\text{current}} = \mathbb{E} \left[ V(z', S') \right] = \max_{z \in Z} \mathbb{E} \left[ V(z, S') \right] \\
#' \mathbb{E} \left[ V(z, S') \right] = \sum_{s \in S'} V(z, s) \times {P'}_s \\
#' }
#'
#' We use combinatorial optimization to identify the (\eqn{z'})
#' best management action given current information. Specifically, we
#' formulate an integer programming problem -- using piecewise linear
#' constraints to linearise the objective function (precision controlled using
#' \code{n_approx_n_approx_obj_fun_points}) -- and solve it to (near)
#' optimality (controlled using
#' \code{optimality_gap}) with the \pkg{gurobi} package. Let \eqn{{X'}_j}
#' indicate if each planning unit \eqn{j \in J} is selected in the \eqn{z'} best
#' management action given current information (using zeros and ones):
#'
#' \deqn{
#' \text{maximize } \sum_{i \in I} \alpha_i \left( \sum_{j \in J} {X'}_j Q_{ij} \right) ^{\gamma_i} \\
#' \text{subject to } \sum_{j \in J} {{X'}_j} A_j \leq b \\
#' {{X'}_j} \geq T_j \text{ } \forall j \in J \\
#' {{X'}_j} \in \{ 0, 1 \} \text{ } \forall j \in J
#' }
#'
#' The accuracy of the approximate method depends on which subset of states
#' \eqn{s \in S'} are used for the calculations. As such, this function
#' calculates multiple estimates of the  \emph{approximate expected value of
#' the management decision given current information}
#' (\eqn{\text{EV'}_{\text{current}}}) using multiple subsets of states \eqn{s
#' \in S'} and reports the mean and standard error of these estimates. The
#' number of estimates is controlled using the \code{n_replicates} parameter,
#' and the number of states in \eqn{S'} is controlled using the
#' \code{n_states_per_replicate} parameter. For a given replicate, the states
#' are sampled randomly without replacement. This means that the
#' \emph{approximation method} is equivalent to the \emph{exact method}
#' when the argument to \code{n_states_per_replicate} is equal to the total
#' number of states.
#'
#' @references
#' Moilanen A (2007) Landscape zonation, benefit functions and target-based
#' planning: unifying reserve selection strategies.
#' \emph{Biological Conservation}, \strong{134}: 571--579.
#'
#' Vane-Wright RI, Humphries CJ, & Williams PH (1991) What to
#' protect? Systematics and the agony of choice.
#' \emph{Biological Conservation}, \strong{55}: 235--254.
#'
#' @return \code{numeric} \code{vector} containing the mean estimate and
#'   standard error.
#'
#' @seealso \code{\link{prior_probability_matrix}},
#' \code{\link{expected_value_of_decision_given_current_information}}.
#'
#' @examples
#' # set seeds for reproducibility
#' library(RandomFields)
#' set.seed(123)
#' RFoptions(seed = 123)
#'
#' # simulate data
#' site_data <- simulate_site_data(n_sites = 5, n_features = 2, prop = 0.5)
#' feature_data <- simulate_feature_data(n_features = 2, prop = 1)
#'
#' # preview simulated data
#' print(site_data)
#' print(feature_data)
#'
#' # set total budget for managing sites for conservation
#  # (i.e. 50% of the cost of managing all sites)
#' total_budget <- sum(site_data$management_cost) * 0.5
#'
#' # calculate expected value of management decision given current information
#' # using approximate method with 100 replicates and 50 states per replicate
#' ev_prime_current <-
#'   approx_expected_value_of_decision_given_current_information(
#'     site_data, feature_data, c("f1", "f2"), c("p1", "p2"),
#'     "management_cost", "survey_sensitivity",
#'     "survey_specificity", "model_sensitivity", "alpha", "gamma",
#'     total_budget, n_approx_replicates = 100,
#'     n_approx_states_per_replicate = 50)
#'
#' # print approximate value
#' print(ev_prime_current)
#'
#' @export
approx_expected_value_of_decision_given_current_information <- function(
  site_data,
  feature_data,
  site_occupancy_columns,
  site_probability_columns,
  site_management_cost_column,
  feature_survey_sensitivity_column,
  feature_survey_specificity_column,
  feature_model_sensitivity_column,
  feature_alpha_column,
  feature_gamma_column,
  total_budget,
  site_management_locked_in_column = NULL,
  prior_matrix = NULL,
  n_approx_obj_fun_points = 1000,
  optimality_gap = 0,
  n_approx_replicates = 100,
  n_approx_states_per_replicate =
    min(1000, n_states(nrow(site_data), nrow(feature_data))),
  seed = 500) {
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
    ## site_probability_columns
    is.character(site_probability_columns),
    identical(nrow(feature_data), length(site_probability_columns)),
    assertthat::noNA(site_probability_columns),
    all(assertthat::has_name(site_data, site_probability_columns)),
    ## site_management_cost_column
    assertthat::is.string(site_management_cost_column),
    all(assertthat::has_name(site_data, site_management_cost_column)),
    is.numeric(site_data[[site_management_cost_column]]),
    assertthat::noNA(site_data[[site_management_cost_column]]),
    ## feature_survey_sensitivity_column
    assertthat::is.string(feature_survey_sensitivity_column),
    all(assertthat::has_name(feature_data, feature_survey_sensitivity_column)),
    is.numeric(feature_data[[feature_survey_sensitivity_column]]),
    assertthat::noNA(
      feature_data[[feature_survey_sensitivity_column]]),
    all(feature_data[[feature_survey_sensitivity_column]] >= 0),
    all(feature_data[[feature_survey_sensitivity_column]] <= 1),
    ## feature_survey_specificity_column
    assertthat::is.string(feature_survey_specificity_column),
    all(assertthat::has_name(feature_data, feature_survey_specificity_column)),
    is.numeric(feature_data[[feature_survey_specificity_column]]),
    assertthat::noNA(feature_data[[feature_survey_specificity_column]]),
    all(feature_data[[feature_survey_specificity_column]] >= 0),
    all(feature_data[[feature_survey_specificity_column]] <= 1),
    ## feature_model_sensitivity_column
    assertthat::is.string(feature_model_sensitivity_column),
    all(assertthat::has_name(feature_data, feature_model_sensitivity_column)),
    is.numeric(feature_data[[feature_model_sensitivity_column]]),
    assertthat::noNA(feature_data[[feature_model_sensitivity_column]]),
    all(feature_data[[feature_model_sensitivity_column]] >= 0),
    all(feature_data[[feature_model_sensitivity_column]] <= 1),
    ## feature_alpha_column
    assertthat::is.string(feature_alpha_column),
    all(assertthat::has_name(feature_data, feature_alpha_column)),
    is.numeric(feature_data[[feature_alpha_column]]),
    assertthat::noNA(feature_data[[feature_alpha_column]]),
    all(feature_data[[feature_alpha_column]] >= 0),
    ## feature_gamma_column
    assertthat::is.string(feature_gamma_column),
    all(assertthat::has_name(feature_data, feature_gamma_column)),
    is.numeric(feature_data[[feature_gamma_column]]),
    assertthat::noNA(feature_data[[feature_gamma_column]]),
    all(feature_data[[feature_gamma_column]] >= 0),
    ## total_budget
    assertthat::is.number(total_budget), assertthat::noNA(total_budget),
    isTRUE(total_budget > 0),
    ## prior_matrix
    inherits(prior_matrix, c("matrix", "NULL")),
    ## n_approx_obj_fun_points
    assertthat::is.number(n_approx_obj_fun_points),
    assertthat::noNA(n_approx_obj_fun_points),
    isTRUE(n_approx_obj_fun_points > 0),
    ## n_approx_replicates
    inherits(n_approx_replicates, c("numeric", "NULL")),
    ## n_approx_states_per_replicate
    inherits(n_approx_states_per_replicate, c("numeric", "NULL")),
    identical(class(n_approx_replicates), class(n_approx_states_per_replicate)),
    ## optimality_gap
    assertthat::is.number(optimality_gap),
    assertthat::noNA(optimality_gap),
    isTRUE(optimality_gap >= 0),
    ## seed
    assertthat::is.number(seed))
  ## site_management_locked_in_column
  if (!is.null(site_management_locked_in_column)) {
    assertthat::assert_that(
      assertthat::is.string(site_management_locked_in_column),
      all(assertthat::has_name(site_data, site_management_locked_in_column)),
      is.logical(site_data[[site_management_locked_in_column]]),
      assertthat::noNA(site_data[[site_management_locked_in_column]]))
    assertthat::assert_that(
      sum(site_data[[site_management_locked_in_column]] *
          site_data[[site_management_cost_column]]) <=
      total_budget,
      msg = "cost of managing locked in sites exceeds total budget")
  }
  ## validate n_approx_states_per_replicate
  if (!is.null(n_approx_states_per_replicate)) {
    assertthat::assert_that(
      assertthat::is.count(n_approx_states_per_replicate),
      assertthat::noNA(n_approx_states_per_replicate),
      isTRUE(n_approx_states_per_replicate <=
             rcpp_n_states(nrow(site_data) * nrow(feature_data))))
  }
  ## validate n_approx_replicates
  if (!is.null(n_approx_replicates)) {
    assertthat::assert_that(
      assertthat::is.count(n_approx_replicates),
      assertthat::noNA(n_approx_replicates))
  }
  ## validate rij values
  validate_site_occupancy_data(site_data, site_occupancy_columns)
  ## validate pij values
  validate_site_prior_data(site_data, site_probability_columns)

  # drop spatial information
  if (inherits(site_data, "sf"))
    site_data <- sf::st_drop_geometry(site_data)

  # calculate prior matrix
  if (is.null(prior_matrix)) {
    pij <- prior_probability_matrix(
      site_data, feature_data, site_occupancy_columns, site_probability_columns,
      feature_survey_sensitivity_column, feature_survey_specificity_column,
      feature_model_sensitivity_column)
  } else {
    validate_prior_data(prior_matrix, nrow(site_data), nrow(feature_data))
    pij <- prior_matrix
  }

  ## prepare locked in data
  if (!is.null(site_management_locked_in_column)) {
    site_management_locked_in <- site_data[[site_management_locked_in_column]]
  } else {
    site_management_locked_in <- rep(FALSE, nrow(site_data))
  }

  # set the seed
  rng_state <- .Random.seed
  set.seed(seed)

  # main calculation
  out <- rcpp_approx_expected_value_of_decision_given_current_info_n_states(
    pij = pij,
    pu_costs = site_data[[site_management_cost_column]],
    pu_locked_in = site_management_locked_in,
    alpha = feature_data[[feature_alpha_column]],
    gamma = feature_data[[feature_gamma_column]],
    n_approx_obj_fun_points = n_approx_obj_fun_points,
    budget = total_budget,
    gap = optimality_gap,
    n_approx_replicates = n_approx_replicates,
    n_approx_states_per_replicate = n_approx_states_per_replicate)

  # restore the previous state
  set.seed(rng_state)

  # return result
  names(out) <- c("mean", "se")
  out
}
