#' Expected value of the decision given current information
#'
#' Calculate the \emph{expected value of the conservation management decision
#' given current information}. This metric describes the value of the management
#' decision that is expected when the decision maker is limited to
#' existing biodiversity data (i.e. survey data and environmental niche models).
#'
#' @inheritParams approx_evdci
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
#' and \eqn{{N'}_i} denote the sensitivity and specificity (respectively) of
#' the models for features \eqn{i \in I} (specified
#' via \code{feature_model_sensitivity_column} and
#' \code{feature_model_specificity_column} respectively).
#' We can calculate the prior probability of each feature \eqn{i \in I}
#' occupying each site \eqn{j \in J} following
#' (or manually specified via \code{prior_matrix}):
#'
#' \deqn{
#' P_{ij} = \\
#' S_i, \text{ if } D_j = 1, H_{ij} = 1 \space (i \text{ detected in } j) \\
#' 1 - N_i, \text{ else if } D_j = 1, H_{ij} = 0 \space (i \text{ not detected in } j) \\
#' {S'}_i, \text{ else if } D_j = 0, {H'}_{ij} \geq 0.5 \space (j \text{ not surveyed and } i \text{ predicted present in } j \text{)} \\
#' 1 - {N'}_i, \text{ else if } D_j = 0, {H'}_{ij} \geq 0.5 \space (j \text{ not surveyed and } i \text{ predicted absent in } j \text{)} \\
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
#' \emph{expected value of the management decision given
#' current information} (\eqn{\text{EV}_{\text{current}}}). To achieve this, we
#' assume that the decision maker will act optimally. Thus the
#' \emph{expected value of the management decision given current information}
#' is the \emph{expected value of the best management action given current
#' information}. Let \eqn{z'} denote an optimal management action
#' \eqn{z \in Z} given current information. This can be expressed as follows:
#'
#' \deqn{
#' \text{EV}_{\text{current}} = \mathbb{E} \left[ V(z', S) \right] = \max_{z \in Z} \mathbb{E} \left[ V(z, S) \right] \\
#' \mathbb{E} \left[ V(z, S) \right] = \sum_{s \in S} V(z, s) \times P_s \\
#' }
#'
#' Typically, the emph{expected value of the best management action given
#' current information} is calculated by iterating over every possible
#' management action \eqn{z \in Z}. However, this approach would be
#' computationally expensive for even a small conservation planning problem.
#' For example, a conservation planning problem involving 60 planning units has
#' approximately \eqn{1.52 \times 10^{18}} different solutions (i.e.
#' combinations of different planning units being selected). Therefore, we
#' applied combinatorial optimization instead. Specifically, we identify the
#' best management action by formulating an
#' integer programming problem -- using piecewise linear constraints to
#' linearise the objective function (precision controlled using
#' \code{n_approx_n_approx_obj_fun_points}) -- and solving it to (near)
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
#' @references
#' Moilanen A (2007) Landscape zonation, benefit functions and target-based
#' planning: unifying reserve selection strategies.
#' \emph{Biological Conservation}, \strong{134}: 571--579.
#'
#' Vane-Wright RI, Humphries CJ, & Williams PH (1991) What to
#' protect? Systematics and the agony of choice.
#' \emph{Biological Conservation}, \strong{55}: 235--254.
#'
#' @return \code{numeric} value.
#'
#' @seealso \code{\link{prior_probability_matrix}},
#' \code{\link{approx_evdci}}.
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
#' # using exact method
#' ev_current <- evdci(
#'   site_data, feature_data, c("f1", "f2"), c("p1", "p2"),
#'   "management_cost", "survey_sensitivity", "survey_specificity",
#'   "model_sensitivity", "model_specificity", "alpha", "gamma", total_budget)
#'
#' # print exact value
#' print(ev_current)
#'
#' @export
evdci <- function(
  site_data,
  feature_data,
  site_occupancy_columns,
  site_probability_columns,
  site_management_cost_column,
  feature_survey_sensitivity_column,
  feature_survey_specificity_column,
  feature_model_sensitivity_column,
  feature_model_specificity_column,
  feature_alpha_column,
  feature_gamma_column,
  total_budget,
  site_management_locked_in_column = NULL,
  prior_matrix = NULL,
  n_approx_obj_fun_points = 1000,
  optimality_gap = 0) {
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
    ## feature_model_specificity_column
    assertthat::is.string(feature_model_specificity_column),
    all(assertthat::has_name(feature_data, feature_model_specificity_column)),
    is.numeric(feature_data[[feature_model_specificity_column]]),
    assertthat::noNA(feature_data[[feature_model_specificity_column]]),
    all(feature_data[[feature_model_specificity_column]] >= 0),
    all(feature_data[[feature_model_specificity_column]] <= 1),
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
    ## optimality_gap
    assertthat::is.number(optimality_gap),
    assertthat::noNA(optimality_gap),
    isTRUE(optimality_gap >= 0))
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
      feature_model_sensitivity_column, feature_survey_specificity_column)
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

  # main calculation
  out <- rcpp_expected_value_of_decision_given_current_info(
    pij = pij,
    pu_costs = site_data[[site_management_cost_column]],
    pu_locked_in = site_management_locked_in,
    alpha = feature_data[[feature_alpha_column]],
    gamma = feature_data[[feature_gamma_column]],
    n_approx_obj_fun_points = n_approx_obj_fun_points,
    budget = total_budget,
    gap = optimality_gap)

  # return result
  out
}
