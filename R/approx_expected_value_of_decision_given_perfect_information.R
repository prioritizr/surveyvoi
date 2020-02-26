#' Approximate expected value of decision given perfect information
#'
#' Calculate the \emph{expected value of the conservation management decision
#' given perfect information}. This metric describes the value of the management
#' decision that is expected when the decision maker is omniscient (i.e. they
#' know exactly which features occur in which sites). Although real-world
#' conservation planing scenarios rarely involve such omniscient decision
#' makers, this metric is useful to provide an upper bound on the expected
#' value of management decisions following additional data collection.
#'
#' @inheritParams approx_expected_value_of_decision_given_current_information
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
#' We can now calculate the \emph{approximate expected value of the management
#' decision given perfect information} (\eqn{\text{EV'}_{\text{certainty}}}).
#' To achieve this, we assume that the decision maker will act optimally based
#' on their perfect knowledge. Thus the \emph{approximate expected value of the
#' management decision given perfect information} is the \emph{approximate
#' expected value of the best management action given perfect information}. Let
#' \eqn{z''} denote an optimal
#' management action \eqn{z \in Z} given perfect information. This can be
#' expressed as follows:
#'
#' \deqn{
#' \text{EV'}_{\text{certainty}} = \mathbb{E} \left[ \max_{z \in Z} V(z, S') \right] = \sum_{s \in S'} \left[ \left\{ V({{z''}_s}, s) \right\} \times {P'}_s \right]
#' }
#'
#' We use combinatorial optimization to identify the (\eqn{{z''}_s})
#' best management action for each state \eqn{s \in S'}. Specifically, we
#' formulate an integer programming problem -- using piecewise linear
#' constraints to linearise the objective function (precision controlled using
#' \code{n_approx_n_approx_obj_fun_points}) -- and solve it to (near)
#' optimality (controlled using
#' \code{optimality_gap}) with the \pkg{gurobi} package. Let \eqn{{X''}_{js}}
#' indicate if each planning unit \eqn{j \in J} is selected in the
#' \eqn{{z''}_s} best management action for state eqn{s} (using zeros and ones):
#'
#' \deqn{
#' \text{maximize } \sum_{i \in I} \alpha_i \left( \sum_{j \in J} {X''}_{js}
#'  G_{ijs} \right) ^{\gamma_i} \\
#' \text{subject to } \sum_{j \in J} {{X''}_{js}} A_j \leq b \\
#' {{X''}_{js}} \geq T_j \text{ } \forall j \in J \\
#' {{X''}_{js}} \in \{ 0, 1 \} \text{ } \forall j \in J
#' }
#'
#' The accuracy of the approximate method depends on which subset of states
#' \eqn{s \in S'} are used for the calculations. As such, this function
#' calculates multiple estimates of the  \emph{approximate expected value of
#' the management decision given perfect information}
#' (\eqn{\text{EV'}_{\text{certainty}}}) using multiple subsets of states \eqn{s
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
#' @inherit approx_expected_value_of_decision_given_current_information return
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
#' # calculate expected value of management decision given perfect information
#' # using approximate method with 100 replicates and 50 states per replicate
#' ev_prime_certainty <-
#'   approx_expected_value_of_decision_given_perfect_information(
#'     site_data, feature_data, c("f1", "f2"), c("p1", "p2"),
#'     "management_cost", "survey_sensitivity",
#'     "survey_specificity", "model_sensitivity", "alpha", "gamma",
#'     total_budget, n_approx_replicates = 100,
#'     n_approx_states_per_replicate = 50)
#'
#' # print approximate value
#' print(ev_prime_certainty)
#'
#' @seealso \code{\link{prior_probability_matrix}},
#' \code{\link{expected_value_of_decision_given_perfect_information}}.
#'
#' @export
approx_expected_value_of_decision_given_perfect_information <- function(
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
  out <- rcpp_approx_expected_value_of_decision_given_perfect_info_n_states(
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
