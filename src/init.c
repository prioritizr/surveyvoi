#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _surveyvoi_rcpp_appproximate_expected_value_of_management_decision_given_current_information(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_approximate_expected_value_of_management_decision_given_perfect_information(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_approximate_expected_value_of_prioritization(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_expected_value_of_management_action(SEXP, SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_expected_value_of_management_decision_given_current_information(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_expected_value_of_management_decision_given_perfect_information(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_feasible_actions_ilp_matrix(SEXP);
extern SEXP _surveyvoi_rcpp_fit_xgboost_models_and_assess_performance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_log_sum(SEXP);
extern SEXP _surveyvoi_rcpp_n_states(SEXP);
extern SEXP _surveyvoi_rcpp_nth_state(SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_nth_state_sparse(SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_pmedian_constraint_matrix(SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_prioritization(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_probability_of_outcome(SEXP, SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_probability_of_state(SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_total_probability_of_negative_result(SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_total_probability_of_positive_result(SEXP, SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_which_state_sparse(SEXP, SEXP);
extern SEXP _surveyvoi_rcpp_xgboost(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_surveyvoi_rcpp_appproximate_expected_value_of_management_decision_given_current_information", (DL_FUNC) &_surveyvoi_rcpp_appproximate_expected_value_of_management_decision_given_current_information, 9},
    {"_surveyvoi_rcpp_approximate_expected_value_of_management_decision_given_perfect_information",  (DL_FUNC) &_surveyvoi_rcpp_approximate_expected_value_of_management_decision_given_perfect_information,  9},
    {"_surveyvoi_rcpp_approximate_expected_value_of_prioritization",                                 (DL_FUNC) &_surveyvoi_rcpp_approximate_expected_value_of_prioritization,                                 5},
    {"_surveyvoi_rcpp_expected_value_of_management_action",                                          (DL_FUNC) &_surveyvoi_rcpp_expected_value_of_management_action,                                          4},
    {"_surveyvoi_rcpp_expected_value_of_management_decision_given_current_information",              (DL_FUNC) &_surveyvoi_rcpp_expected_value_of_management_decision_given_current_information,              8},
    {"_surveyvoi_rcpp_expected_value_of_management_decision_given_perfect_information",              (DL_FUNC) &_surveyvoi_rcpp_expected_value_of_management_decision_given_perfect_information,              8},
    {"_surveyvoi_rcpp_feasible_actions_ilp_matrix",                                                  (DL_FUNC) &_surveyvoi_rcpp_feasible_actions_ilp_matrix,                                                  1},
    {"_surveyvoi_rcpp_fit_xgboost_models_and_assess_performance",                                    (DL_FUNC) &_surveyvoi_rcpp_fit_xgboost_models_and_assess_performance,                                    8},
    {"_surveyvoi_rcpp_log_sum",                                                                      (DL_FUNC) &_surveyvoi_rcpp_log_sum,                                                                      1},
    {"_surveyvoi_rcpp_n_states",                                                                     (DL_FUNC) &_surveyvoi_rcpp_n_states,                                                                     1},
    {"_surveyvoi_rcpp_nth_state",                                                                    (DL_FUNC) &_surveyvoi_rcpp_nth_state,                                                                    2},
    {"_surveyvoi_rcpp_nth_state_sparse",                                                             (DL_FUNC) &_surveyvoi_rcpp_nth_state_sparse,                                                             3},
    {"_surveyvoi_rcpp_pmedian_constraint_matrix",                                                    (DL_FUNC) &_surveyvoi_rcpp_pmedian_constraint_matrix,                                                    2},
    {"_surveyvoi_rcpp_prioritization",                                                               (DL_FUNC) &_surveyvoi_rcpp_prioritization,                                                               9},
    {"_surveyvoi_rcpp_probability_of_outcome",                                                       (DL_FUNC) &_surveyvoi_rcpp_probability_of_outcome,                                                       4},
    {"_surveyvoi_rcpp_probability_of_state",                                                         (DL_FUNC) &_surveyvoi_rcpp_probability_of_state,                                                         2},
    {"_surveyvoi_rcpp_total_probability_of_negative_result",                                         (DL_FUNC) &_surveyvoi_rcpp_total_probability_of_negative_result,                                         3},
    {"_surveyvoi_rcpp_total_probability_of_positive_result",                                         (DL_FUNC) &_surveyvoi_rcpp_total_probability_of_positive_result,                                         3},
    {"_surveyvoi_rcpp_which_state_sparse",                                                           (DL_FUNC) &_surveyvoi_rcpp_which_state_sparse,                                                           2},
    {"_surveyvoi_rcpp_xgboost",                                                                      (DL_FUNC) &_surveyvoi_rcpp_xgboost,                                                                      5},
    {NULL, NULL, 0}
};

void R_init_surveyvoi(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
