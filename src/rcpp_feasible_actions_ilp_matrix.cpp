#include "package.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_feasible_actions_ilp_matrix(Rcpp::NumericMatrix x) {
  // initialization
  const std::size_t nc = x.ncol();
  const std::size_t nr = x.nrow();
  Rcpp::NumericMatrix out(nc + 1, nr * nc);

  // fill matrix with constraints to ensure that each node is allocated
  // to a single action
  std::size_t counter = -1;
  for (std::size_t c = 0; c < nc; ++c) {
    ++counter;
    for (std::size_t r = 0; r < nr; ++r) {
      out(counter, (c * nr) + r) = 1.0;
    }
  }

  // add constraint to ensure that the total cost of the set of actions
  // does not exceed the budget
  out(nc, _) = Rcpp::as<Rcpp::NumericVector>(x);

  // return result
  return out;
}
