# include "rcpp_model_performance.h"

void model_sensitivity_and_specificity(
  Eigen::VectorXf &y, Eigen::VectorXf &yhat, Eigen::VectorXf &w,
  double data_sens, double data_spec,
  double &model_sens, double &model_spec) {

  // check if preds contains at least one predicted absence and presence,
  // if not return very poor performance
  double max_wt_yhat_pres =
    static_cast<double>((w.array() *
                         (yhat.array() >= 0.5).cast<float>()).maxCoeff());
  double max_wt_yhat_abs =
    static_cast<double>((w.array() *
                         (yhat.array() < 0.5).cast<float>()).maxCoeff());
  if ((max_wt_yhat_pres < 1.0e-5) || (max_wt_yhat_abs < 1.0e-5)) {
    model_sens = 0.0;
    model_spec = 0.0;
    return;
  }

  // generate confusion table
  double total_positive = static_cast<double>((y.array() * w.array()).sum());
  double total_negative = static_cast<double>(w.sum()) - total_positive;
  double true_positive = static_cast<double>(
    (w.array() *
    ((yhat.array() >= 0.5) && (y.array() >= 0.5)).cast<float>()).sum()
  );
  double true_negative = static_cast<double>(
    (w.array() *
     ((yhat.array() < 0.5) && (y.array() < 0.5)).cast<float>()).sum()
  );
  double false_negative = total_positive - true_positive;
  double false_positive = total_negative - true_negative;

  // estimate model performance using Staquet formulas
  formula_sensitivity_and_specificity(
    true_positive, false_positive, false_negative, true_negative,
    data_sens, data_spec, model_sens, model_spec);

  // if the formula approach gives invalid estimates, then
  // re-estimate model performance using maximum likelihood approach,
  if ((model_sens < 0.0) || (model_sens > 1.0) ||
      (model_spec < 0.0) || (model_spec > 1.0)) {
  maxlik_sensitivity_and_specificity(
    true_positive, false_positive, false_negative, true_negative,
    data_sens, data_spec, model_sens, model_spec);
  }

  // clamp values to (1e-10) and (1 - 1e-10) to avoid numerical issues
  // with probabilities that are exactly zero and one
  model_sens = std::max(model_sens, 1.0e-10);
  model_sens = std::min(model_sens, 1.0 - 1.0e-10);
  model_spec = std::max(model_spec, 1.0e-10);
  model_spec = std::min(model_spec, 1.0 - 1.0e-10);

  // return void
  return;
}

void formula_sensitivity_and_specificity(
  double true_positive, double false_positive,
  double false_negative, double true_negative,
  double data_sens, double data_spec,
  double &model_sens, double &model_spec) {
  // calculate model sensitivity and specificity
  double n = true_positive + false_positive + false_negative + true_negative;
  model_sens =
    (((true_positive + false_positive) * data_spec) - false_positive) /
    ((n * (data_spec - 1.0)) + (true_positive + false_negative));
  model_spec =
    (((true_negative + false_negative) * data_sens) - false_negative) /
    ((n * data_sens) - (true_positive + false_negative));
  // return void
  return;
}

void maxlik_sensitivity_and_specificity(
  double true_positive, double false_positive,
  double false_negative, double true_negative,
  double data_sens, double data_spec,
  double &model_sens, double &model_spec) {

  // prepare values for optimization
  std::vector<double> par = {0.9, 0.9, 0.5};
  std::vector<double> lb(3, 1.0e-10);
  std::vector<double> ub(3, 1.0 - 1.0e-10);
  double value;

  // prepare data for optimization
  nll_data f_data;
  mpfr_inits2(
    1000, f_data.tp, f_data.fn, f_data.fp, f_data.tn,
    f_data.t1, f_data.t2, f_data.t3, f_data.t4, (mpfr_ptr) 0);
  mpfr_set_d(f_data.tp, true_positive, MPFR_RNDD);
  mpfr_set_d(f_data.fn, false_negative, MPFR_RNDD);
  mpfr_set_d(f_data.fp, false_positive, MPFR_RNDD);
  mpfr_set_d(f_data.tn, true_negative, MPFR_RNDD);
  f_data.dse = data_sens,
  f_data.dsp = data_spec;

  // run optimization
  nlopt::opt opt(nlopt::LN_BOBYQA, static_cast<unsigned>(3));
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_min_objective(nll, &f_data);
  opt.set_xtol_rel(1.0e-5);
  try {
    nlopt::result res = opt.optimize(par, value);
  } catch (std::exception &e) {

  }

  // clean up
  mpfr_clears(
    f_data.t1, f_data.t2, f_data.t3, f_data.t4,
    f_data.tp, f_data.fn, f_data.fp, f_data.tn, (mpfr_ptr) 0);
  mpfr_free_cache();

  // export result
  model_sens = par[0];
  model_spec = par[1];

  // return void;
  return;
}

double nll (
  const std::vector<double> &x, std::vector<double> &grad, void *f_data) {
  // calculate negative log likelihood
  nll_data *d = reinterpret_cast<nll_data*>(f_data);
  mpfr_set_d(
    d->t1,
    (d->dse * x[0] * x[2]) +
    ((1.0 - d->dsp) * (1.0 - x[1]) * (1.0 - x[2])),
   MPFR_RNDD);
  mpfr_set_d(
    d->t2,
    (d->dse * (1.0 - x[0]) * x[2]) +
    ((1.0 - d->dsp) * x[1] * (1.0 - x[2])),
    MPFR_RNDD);
  mpfr_set_d(
    d->t3,
    ((1.0 - d->dse) * x[0] * x[2]) +
    (d->dsp * (1.0 - x[1]) * (1.0 - x[2])),
    MPFR_RNDD);
  mpfr_set_d(
    d->t4,
    ((1.0 - d->dse) * (1.0 - x[0]) * x[2]) +
    (d->dsp * x[1] * (1.0 - x[2])),
    MPFR_RNDD);
  mpfr_pow(d->t1, d->t1, d->tp, MPFR_RNDD);
  mpfr_pow(d->t2, d->t2, d->fn, MPFR_RNDD);
  mpfr_pow(d->t3, d->t3, d->fp, MPFR_RNDD);
  mpfr_pow(d->t4, d->t4, d->tn, MPFR_RNDD);
  mpfr_log(d->t1, d->t1, MPFR_RNDD);
  mpfr_log(d->t2, d->t2, MPFR_RNDD);
  mpfr_log(d->t3, d->t3, MPFR_RNDD);
  mpfr_log(d->t4, d->t4, MPFR_RNDD);
  // return result
  return -1.0 * (mpfr_get_d(d->t1, MPFR_RNDD) +
                 mpfr_get_d(d->t2, MPFR_RNDD) +
                 mpfr_get_d(d->t3, MPFR_RNDD) +
                 mpfr_get_d(d->t4, MPFR_RNDD));
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_model_performance(
  Eigen::VectorXd y, Eigen::VectorXd yhat, Eigen::VectorXd w,
  double data_sens, double data_spec) {
  // prepare data for calculations
  double model_sens, model_spec;
  Eigen::VectorXf yf = y.cast<float>();
  Eigen::VectorXf yhatf = yhat.cast<float>();
  Eigen::VectorXf wf = w.cast<float>();
  // main calculations
  model_sensitivity_and_specificity(
    yf, yhatf, wf, data_sens, data_spec, model_sens, model_spec);
  // return result
  Rcpp::NumericVector out(3);
  out[0] = model_sens + model_spec - 1.0;
  out[1] = model_sens;
  out[2] = model_spec;
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_formula_sensitivity_and_specificity(
  Rcpp::NumericMatrix x, double data_sens, double data_spec) {
  // prepare data for calculations
  double model_sens, model_spec;
  // main calculations
  formula_sensitivity_and_specificity(
    x(0, 0), x(0, 1), x(1, 0), x(1, 1),
    data_sens, data_spec,
    model_sens, model_spec);
  // return result
  Rcpp::NumericVector out(2);
  out[0] = model_sens;
  out[1] = model_spec;
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_maxlik_sensitivity_and_specificity(
  Rcpp::NumericMatrix x, double data_sens, double data_spec) {
  // prepare data for calculations
  double model_sens, model_spec;
  // main calculations
  maxlik_sensitivity_and_specificity(
    x(0, 0), x(0, 1), x(1, 0), x(1, 1),
    data_sens, data_spec,
    model_sens, model_spec);
  // return result
  Rcpp::NumericVector out(2);
  out[0] = model_sens;
  out[1] = model_spec;
  return out;
}
