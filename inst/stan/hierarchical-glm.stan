data {
  // constants
  int<lower=1> n_vars;
  int<lower=1> train_n_sites;
  int<lower=1> predict_n_sites;

  // feature data
  real sens;
  real spec;

  // training variables
  real[train_n_sites, n_vars] train_model_matrix;
  vector[train_n_sites] train_n_det;
  vector[train_n_sites] train_n_non_det;

  // prediction variables
  real[predict_n_sites, n_vars] predict_model_matrix;
}

transformed data {
  real sens_log = log(sens);
  real spec_log = log(spec);
}

parameters {
  vector[n_vars] train_coef;
  vector[train_n_sites] train_occ;
}

transformed parameters {
  predict_model_matrix

}


model {
  // priors
  train_coef ~ normal(0, 5);

  // likelihood
  y ~ normal(intercept + beta * x, sigma);
}
