data {
  int<lower=0> nobs;        // number of observed Ct (< ct_lim)
  int<lower=0> ncen;        // number of censored Ct (>= ct_lim or NA)

  vector[nobs] t_obs;       // times for observed Ct
  vector[nobs] ct_obs;      // observed Ct values

  vector[ncen] t_cen;       // times for censored Ct
  real ct_lim;              // Ct right-censoring threshold (e.g., 40)
}

parameters {
  real ct0;                 // Ct at day 0
  real<lower=0> k;          // slope: Ct increases over time
  real<lower=0> sig_ct;     // Ct measurement noise
}

transformed parameters {
  vector[nobs] mu_obs;
  vector[ncen] mu_cen;

  for (i in 1:nobs) mu_obs[i] = ct0 + k * t_obs[i];
  for (j in 1:ncen) mu_cen[j] = ct0 + k * t_cen[j];
}

model {
  // Priors (reasonable defaults; adjust if your Ct range is very different)
  ct0 ~ normal(30, 10);
  k ~ normal(0.2, 0.2);        // positive, typically small
  sig_ct ~ exponential(1);

  // Likelihood for observed Ct
  ct_obs ~ normal(mu_obs, sig_ct);

  // Right-censoring: Ct >= ct_lim
  target += normal_lccdf(rep_vector(ct_lim, ncen) | mu_cen, sig_ct);
}
