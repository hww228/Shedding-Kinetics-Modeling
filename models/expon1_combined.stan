data {
  // --- Viral part (same structure as expon1cen.stan) ---
  int<lower=0> nobs;
  int<lower=0> ncen;
  vector[nobs] t_obs;
  vector[nobs] c_obs;          // observed log10 VL
  vector[ncen] t_cen;          // censored VL times
  real<upper=min(c_obs)> censlim;

  // --- Ct part ---
  int<lower=0> nobs_ct;
  int<lower=0> ncen_ct;
  vector[nobs_ct] t_obs_ct;
  vector[nobs_ct] ct_obs;
  vector[ncen_ct] t_cen_ct;
  real ct_lim;                  // Ct censoring threshold (e.g., 40)
}

parameters {
  // --- Viral parameters (same as expon1cen.stan) ---
  real<lower=0> sig_obs;
  real a0;
  real<lower=0> c0;

  // --- Ct forward calibration + noise ---
  real alpha0;
  real<upper=0> alpha1;         // negative: higher VL -> lower Ct
  real<lower=0> sig_ct;
}

transformed parameters {
  // Viral means
  vector[nobs] gene_obs;
  vector[ncen] gene_cen;

  // Ct means
  vector[nobs_ct] mu_ct_obs;
  vector[ncen_ct] mu_ct_cen;

  for (i in 1:nobs)
    gene_obs[i] = c0 / log(10) - a0 * t_obs[i] / log(10);

  for (j in 1:ncen)
    gene_cen[j] = c0 / log(10) - a0 * t_cen[j] / log(10);

  for (i in 1:nobs_ct) {
    real mu_vl = c0 / log(10) - a0 * t_obs_ct[i] / log(10);
    mu_ct_obs[i] = alpha0 + alpha1 * mu_vl;
  }

  for (j in 1:ncen_ct) {
    real mu_vl = c0 / log(10) - a0 * t_cen_ct[j] / log(10);
    mu_ct_cen[j] = alpha0 + alpha1 * mu_vl;
  }
}

model {
  // --- Priors (keep your original viral priors) ---
  sig_obs ~ uniform(0, 10);
  a0 ~ normal(0, 100);
  c0 ~ normal(0, 100);

  // --- Ct priors (typical scale: Ct decreases ~3 cycles per +1 log10 VL) ---
  alpha0 ~ normal(40, 10);
  alpha1 ~ normal(-3, 1);
  sig_ct  ~ exponential(1);

  // --- Viral likelihood (same as expon1cen.stan) ---
  c_obs ~ normal(gene_obs, sig_obs);
  target += normal_lcdf(censlim | gene_cen, sig_obs);

  // --- Ct likelihood ---
  ct_obs ~ normal(mu_ct_obs, sig_ct);
  target += normal_lccdf(rep_vector(ct_lim, ncen_ct) | mu_ct_cen, sig_ct);
}

