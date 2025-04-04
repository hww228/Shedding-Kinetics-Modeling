data {
  int<lower=0> nobs;    // number of noncensored observation
  int<lower=0> ncen;    // number of censored observation
  vector[nobs] t_obs;   // noncensored observations time vector
  vector[nobs] c_obs;   // noncensored observations log10 concentration vector
  vector[ncen] t_cen;   // censored observations time vector
  real<upper=min(c_obs)> censlim; //log10 scale limit of quantification
}
parameters {
  real<lower=0> sig_obs;
  real a0; //decay rate;
  real<lower=0> c0; //log scale concentration at day 0;
}

transformed parameters {
  vector[nobs] gene_obs;
  vector[ncen] gene_cen;
  for (k_obs in 1:nobs)
    gene_obs[k_obs] = c0/log(10)-a0*t_obs[k_obs]/log(10);
  for (k_cen in 1:ncen)
    gene_cen[k_cen] = c0/log(10)-a0*t_cen[k_cen]/log(10);
}

model {
  // priors
  sig_obs ~ uniform(0,10);
  a0 ~ normal(0,100);
  c0 ~ normal(0,100);
  // likelihood
  c_obs ~ normal(gene_obs, sig_obs);
  target += normal_lcdf(censlim | gene_cen, sig_obs);
}