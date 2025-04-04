data {
  int<lower=0> nobs;// number of observation
  vector[nobs] t;   // noncensored observations time vector
  vector[nobs] c;   // noncensored observations log10 concentration vector
}

parameters {
  real<lower=0> sig_obs;
  real a0; //rate;
  real b0; //shape;
  real<lower=0> c0; //log scale concentration at day 0;
}

transformed parameters {
  vector[nobs] gene_obs;
  for (k_obs in 1:nobs)
    gene_obs[k_obs] = c0/log(10)+log(t[k_obs]^b0)/log(10)-a0*t[k_obs]/log(10);
}

model {
  // priors
  sig_obs ~ uniform(0,10);
  a0 ~ normal(0,100);
  b0 ~ normal(0,100);
  c0 ~ normal(0,100);
  // likelihood
  c ~ normal(gene_obs, sig_obs);
}