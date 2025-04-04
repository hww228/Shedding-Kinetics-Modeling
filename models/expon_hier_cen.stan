data {
  int<lower=0> nsubj;   // number of subjects
  int<lower=0> nobs;    // number of noncensored observation
  int<lower=0> ncen;    // number of censored observation
  vector[nobs] t_obs;   // noncensored observations time vector
  vector[nobs] c_obs;   // noncensored observations log10 concentration vector
  int p_obs[nobs];      // noncensored observations subject vector
  vector[ncen] t_cen;   // censored observations time vector
  int p_cen[ncen];      // censored observations subject vector
  real<upper=min(c_obs)> censlim;
}

parameters {
  real<lower=0> sig_obs; // prediction error scale;
  cholesky_factor_corr[2] L_Omega; // Cholesky factorization for LKJ prior correlation;
  vector<lower=0>[2] tau;
  vector[2] mu;
  matrix[nsubj,2] par;
}

transformed parameters {
  vector[nsubj] a0; //decay rate;
  vector<lower=0>[nsubj] c0; //log scale concentration at day 0;
  for (k_subj in 1:nsubj){
    a0[k_subj] = exp(par[k_subj,1]);
    c0[k_subj] = exp(par[k_subj,2]);
  }
}

model {
  // priors
  sig_obs ~ uniform(0,10);
  tau ~ cauchy(0, 2.5);
  L_Omega ~ lkj_corr_cholesky(2);
  mu ~ normal(0,10);
  for (k_subj in 1:nsubj)
    //Using Cholesky parameterization;
    par[k_subj,] ~ multi_normal_cholesky(mu,diag_pre_multiply(tau, L_Omega));
  vector[nobs] gene_obs;
  vector[ncen] gene_cen;
  for (k_obs in 1:nobs)
    gene_obs[k_obs] = c0[p_obs[k_obs]]/log(10)-a0[p_obs[k_obs]]*t_obs[k_obs]/log(10);
  for (k_cen in 1:ncen)
    gene_cen[k_cen] = c0[p_cen[k_cen]]/log(10)-a0[p_cen[k_cen]]*t_cen[k_cen]/log(10);
  // likelihood
  c_obs ~ normal(gene_obs, sig_obs);
  target += normal_lcdf(censlim | gene_cen, sig_obs);
}