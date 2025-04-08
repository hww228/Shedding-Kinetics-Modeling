data {
  int<lower=0> nobs;         
  vector[nobs] t;            
  vector[nobs] ct_value;     
}

parameters {
  real<lower=0> sigma;     
  real a0;                   
  real<lower=0> c0;          
  real<lower=0> beta0;                
  real<upper=0> beta1;                
}

transformed parameters {
  vector[nobs] mu_obs;     
  for (i in 1:nobs)
    mu_obs[i] = (c0 / log(10) - a0 * t[i] / log(10) - beta0) / beta1;
}

model {
  sigma ~ normal(1, 0.5);
  a0 ~ normal(0.5, 0.5);
  c0 ~ lognormal(2, 0.5);
  beta0 ~ normal(2, 1);
  beta1 ~ normal(-0.3, 0.05);
  
  ct_value ~ normal(mu_obs, sigma/(-beta1)); #observation distribution
}
// model {
//   sigma ~ uniform(0, 10);
//   a0 ~ normal(0, 100);
//   c0 ~ normal(0, 100);
//   beta0 ~ uniform(0, 10);
//   beta1 ~ normal(0, 10);
//   
//   ct_value ~ normal(mu_obs, -sigma/beta1);
// }


