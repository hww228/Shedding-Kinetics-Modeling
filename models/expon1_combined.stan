data {
    int<lower=0> nobs_ct;           // number of Ct observations
    int<lower=0> nobs_viral;        // number of viral load observations

    vector[nobs_ct] t_obs_ct;       // time points for Ct
    vector[nobs_viral] t_obs_viral; // time points for viral load

    vector[nobs_ct] c_obs_ct;       // Ct values
    vector[nobs_viral] c_obs_viral; // viral load values (log10)
}

parameters {
    real<lower=0> a0;
    real<lower=0> c0;           
    real beta0;
    real<upper=0> beta1;
    real<lower=0> sigma_viral;
}

transformed parameters {
    vector[nobs_ct] mu_ct;
    vector[nobs_viral] mu_viral;

    for (i in 1:nobs_ct) {
        mu_ct[i] = (c0 / log(10) - a0 * t_obs_ct[i] / log(10) - beta0) / beta1;
    }
    for (j in 1:nobs_viral) {
        mu_viral[j] = c0 / log(10) - a0 * t_obs_viral[j] / log(10);
    }
}

model {
    a0           ~ normal(0.05, 0.03);         
    c0           ~ lognormal(2.3, 0.2);        
    beta0        ~ normal(6.5, 1.5);           
    beta1        ~ normal(-0.1, 0.05);         
    sigma_viral  ~ exponential(1);   

    c_obs_ct ~ normal(mu_ct, sigma_viral / (-beta1));
    c_obs_viral ~ normal(mu_viral, sigma_viral);

}
