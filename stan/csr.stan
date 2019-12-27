// manually written CSR method


data {
  int n_origin;
  int n_dev;
  int n_data;
  real logprem[n_data];
  real logloss[n_data];
  int<lower=1,upper=n_origin> origin_lag[n_data];
  int<lower=1,upper=n_dev> dev_lag[n_data];
}

parameters {
  real r_alpha[n_origin - 1];
  real r_beta[n_dev - 1];
  real<lower=-4.0,upper=4.0> logelr;
  real<lower=0,upper=100000> a_ig[10];
  real gamma;
}
transformed parameters{
  real alpha[n_origin];
  real beta[n_dev];
  real speedup[n_origin];
  real sig2[n_dev];
  real sig[n_dev];
  real mu[n_data];
  
  alpha[1] = 0;
  for (i in 2:n_origin) alpha[i] = r_alpha[i-1];
  
  for (i in 1:(n_dev - 1)) beta[i] = r_beta[i];
  beta[n_dev] = 0;
  
  speedup[1] = 1;
  for (i in 2:n_origin) speedup[i] = speedup[i-1]*(1-gamma);
  
  sig2[n_dev] = gamma_cdf(1/a_ig[n_dev],1,1);
  for (i in 1:(n_dev - 1)) sig2[n_dev-i] = sig2[n_dev + 1 -i] + gamma_cdf(1/a_ig[i], 1, 1);
  sig = sqrt(sig2);
  
  for (i in 1:n_data){
    mu[i] = logprem[i] + logelr + alpha[origin_lag[i]] + beta[dev_lag[i]]*speedup[origin_lag[i]];
  }
}

model {
  r_alpha ~ normal(0, 3.162);
  r_beta ~ normal(0, 3.162);
  a_ig ~ inv_gamma(1,1);
  logelr ~ normal(-.4, 3.162);
  gamma ~ normal(0, 0.05);
  
  logloss ~ normal(mu, sig[dev_lag]);
}

generated quantities{
  vector[n_data] log_lik;
  vector[n_origin] sim_ay;
  real sim_total;

  for (i in 1:n_data){
    log_lik[i] = normal_lpdf(logloss[i] | mu[i], sig[dev_lag[i]]);
  }
  //log_lik = normal_lpdf(logloss | mu, sig[d]);

  for (i in 1:n_origin){
    sim_ay[i] = lognormal_rng(logprem[i] + logelr + alpha[i], sig[n_dev]);
  }
  //sim_ay = lognormal_rng(logprem + logelr + alpha, sig[n_dev]);
  sim_total = sum(sim_ay);

}

