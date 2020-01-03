// manually written CSR method


data {
  int n_origin;
  int n_dev;
  int n_data;
  real logprem[n_data];
  real logloss[n_data];
  int<lower=1,upper=n_origin> origin_lag[n_data];
  int<lower=1,upper=n_dev> dev_lag[n_data];
  int prior_only;  // should the likelihood be ignored?
  real logelrparams[2];
  real alphasigma;
  vector[2] betaparams[n_dev];
  real sigparams;
  real gammasig;
}

parameters {
  real<lower=-5.0,upper=5.0> r_alpha[n_origin - 1];
  real<lower=-5.0,upper=8.0> r_beta[n_dev - 1];
  real<lower=-4.0,upper=4.0> logelr;
  real<lower=0,upper=100000> a_ig[n_dev];
  real<lower=-3.0,upper=3.0> gamma;
}
transformed parameters{
  real alpha[n_origin];
  real beta[n_dev];
  real speedup[n_origin];
  real sig[n_dev];
  real mu[n_data];
  
  alpha[1] = 0;
  alpha[2:n_origin] = r_alpha[1:(n_origin - 1)];
  
  beta[1:(n_dev - 1)] = r_beta[1:(n_dev - 1)];
  beta[n_dev] = betaparams[n_dev, 1];
  
  speedup[1] = 1;
  for (i in 2:n_origin) speedup[i] = speedup[i-1] * (1 - gamma);
  
  sig[n_dev] = gamma_cdf(1/a_ig[n_dev], sigparams , sigparams);
  for (i in 1:(n_dev - 1)) sig[n_dev-i] = sig[n_dev + 1 - i] + gamma_cdf(1 / a_ig[i], sigparams, sigparams);
  sig = sqrt(sig);
  
  for (i in 1:n_data){
    mu[i] = logprem[i] + logelr + alpha[origin_lag[i]] + beta[dev_lag[i]] * speedup[origin_lag[i]];
  }
}

model {
  target += student_t_lpdf(r_alpha | 3, 0, alphasigma);
  target += student_t_lpdf(r_beta | 3, betaparams[1:(n_dev-1), 1], betaparams[1:(n_dev-1), 2]);
  target += inv_gamma_lpdf(a_ig | sigparams , sigparams);
  target += normal_lpdf(logelr | logelrparams[1], logelrparams[2]);
  target += student_t_lpdf(gamma | 4, 0, gammasig);
  
  // likelihood including all constants
  if (!prior_only) {
    target += normal_lpdf(logloss | mu, sig[dev_lag]);
  }
}

generated quantities{
  vector[n_data] log_lik;

  for (i in 1:n_data){
    log_lik[i] = normal_lpdf(logloss[i] | mu[i], sig[dev_lag[i]]);
  }
}

