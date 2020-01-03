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
}

parameters {
  real r_alpha[n_origin - 1];
  real r_beta[n_dev - 1];
  real<lower=-4.0,upper=4.0> logelr;
  real<lower=0,upper=100000> a_ig[n_dev];
  real gamma;
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
  beta[n_dev] = 0;
  
  speedup[1] = 1;
  for (i in 2:n_origin) speedup[i] = speedup[i-1] * (1 - gamma);
  
  sig[n_dev] = gamma_cdf(1/a_ig[n_dev],1,1);
  for (i in 1:(n_dev - 1)) sig[n_dev-i] = sig[n_dev + 1 - i] + gamma_cdf(1 / a_ig[i], 1, 1);
  sig = sqrt(sig);
  
  for (i in 1:n_data){
    mu[i] = logprem[i] + logelr + alpha[origin_lag[i]] + beta[dev_lag[i]] * speedup[origin_lag[i]];
  }
}

model {
  target += student_t_lpdf(r_alpha | 3, 0, 1);
  target += student_t_lpdf(r_beta | 3, 0, 1);
  target += inv_gamma_lpdf(a_ig | 1, 1);
  target += normal_lpdf(logelr | -0.4, 3.162);
  target += student_t_lpdf(gamma | 4, 0, 0.1);
  
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

