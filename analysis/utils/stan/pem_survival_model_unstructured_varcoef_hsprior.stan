/*  Variable naming:
 // dimensions
 N          = total number of observations (length of data)
 S          = number of sample ids
 T          = max timepoint (number of timepoint ids)
 M          = number of covariates
 G          = number of groups
 
 // main data matrix (per observed timepoint*record)
 s          = sample id for each obs
 t          = timepoint id for each obs
 event      = integer indicating if there was an event at time t for sample s
 x          = matrix of real-valued covariates at time t for sample n [N, X]
 g          = group id for each obs
 
 // timepoint-specific data (per timepoint, ordered by timepoint id)
 t_obs      = observed time since origin for each timepoint id (end of period)
 t_dur      = duration of each timepoint period (first diff of t_obs)
 
*/
// Jacqueline Buros Novik <jackinovik@gmail.com>

functions {
  // function borrowed from https://github.com/to-mi/stan-survival-shrinkage/blob/master/wei_hs.stan
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;

    for (m in 1:dims(x)[1]){
      res[m] <- sqrt(x[m]);
    }

    return res;
  }

  // function borrowed from https://github.com/to-mi/stan-survival-shrinkage/blob/master/wei_hs.stan
  vector hs_prior_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, real nu) {
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);

    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu, 0.5 * nu);

    return (r1_global * sqrt(r2_global)) * r1_local .* sqrt_vec(r2_local);
  }
}
data {
  // dimensions
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> M;
  int<lower=1> G;
  
  // data matrix
  int<lower=1, upper=N> s[N];     // sample id
  int<lower=1, upper=T> t[N];     // timepoint id
  int<lower=0, upper=1> event[N]; // 1: event, 0:censor
  matrix[N, M] x;                 // explanatory vars
  int<lower=1, upper=G> g[N];     // group ids
  
  // timepoint data
  vector<lower=0>[T] t_obs;
  vector<lower=0>[T] t_dur;
  
  // param to regularizing prior
  real<lower=1> nu;
}
transformed data {
  vector[T] log_t_dur;  // log-duration for each timepoint
  
  log_t_dur = log(t_obs);
}
parameters {
  vector[M] beta;             // overall beta for each covariate
  vector[M] grp_beta_raw[G];  // group-level beta for each covariate  
  real<lower=0> grp_beta_sigma[G];
  
  vector[T] log_baseline_raw; // unstructured baseline hazard for each timepoint t
  real<lower=0> baseline_sigma;
  real log_baseline_mu;
  
  // group-level alphas
  real grp_alpha[G];
  real<lower=0> grp_alpha_sigma;
  
  // parameters to regularization of betas
  real<lower=0> tau_s1_raw[G];
  real<lower=0> tau_s2_raw[G];
  vector<lower=0>[M] tau1_raw[G];
  vector<lower=0>[M] tau2_raw[G];
}
transformed parameters {
  vector[N] log_hazard;
  vector[T] log_baseline;       // unstructured baseline hazard for each timepoint t
  vector[M] grp_beta_unreg[G];
  vector[M] grp_beta[G];
  
  for (grp in 1:G) {
      grp_beta_unreg[grp] = beta + grp_beta_raw[grp];
      grp_beta[grp] <- hs_prior_lp(tau_s1_raw[grp], tau_s2_raw[grp], tau1_raw[grp], tau2_raw[grp], nu) .* grp_beta_unreg[grp];
  }
  
  log_baseline = log_baseline_raw + log_t_dur;
  
  for (n in 1:N) {
    log_hazard[n] = log_baseline_mu + grp_alpha[g[n]] + log_baseline[t[n]] + x[n,]*grp_beta[g[n]];
  }
}
model {
  beta ~ cauchy(0, 2);
  event ~ poisson_log(log_hazard);
  log_baseline_mu ~ normal(0, 1);
  grp_alpha_sigma ~ cauchy(0, 1);
  grp_beta_sigma ~ cauchy(0, 1);
  grp_alpha ~ normal(0, grp_alpha_sigma);
  for (grp in 1:G)
      grp_beta_raw[grp] ~ cauchy(0, grp_beta_sigma[grp]);
  baseline_sigma ~ normal(0, 1);
  log_baseline_raw ~ normal(0, baseline_sigma);
}
generated quantities {
  real log_lik[N];
  //int yhat_uncens[N];
  vector[T] baseline_raw;
  // matrix<lower=0>[T, G] grp_baseline;
  
  // compute raw baseline hazard, for summary/plotting
  baseline_raw = exp(log_baseline_raw);
  
  // prepare yhat_uncens & log_lik
  for (n in 1:N) {
      //yhat_uncens[n] = poisson_log_rng(log_hazard[n]);
      log_lik[n] = poisson_log_log(event[n], log_hazard[n]);
  }
}