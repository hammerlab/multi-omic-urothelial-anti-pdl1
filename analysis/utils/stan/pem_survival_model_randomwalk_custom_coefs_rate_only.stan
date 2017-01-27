/*  Variable naming:
 // dimensions
 N          = total number of observations (length of data)
 S          = number of sample ids
 T          = max timepoint (number of timepoint ids)
 M          = number of covariates
 
 // data
 s          = sample id for each obs
 t          = timepoint id for each obs
 event      = integer indicating if there was an event at time t for sample s
 x          = matrix of real-valued covariates at time t for sample n [N, X]
 obs_t      = observed end time for interval for timepoint for that obs
 
*/
// Jacqueline Buros Novik <jackinovik@gmail.com>

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1> M;
  int<lower=1, upper=N> s[N];     // sample id
  int<lower=1, upper=T> t[N];     // timepoint id
  int<lower=0, upper=1> event[N]; // 1: event, 0:censor
  real<lower=0> obs_t[N];         // observed end time for each obs
  int<lower=0> observed_missense_count[N];
  int<lower=0> observed_snv_count[N];
  matrix[N, M] x;                 // explanatory vars
}
transformed data {
  real<lower=0> t_dur[T];  // duration for each timepoint
  real<lower=0> t_obs[T];  // observed end time for each timepoint
  int<lower=0> s_observed_snv_count[S];
  int<lower=0> s_observed_missense_count[S];
  
  // capture per-sample observed counts 
  // since these don't vary over time
  for (n in 1:N) {
      s_observed_snv_count[s[n]] = observed_snv_count[n];
      s_observed_missense_count[s[n]] = observed_missense_count[n];
  }

  // capture observation time for each timepoint id t
  for (i in 1:N) {
      // assume these are constant per id across samples
      t_obs[t[i]] <- obs_t[i];  
  }
  
  // duration of each timepoint
  // duration at first timepoint = t_obs[1] ( implicit t0 = 0 )
  t_dur[1] <- t_obs[1];
  for (i in 2:T) {
      t_dur[i] <- t_obs[i] - t_obs[i-1];
  }
}
parameters {
  // baseline hazard 
  vector<lower=0>[T] baseline_raw; // unstructured baseline hazard for each timepoint t
  real<lower=0> baseline_sigma;
  real<lower=0> baseline_loc;
  vector[M] beta;
  vector[1] beta_missense;
  
  // estimate missense_snv_rate from data
  real<lower=0, upper=1> missense_snv_rate[S];
}
transformed parameters {
  vector<lower=0>[N] hazard;
  vector<lower=0>[T] baseline;
  
  for (i in 1:T) {
    baseline[i] <- baseline_raw[i]*t_dur[i];
  }
  
  for (n in 1:N) {
    hazard[n] <- exp(x[n,]*beta + 
                    missense_snv_rate[s[n]]*beta_missense[1]
                    )*baseline[t[n]];
  }
}
model {
  // priors on missense rate 
  missense_snv_rate ~ normal(0, 1);
  beta_missense ~ normal(0, 1);
  beta ~ normal(0, 1);
  
  // priors on baseline hazards
  baseline_loc ~ normal(0, 1);
  baseline_raw[1] ~ normal(baseline_loc, 1);
  for (i in 2:T) {
      baseline_raw[i] ~ normal(baseline_raw[i-1], baseline_sigma);
  }
  
  // models
  s_observed_missense_count ~ binomial(s_observed_snv_count, missense_snv_rate);
  event ~ poisson(hazard);
}
generated quantities {
  real log_lik[N];

  for (n in 1:N) {
      log_lik[n] <- poisson_log(event[n], hazard[n]);
  }
}