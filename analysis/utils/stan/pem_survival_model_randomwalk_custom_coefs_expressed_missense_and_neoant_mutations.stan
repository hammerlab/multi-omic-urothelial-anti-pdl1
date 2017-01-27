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
  int<lower=0> observed_neoant_count[N];
  int<lower=0> observed_snv_count[N];
  int<lower=0> expressed_missense_count[N];
  int<lower=0> expressed_neoant_count[N];
  matrix[N, M] x;                 // explanatory vars
}
transformed data {
  real<lower=0> t_dur[T];  // duration for each timepoint
  real<lower=0> t_obs[T];  // observed end time for each timepoint
  int<lower=0> measurement_inputs[S, 5]; // matrix of observed measurements
  vector[1] one;
  
  one[1] <- 1;
  
  // capture per-sample observed counts 
  // since these don't vary over time
  for (n in 1:N) {
      measurement_inputs[s[n], 1] = observed_snv_count[n];
      measurement_inputs[s[n], 2] = observed_missense_count[n];
      measurement_inputs[s[n], 3] = observed_neoant_count[n];
      measurement_inputs[s[n], 4] = expressed_missense_count[n];
      measurement_inputs[s[n], 5] = expressed_neoant_count[n];
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
  
  // fixed covariate betas for outcome model
  vector[M] beta;
  vector[1] beta_antigen_potential;
  
  // terms for submodels estimating latent effect of antigenic potential 
  vector[4] raw_antigen_submodel_alpha;
  vector[4] raw_antigen_submodel_beta;
  vector<lower=0>[5] antigen_submodel_phi;
  
  // latent parameters estimated from data
  real antigenic_potential[S]; // subject-level parameter for antigenic load
}
transformed parameters {
  vector<lower=0>[N] hazard;
  vector<lower=0>[T] baseline;
  vector[5] antigen_submodel_alpha;
  vector[5] antigen_submodel_beta;
  vector<lower=0>[5] lp_antigen_load[S];
  
  // constrain first loading factors to be one, to ensure identifiability
  antigen_submodel_alpha = append_row(one, raw_antigen_submodel_alpha);
  antigen_submodel_beta = append_row(one, raw_antigen_submodel_beta);
  
  // prep linear predictor for antigenic load measurements 
  for (sample in 1:S) {
      lp_antigen_load[sample] = exp(antigen_submodel_alpha + antigen_submodel_beta * antigenic_potential[sample]);
  }
  
  // transform baseline hazard depending on duration of interval
  for (i in 1:T) {
    baseline[i] <- baseline_raw[i]*t_dur[i];
  }

  // linear predictor for hazard
  for (n in 1:N) {
    hazard[n] <- exp(x[n,]*beta + 
                    antigenic_potential[s[n]]*beta_antigen_potential[1]
                    )*baseline[t[n]];
  }
}
model {
  // priors on latent submodel parameters
  raw_antigen_submodel_alpha ~ normal(0, 1);
  raw_antigen_submodel_beta ~ normal(0, 1);
  antigen_submodel_phi ~ normal(0, 1);
  
  // priors on linear predictors of hazard
  beta ~ normal(0, 1);
  beta_antigen_potential ~ normal(0, 1);
  
  // priors on baseline hazards
  baseline_loc ~ normal(0, 1);
  baseline_raw[1] ~ normal(baseline_loc, 1);
  for (i in 2:T) {
      baseline_raw[i] ~ normal(baseline_raw[i-1], baseline_sigma);
  }
  
  // measurement model for antigenic_potential
  for (sample in 1:S) {
      measurement_inputs[sample,] ~ neg_binomial_2(lp_antigen_load[sample], antigen_submodel_phi);
  }
  
  // model for outcome
  event ~ poisson(hazard);
}
generated quantities {
  real log_lik[N];

  for (n in 1:N) {
      log_lik[n] <- poisson_log(event[n], hazard[n]);
  }
}