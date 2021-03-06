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

functions {
  matrix spline(vector x, int N, int H, vector xi, int P) {
    matrix[N, H + P] b_x;         // expanded predictors
    for (n in 1:N) {
        for (p in 1:P) {
            b_x[n,p] <- pow(x[n],p-1);  // x[n]^(p-1)
        }
        for (h in 1:H)
          b_x[n, h + P] <- fmax(0, pow(x[n] - xi[h],P-1)); 
    }
    return b_x;
  }
}
data {
  // dimensions
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> M;
  
  // data matrix
  int<lower=1, upper=N> s[N];     // sample id
  int<lower=1, upper=T> t[N];     // timepoint id
  int<lower=0, upper=1> event[N]; // 1: event, 0:censor
  matrix[N, M] x;                 // explanatory vars

  // timepoint data
  vector<lower=0>[T] t_obs;
  vector<lower=0>[T] t_dur;
      
  // spline params (for tvc)
  int<lower=1> H;                 // number of knots (fixed)
  int<lower=0> power;             // power of spline (1:linear, 2:quad, 3:cubic)
}
transformed data {
  int<lower=1> P;
  vector[H+1] xi_prior;
  vector[T] log_t_dur;
  
  log_t_dur = log(t_obs);
  
  for (h in 1:(H+1)) {
      xi_prior[h] <- 1;
  }
  
  P <- 1+power;
}
parameters {
  vector[T] log_baseline_raw;    // unstructured baseline hazard for each timepoint t
  real<lower=0> baseline_sigma;
  real log_baseline_mu;
  
  vector[H+P] beta_time_spline[M];    // time-spline coefficients for each beta
  simplex[H+1] xi_proportions; // time-segment proportions
}
transformed parameters {
  vector[N] log_hazard;
  vector[T] log_baseline;
  matrix[T, H+P] time_spline;
  vector[T] beta_time[M];
  positive_ordered[H] est_xi;          // locations of knots
  
  // adjust baseline hazard for duration of each period
  log_baseline = log_baseline_raw + log_t_dur;
  
  // xi as proportions of max observed time
  est_xi <- cumulative_sum(xi_proportions[1:H])*max(t_obs);
  
  // coefficients for each timepoint T
  time_spline <- spline(t_obs, T, H, est_xi, P);
  for (m in 1:M) {
      beta_time[m] <- time_spline*beta_time_spline[m];
  }

  for (n in 1:N) {
    real linpred;
    linpred <- 0;
    for (m in 1:M) {
      // for now, handle each M separately
      // (to be sure we pull out the "right" beta.. )
      linpred <- linpred + x[n, m] * beta_time[m][t[n]]; 
    }
    // log hazard for each obs N
    log_hazard[n] = log_baseline_mu + log_baseline[t[n]] + linpred;
  }
}
model {
  xi_proportions ~ dirichlet(xi_prior);
  for (m in 1:M) {
    beta_time_spline[m] ~ normal(0, 1);
  }
  event ~ poisson_log(log_hazard);
  log_baseline_mu ~ normal(0, 1);
  baseline_sigma ~ normal(0, 1);
  log_baseline_raw[1] ~ normal(0, 1);
  for (i in 2:T) {
      log_baseline_raw[i] ~ normal(log_baseline_raw[i-1], baseline_sigma);
  }
}
generated quantities {
  real log_lik[N];
  real beta[M];
  vector[T] baseline_raw;
  
  // compute raw baseline hazard, for summary/plotting
  baseline_raw = exp(log_baseline_raw);
  
  // log_lik for use with loo
  for (n in 1:N) {
      //yhat_uncens[n] = poisson_log_rng(log_hazard[n]);
      log_lik[n] <- poisson_log_log(event[n], log_hazard[n]);
  }
  
  // pull out linear component of spline as primary beta coef
  for (m in 1:M) {
    beta[m] <- beta_time_spline[m][1];
  }
  
}