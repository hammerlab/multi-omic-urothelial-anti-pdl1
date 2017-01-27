data {
    int N; // number of obs (pregnancies)
    int M; // number of predictors
  
    int<lower=0, upper=1> event[N]; // outcome
    matrix[N,M] x; // predictors
}
parameters {
    real alpha;
    vector[M] beta;
}
model {
    alpha ~ normal(0, 1);
    beta ~ normal(0, 1);
    event ~ bernoulli_logit(alpha + x*beta);
}
generated quantities {
    real log_lik[N];
    int y_rep[N];
    
    for (n in 1:N) {
        log_lik[n] = bernoulli_logit_lpmf(event[n] | alpha + x[n]*beta);
        y_rep[n] = bernoulli_rng(inv_logit(alpha + x[n]*beta));
    }
}