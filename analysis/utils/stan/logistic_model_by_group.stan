data {
    int N; // number of obs (pregnancies)
    int M; // number of predictors
    int G; // number of strata
    
    int<lower=0, upper=1> event[N]; // outcome
    matrix[N,M] x; // predictors
    int<lower=0, upper=G> g[N]; // group id for each obs
}
parameters {
    real alpha;
    vector[M] beta;
    vector[G] grp_alpha;
    matrix[M, G] grp_beta;
    real<lower=0> grp_alpha_sigma;
    real<lower=0> grp_beta_sigma;
}
transformed parameters {
    vector[N] linear_pred;
    
    for (n in 1:N) {
        linear_pred[n] = grp_alpha[g[n]] + x[n]*grp_beta[,g[n]];
    }
}
model {
    alpha ~ normal(0, 1);
    beta ~ normal(0, 1);
    grp_alpha ~ normal(alpha, grp_alpha_sigma);
    for (grp in 1:G) {
        grp_beta[,grp] ~ normal(beta, grp_beta_sigma);
    }
    grp_alpha_sigma ~ cauchy(0, 2.5);
    grp_beta_sigma ~ cauchy(0, 2.5);
    event ~ bernoulli_logit(linear_pred);
}
generated quantities {
    real log_lik[N];
    int y_rep[N];
    
    for (n in 1:N) {
        log_lik[n] = bernoulli_logit_lpmf(event[n] | linear_pred[n]);
        y_rep[n] = bernoulli_rng(inv_logit(linear_pred[n]));
    }
}