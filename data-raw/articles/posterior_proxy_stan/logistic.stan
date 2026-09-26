// Logistic regression with independent normal(0, prior_sd) priors on every
// coefficient, the intercept included.
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] X;
  array[n] int<lower=0, upper=1> y;
  real<lower=0> prior_sd;
}
parameters {
  vector[p] beta;
}
model {
  beta ~ normal(0, prior_sd);
  y ~ bernoulli_logit(X * beta);
}
