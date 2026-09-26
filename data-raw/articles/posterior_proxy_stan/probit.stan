// Probit regression with independent normal(0, prior_sd) priors on every
// coefficient, the intercept included; the log likelihood goes through the
// standard normal log CDF so that extreme linear predictors stay finite.
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
  vector[n] eta = X * beta;
  beta ~ normal(0, prior_sd);
  for (i in 1:n) {
    target += std_normal_lcdf(y[i] == 1 ? eta[i] : -eta[i]);
  }
}
