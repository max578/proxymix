// Kernel density estimate with a common Gaussian kernel, as a target for NUTS.
data {
  int<lower=1> K;                 // number of kernel centres
  int<lower=1> D;                 // dimension
  array[K] vector[D] centre;      // kernel centres, the training data
  cov_matrix[D] H;                // kernel bandwidth matrix
}
transformed data {
  matrix[D, D] L = cholesky_decompose(H);
}
parameters {
  vector[D] x;
}
model {
  vector[K] lp;
  for (k in 1:K) lp[k] = multi_normal_cholesky_lpdf(x | centre[k], L);
  target += log_sum_exp(lp) - log(K);
}
