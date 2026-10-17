# Mixture of probabilistic PCA fitted to samples

Fits a Gaussian mixture whose component covariances are restricted to
\\\Sigma_k = W_k W_k^\top + \sigma_k^2 I\\, where \\W_k\\ has `q`
columns (Tipping and Bishop, 1999). Each component varies mainly in a
`q`-dimensional subspace, with equal noise variance in every other
direction. With `p` variables a component has \\pq - q(q-1)/2 + 1\\
covariance parameters instead of \\p(p+1)/2\\, so the model can be
fitted when a component holds fewer points than there are variables,
where a full-covariance fit
([`fit_em_samples()`](https://max578.github.io/proxymix/dev/reference/fit_em_samples.md))
is singular.

## Usage

``` r
fit_mppca_samples(
  target,
  N = 2L,
  q = 1L,
  n_starts = 5L,
  max_iter = 500L,
  tol = 1e-08,
  seed = NULL,
  canonicalise = TRUE
)
```

## Arguments

- target:

  A
  [gmm_target](https://max578.github.io/proxymix/dev/reference/gmm_target.md)
  carrying an `n` by `p` `samples` matrix.

- N:

  Number of mixture components.

- q:

  Dimension of each component's principal subspace, an integer with
  `1 <= q < p`. Compare values of `q` (and `N`) with the BIC in the
  diagnostics.

- n_starts:

  Number of k-means starts. The fit with the largest final
  log-likelihood is returned.

- max_iter:

  Maximum number of EM iterations per start.

- tol:

  Relative log-likelihood convergence tolerance.

- seed:

  Optional integer seed for the k-means starts.

- canonicalise:

  Logical. If `TRUE` (the default), components are sorted by
  [`gmm_canonicalise()`](https://max578.github.io/proxymix/dev/reference/gmm_canonicalise.md).

## Value

A [gmm_fit](https://max578.github.io/proxymix/dev/reference/gmm_fit.md)
with `regime = "sample"`. The diagnostics carry the log-likelihood trace
and final value, the number of free parameters, BIC and AIC.
`metadata$mppca` holds `q`, the loading matrices `loadings` and the
noise variances `noise_var`, in component order.

## Details

The fit is by EM. Given the responsibilities, the update of \\(W_k,
\sigma_k^2)\\ is the closed-form maximum-likelihood solution of
probabilistic PCA on the responsibility-weighted covariance \\S_k\\: the
columns of \\W_k\\ are the leading `q` eigenvectors of \\S_k\\, scaled
by \\(\lambda_j - \sigma_k^2)^{1/2}\\, and \\\sigma_k^2\\ is the mean of
the remaining \\p - q\\ eigenvalues. Starts are k-means partitions.

The result is an ordinary
[gmm_fit](https://max578.github.io/proxymix/dev/reference/gmm_fit.md),
so every operator of the package (conditioning, marginals, imputation,
sampling) applies to it.

## References

Tipping, M. E. and Bishop, C. M. (1999). Mixtures of probabilistic
principal component analyzers. *Neural Computation*, 11(2), 443-482.
[doi:10.1162/089976699300016728](https://doi.org/10.1162/089976699300016728)

## See also

Other fitting:
[`fit_em_samples()`](https://max578.github.io/proxymix/dev/reference/fit_em_samples.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/dev/reference/fit_kld_em.md),
[`fit_moment_match()`](https://max578.github.io/proxymix/dev/reference/fit_moment_match.md),
[`from_cf()`](https://max578.github.io/proxymix/dev/reference/from_cf.md),
[`from_cf_elliptical()`](https://max578.github.io/proxymix/dev/reference/from_cf_elliptical.md),
[`from_cf_sources()`](https://max578.github.io/proxymix/dev/reference/from_cf_sources.md),
[`from_kde()`](https://max578.github.io/proxymix/dev/reference/from_kde.md),
[`from_objective()`](https://max578.github.io/proxymix/dev/reference/from_objective.md),
[`select_N()`](https://max578.github.io/proxymix/dev/reference/select_N.md)

## Examples

``` r
set.seed(1)
p <- 10
w1 <- qr.Q(qr(matrix(rnorm(p * 2), p)))
w2 <- qr.Q(qr(matrix(rnorm(p * 2), p)))
x <- rbind(matrix(rnorm(60 * 2), 60) %*% t(3 * w1) + matrix(rnorm(60 * p, sd = 0.3), 60),
           matrix(rnorm(60 * 2), 60) %*% t(3 * w2) + matrix(rnorm(60 * p, sd = 0.3), 60) + 4)
fit <- fit_mppca_samples(gmm_target_from_samples(x), N = 2L, q = 2L, seed = 1L)
fit@metadata$mppca$noise_var
#> [1] 0.08906958 0.09604513
fit@diagnostics$bic
#> [1] 2100.344
```
