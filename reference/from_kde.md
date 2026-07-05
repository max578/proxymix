# Compile a kernel-density estimate into a Gaussian-mixture proxy

Fits an `N`-component Gaussian-mixture proxy to a (Gaussian, diagonal-
bandwidth) kernel-density estimate over `samples`, via regime (iii)
KLD-EM. The proxy is closed-form marginalisable, conditionable, and
samplable; the KDE is none of those things on its own.

## Usage

``` r
from_kde(
  samples,
  N = 3L,
  bandwidth = "silverman",
  proposal = NULL,
  is_size = 5000L,
  max_iter = 100L,
  tol = 1e-05,
  ridge_eps = 1e-06,
  min_ess = 50L,
  seed = NULL,
  validation_size = 0L,
  validation_proposal = NULL,
  validation_seed = NULL,
  support_warn = TRUE,
  canonicalise = TRUE
)
```

## Arguments

- samples:

  An `n` by `p` numeric matrix of points. `n >= 5`, `p <= 10`.

- N:

  Number of mixture components in the proxy.

- bandwidth:

  Either `"silverman"`, `"scott"`, a positive numeric scalar (absolute
  bandwidth applied to every coordinate), or a length-`p` positive
  numeric vector of per-coordinate absolute bandwidths. Default
  `"silverman"`.

- proposal:

  Optional
  [is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md).
  Default is a multivariate-t centred at `colMeans(samples)`, scale =
  ridge(cov(samples)) + diag(h^2), `df = 5`.

- is_size:

  Importance-sample size for fitting. Default `5000L`.

- max_iter:

  Maximum EM iterations. Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- tol:

  Convergence tolerance. Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- ridge_eps:

  Ridge added to each component covariance at every M-step. Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- min_ess:

  Minimum effective sample size below which a warning is issued.
  Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- seed:

  Optional integer seed for the fitting IS draw.

- validation_size:

  Held-out IS sample size. Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- validation_proposal:

  Optional
  [is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md)
  for the held-out sample. Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- validation_seed:

  Optional integer seed for the held-out IS draw. Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- support_warn:

  Logical. Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

- canonicalise:

  Logical. If `TRUE`, the fitted mixture is post-processed by
  [`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md).
  Forwarded to
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md).

## Value

A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) with
`regime = "kld"` and metadata recording the KDE inputs (`kde_samples_n`,
`bandwidth`, `bandwidth_method`).

## Details

This is a **compression** operation: take an `n`-sample KDE and replace
it with the closest `N`-component mixture in the Kullback-Leibler sense
(which is much smaller than `n` for typical use). Bias inherited from
the KDE is reproduced in the proxy; the bandwidth controls the
bias-variance trade-off.

Dimensional scope. The dimensional guard is `p <= 5` (recommended),
`p <= 10` (allowed with warning), `p > 10` (rejected). Regime-(iii)
KLD-EM is driven by importance sampling, whose effective sample size
collapses sharply in high dimensions.

## See also

Other fitting:
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md),
[`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md),
[`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md),
[`select_N()`](https://max578.github.io/proxymix/reference/select_N.md)

## Examples

``` r
set.seed(1L)
x <- rbind(
  mvnfast::rmvn(120L, mu = c(-2, 0), sigma = diag(2)),
  mvnfast::rmvn(120L, mu = c( 2, 0), sigma = diag(2))
)
fit <- from_kde(x, N = 2L, is_size = 2000L, max_iter = 40L, seed = 1L)
fit
#> <gmm_fit>: regime = "kld", K = 2, p = 2
#>   target     : from_kde
#>   iterations : 16
#>   converged  : TRUE
#>   [1] w = 0.5014, |mu| = 2.0118, tr(Sigma) = 3.1904
#>   [2] w = 0.4986, |mu| = 2.1807, tr(Sigma) = 2.4859
ess_summary(fit)
#> $is_size
#> [1] 2000
#> 
#> $ess
#> [1] 1552.231
#> 
#> $ess_relative
#> [1] 0.7761155
#> 
#> $max_weight
#> [1] 0.001135827
#> 
#> $support_fraction
#> [1] 1
#> 
#> $mc_se_kld
#> [1] 0.005537856
#> 
#> $validation_size
#> [1] 0
#> 
#> $validation_ess
#> [1] NA
#> 
#> $validation_ess_relative
#> [1] NA
#> 
#> $validation_max_weight
#> [1] NA
#> 
#> $validation_support_fraction
#> [1] NA
#> 
#> $validation_kld
#> [1] NA
#> 
```
