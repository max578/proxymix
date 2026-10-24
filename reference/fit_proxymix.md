# Fit a Gaussian-mixture proxy to a target distribution

The unified front door of `proxymix`. Picks a fitting regime (or honours
an explicit choice) and dispatches to the corresponding regime-specific
fitter:

## Usage

``` r
fit_proxymix(
  target,
  N = 1L,
  regime = c("auto", "moment", "sample", "kld", "cf"),
  ...
)
```

## Arguments

- target:

  A
  [gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md).

- N:

  Number of components.

- regime:

  One of `"auto"`, `"moment"`, `"sample"`, `"kld"`, `"cf"`.

- ...:

  Additional arguments forwarded to the regime-specific fitter. The most
  useful pass-throughs are `canonicalise` (whether to apply
  [`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md)
  to the returned fit; default `TRUE`), `validation_size` and
  `validation_proposal` (held-out IS validation for regime `"kld"`),
  `max_iter`, `tol`, `n_starts`, and `seed`. For regime `"cf"` the
  arguments of
  [`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md)
  apply, among them `n_freq`, `n_start`, `seed` and `bracket`; its fit
  is always returned in canonical order, and `canonicalise` is ignored.

## Value

A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md).

## Details

- `"moment"` - closed-form moment matching
  ([`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md)).

- `"sample"` - classical EM on i.i.d. samples
  ([`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md)).

- `"kld"` - importance-sampled KLD-EM
  ([`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md))
  for a target that can be evaluated but not sampled.

- `"cf"` - a fit to the target's characteristic function
  ([`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md)),
  for a target known only through it
  ([`cf_target()`](https://max578.github.io/proxymix/reference/cf_target.md)).

With `regime = "auto"` the choice is made from the shape of the supplied
target:

- `N == 1` and the target carries samples or moments: `"moment"`.

- `N >= 2` and the target carries samples: `"sample"`.

- The target carries `log_density` only (no samples): `"kld"`.

- The target carries `cf` only (no samples, no `log_density`): `"cf"`.

## Examples

``` r
## auto: samples + N=2 -> classical EM.
x <- matrix(stats::rnorm(200), ncol = 2)
tgt_s <- gmm_target_from_samples(x)
fit_proxymix(tgt_s, N = 2L, max_iter = 25L)
#> Auto-selected regime: "sample".
#> <gmm_fit>: regime = "sample", K = 2, p = 2
#>   target     : target_from_samples
#>   iterations : 25
#>   converged  : FALSE
#>   [1] w = 0.6478, |mu| = 0.1673, tr(Sigma) = 1.6547
#>   [2] w = 0.3522, |mu| = 0.1282, tr(Sigma) = 2.4231

## explicit "kld" on a log-density-only target.
fit_proxymix(banana_target(), N = 3L, regime = "kld",
             is_size = 1000L, max_iter = 20L, seed = 1L)
#> <gmm_fit>: regime = "kld", K = 3, p = 2
#>   target     : banana
#>   iterations : 20
#>   converged  : FALSE
#>   [1] w = 0.6103, |mu| = 0.3840, tr(Sigma) = 1.3456
#>   [2] w = 0.3228, |mu| = 0.8867, tr(Sigma) = 2.4049
#>   [3] w = 0.0670, |mu| = 2.4700, tr(Sigma) = 5.6889

## auto on a target known only through its characteristic function:
## a Gamma(3) variable plus N(0, 0.25) noise -> "cf".
tgt_cf <- cf_target(
  function(t) (1 - 1i * t[, 1])^(-3) * exp(-0.125 * t[, 1]^2),
  n_dim = 1L
)
fit_proxymix(tgt_cf, N = 3L, n_freq = 1024L, seed = 1L)
#> Auto-selected regime: "cf".
#> <gmm_fit>: regime = "cf", K = 3, p = 1
#>   target     : cf_target
#>   iterations : NA
#>   converged  : TRUE
#>   cf held-out RMS error : 0.00477
#>   Kolmogorov distance   : at least 0.00398, at most 0.00859 (1 direction)
#>   [1] w = 0.4671, |mu| = 2.9662, tr(Sigma) = 1.4310
#>   [2] w = 0.3043, |mu| = 1.5676, tr(Sigma) = 0.6571
#>   [3] w = 0.2286, |mu| = 4.8617, tr(Sigma) = 3.6102
```
