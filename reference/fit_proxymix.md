# Fit a Gaussian-mixture proxy to a target density

The unified front door of `proxymix`. Picks a fitting regime (or honours
an explicit choice) and dispatches to the corresponding regime-specific
fitter:

## Usage

``` r
fit_proxymix(
  target,
  N = 1L,
  regime = c("auto", "moment", "sample", "kld"),
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

  One of `"auto"`, `"moment"`, `"sample"`, `"kld"`.

- ...:

  Additional arguments forwarded to the regime-specific fitter. The most
  useful pass-throughs are `canonicalise` (whether to apply
  [`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md)
  to the returned fit; default `TRUE`), `validation_size` and
  `validation_proposal` (held-out IS validation for regime `"kld"`),
  `max_iter`, `tol`, `n_starts`, and `seed`.

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

With `regime = "auto"` the choice is made from the shape of the supplied
target:

- `N == 1` and the target carries samples or moments: `"moment"`.

- `N >= 2` and the target carries samples: `"sample"`.

- The target carries `log_density` only (no samples): `"kld"`.

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
#>   [1] w = 0.7928, |mu| = 0.3344, tr(Sigma) = 1.4380
#>   [2] w = 0.2072, |mu| = 1.5704, tr(Sigma) = 1.3301

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
```
