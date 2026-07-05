# Closed-form moment-matching fit

Implements regime (i) of Hoek and Elliott (2024). When `N == 1`, this is
the exact moment match: `mu` is the target mean and `Sigma` is the
target covariance. When `N > 1`, the function returns the deterministic
moment-seed of
[`init_moment_seed()`](https://max578.github.io/proxymix/reference/init_moment_seed.md)
wrapped as a
[gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md),
without iterative refinement — useful as a starting point for the
iterative regimes.

## Usage

``` r
fit_moment_match(target, N = 1L, ridge_eps = 1e-06, canonicalise = TRUE)
```

## Arguments

- target:

  A
  [gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md).

- N:

  Number of components. `N >= 2` returns a moment-seeded mixture without
  iterative refinement.

- ridge_eps:

  Ridge added to the empirical covariance for numerical stability.

- canonicalise:

  Logical. If `TRUE` (the default), the fitted mixture is post-processed
  by
  [`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md)
  so that components are sorted by descending weight and (as a
  tiebreaker) by descending `||mu||`. Set `FALSE` to retain the raw
  component order.

## Value

A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) with
`regime = "moment"`.

## Details

Either the `target` must carry an `n` by `p` `samples` matrix, or its
`metadata` slot must contain pre-computed `moments` of the form
`list(mean = <p-vec>, cov = <p-by-p>)`.

## See also

Other fitting:
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md),
[`from_kde()`](https://max578.github.io/proxymix/reference/from_kde.md),
[`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md),
[`select_N()`](https://max578.github.io/proxymix/reference/select_N.md)

## Examples

``` r
x <- matrix(stats::rnorm(200), ncol = 2)
tgt <- gmm_target_from_samples(x)
fit_moment_match(tgt, N = 1L)
#> <gmm_fit>: regime = "moment", K = 1, p = 2
#>   target     : target_from_samples
#>   iterations : 0
#>   converged  : TRUE
#>   [1] w = 1.0000, |mu| = 0.0897, tr(Sigma) = 2.0274
```
