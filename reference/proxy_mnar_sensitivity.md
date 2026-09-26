# Missing-not-at-random sensitivity analysis for a coordinate mean

Sweeps the missing-not-at-random sensitivity slope `beta` over a grid
and, at each value, multiply-imputes `coord` under the selection model
\\P(\text{missing}\mid y) = g(\alpha + \beta y)\\ and pools its mean by
Rubin's rules. The result traces how the estimate and its confidence
interval move as the assumed dependence of missingness on the unobserved
value strengthens, so an analyst can read off the value of `beta` at
which a conclusion would change. `beta = 0` is the missing-at-random
anchor.

## Usage

``` r
proxy_mnar_sensitivity(
  data,
  coord,
  beta_grid = seq(0, 1, by = 0.25),
  link = c("logit", "probit"),
  N = NULL,
  m = 20L,
  seed = NULL,
  max_iter = 500L,
  ...
)
```

## Arguments

- data:

  A numeric matrix or data frame with `NA` in `coord` only (its other
  columns must be fully observed).

- coord:

  Name or index of the coordinate the mechanism acts on.

- beta_grid:

  Numeric vector of sensitivity slopes. Positive values make larger
  unobserved values more likely to be missing.

- link:

  Selection link, `"logit"` (the default) or `"probit"`.

- N, m, seed, ...:

  Passed to
  [`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md);
  a single `seed` makes the whole sweep reproducible and keeps the curve
  smooth across the grid.

- max_iter:

  Maximum EM iterations per fit, passed to
  [`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md).
  Default `500L`; fits near missing at random can need more than the
  [`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md)
  default of 100.

## Value

A data frame with one row per grid value: `beta`, `estimate`,
`std.error`, `conf.low`, `conf.high`, `fmi`, `loglik` (the observed-data
log-likelihood of the selection model at the fitted mixture and
intercept) and `converged` (whether the point fit converged). A warning
names any slope whose fit did not converge.

## Details

The mixture is refitted at each slope. Under this selection model the
observed data carry information about the slope, but only through the
assumed shape of the outcome distribution; a different shape can fit the
observed data equally well at a different slope. The `loglik` column
shows how strongly the data favour one slope over another under the
mixture. Treat a large difference as a consequence of the model's shape,
not as evidence about the missingness mechanism itself, and report the
whole curve.

## See also

[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`mnar()`](https://max578.github.io/proxymix/reference/mechanism.md).

Other imputation:
[`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md),
[`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md),
[`gmm_imputation()`](https://max578.github.io/proxymix/reference/gmm_imputation.md),
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`mechanism`](https://max578.github.io/proxymix/reference/mechanism.md),
[`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md),
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)

## Examples

``` r
set.seed(1)
x1 <- rnorm(300)
y <- x1 + rnorm(300)
y[runif(300) < plogis(-0.4 + 0.8 * y)] <- NA      # MNAR on y
dat <- data.frame(x1 = x1, y = y)
proxy_mnar_sensitivity(dat, "y", beta_grid = c(0, 0.5, 1), m = 5L, seed = 1L)
#>   beta    estimate  std.error    conf.low    conf.high       fmi    loglik
#> 1  0.0 -0.19648308 0.09406029 -0.38446678 -0.008499378 0.2526332 -857.1654
#> 2  0.5  0.01881364 0.10018178 -0.18144275  0.219070025 0.2538684 -836.8685
#> 3  1.0  0.21915436 0.11313359 -0.01039855  0.448707273 0.3355709 -832.4386
#>   converged
#> 1      TRUE
#> 2      TRUE
#> 3      TRUE
```
