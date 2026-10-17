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
  [`gmm_impute()`](https://max578.github.io/proxymix/dev/reference/gmm_impute.md);
  a single `seed` makes the whole sweep reproducible and keeps the curve
  smooth across the grid.

- max_iter:

  Maximum EM iterations per fit, passed to
  [`gmm_impute()`](https://max578.github.io/proxymix/dev/reference/gmm_impute.md).
  Default `500L`; fits near missing at random can need more than the
  [`gmm_impute()`](https://max578.github.io/proxymix/dev/reference/gmm_impute.md)
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

[`gmm_impute()`](https://max578.github.io/proxymix/dev/reference/gmm_impute.md),
[`mnar()`](https://max578.github.io/proxymix/dev/reference/mechanism.md).

Other imputation:
[`as_mids()`](https://max578.github.io/proxymix/dev/reference/as_mids.md),
[`gmm_complete()`](https://max578.github.io/proxymix/dev/reference/gmm_complete.md),
[`gmm_imputation()`](https://max578.github.io/proxymix/dev/reference/gmm_imputation.md),
[`gmm_impute()`](https://max578.github.io/proxymix/dev/reference/gmm_impute.md),
[`mechanism`](https://max578.github.io/proxymix/dev/reference/mechanism.md),
[`proxy_fmi()`](https://max578.github.io/proxymix/dev/reference/proxy_fmi.md),
[`proxy_pool()`](https://max578.github.io/proxymix/dev/reference/proxy_pool.md)

## Examples

``` r
set.seed(1)
x1 <- rnorm(200)
y <- x1 + rnorm(200)
y[runif(200) < plogis(-0.4 + 0.8 * y)] <- NA      # MNAR on y
dat <- data.frame(x1 = x1, y = y)
proxy_mnar_sensitivity(dat, "y", beta_grid = c(0, 0.5, 1), m = 5L, seed = 1L)
#>   beta   estimate std.error   conf.low   conf.high       fmi    loglik
#> 1  0.0 -0.3588626 0.1471407 -0.6864382 -0.03128708 0.6304986 -557.7112
#> 2  0.5 -0.1233208 0.1261837 -0.3886618  0.14202015 0.4743556 -552.2003
#> 3  1.0  0.1370836 0.1242975 -0.1173614  0.39152868 0.3752030 -552.1423
#>   converged
#> 1      TRUE
#> 2      TRUE
#> 3      TRUE
```
