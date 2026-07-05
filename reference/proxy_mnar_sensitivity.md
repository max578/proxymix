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

## Value

A data frame with one row per grid value: `beta`, `estimate`,
`std.error`, `conf.low`, `conf.high`, `fmi`.

## Details

The slope is a sensitivity parameter, not an estimate: the data do not
identify it. Report the curve, not a single point.

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
proxy_mnar_sensitivity(dat, "y", beta_grid = c(0, 0.4, 0.8, 1.2), m = 10L, seed = 1L)
#>   beta   estimate  std.error    conf.low   conf.high       fmi
#> 1  0.0 -0.2114708 0.09202746 -0.39312678 -0.02981479 0.2294043
#> 2  0.4  0.0105841 0.09663433 -0.17992398  0.20109218 0.2079950
#> 3  0.8  0.1755879 0.10466714 -0.03128688  0.38246275 0.2494541
#> 4  1.2  0.3315785 0.11298771  0.10764312  0.55551396 0.2872016
```
