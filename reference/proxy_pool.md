# Pool a column mean across imputations

Pools the mean of one column over the `m` completed datasets in a
[gmm_imputation](https://max578.github.io/proxymix/reference/gmm_imputation.md),
returning the estimate with a standard error, degrees of freedom,
confidence interval, and fraction of missing information. The default
`method = "analytic"` computes the between-imputation variance in closed
form – the exact \\m \to \infty\\ limit, with no Monte-Carlo noise –
from the mixture conditional, and splits the total variance into
complete-data, imputation, and parameter parts. `method = "rubin"`
instead applies Rubin's rules to the drawn completions (useful as a
check).

## Usage

``` r
proxy_pool(object, column, method = c("analytic", "rubin"))
```

## Arguments

- object:

  A
  [gmm_imputation](https://max578.github.io/proxymix/reference/gmm_imputation.md).

- column:

  Name of a single numeric column whose mean is pooled.

- method:

  `"analytic"` (the default) for the closed-form pooling, or `"rubin"`
  for Rubin's rules over the drawn completions.

## Value

A one-row data frame: `term`, `estimate`, `std.error`, `statistic`,
`df`, `conf.low`, `conf.high`, `fmi`.

## Details

For a regression or any other model estimand, do not pool here: convert
the imputations to a `mice` object with
[`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md)
and pool with
[`mice::pool()`](https://amices.org/mice/reference/pool.html), which is
the established workflow and reports the same diagnostics.

## See also

[`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md),
[`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md) to
pool models with
[`mice::pool()`](https://amices.org/mice/reference/pool.html).

Other imputation:
[`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md),
[`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md),
[`gmm_imputation()`](https://max578.github.io/proxymix/reference/gmm_imputation.md),
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`mechanism`](https://max578.github.io/proxymix/reference/mechanism.md),
[`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md),
[`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md)

## Examples

``` r
set.seed(1)
x1 <- rnorm(150); x2 <- x1 + rnorm(150)
x2[runif(150) < 0.3] <- NA
imp <- gmm_impute(cbind(x1, x2), N = 1L, m = 10L, seed = 1L)
proxy_pool(imp, "x2")                       # analytic column mean
#>   term   estimate std.error statistic       df   conf.low conf.high       fmi
#> 1   x2 0.03396067 0.1171525 0.2898843 3583.756 -0.1957316 0.2636529 0.1765858
```
