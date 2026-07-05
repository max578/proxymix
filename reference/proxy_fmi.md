# Fraction of missing information for a column mean

The share of a column mean's total variance attributable to the missing
data, read from
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md).

## Usage

``` r
proxy_fmi(object, column, method = c("analytic", "rubin"))
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

A named numeric scalar.

## See also

[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md).

Other imputation:
[`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md),
[`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md),
[`gmm_imputation()`](https://max578.github.io/proxymix/reference/gmm_imputation.md),
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`mechanism`](https://max578.github.io/proxymix/reference/mechanism.md),
[`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md),
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)

## Examples

``` r
set.seed(1)
x1 <- rnorm(150); x2 <- x1 + rnorm(150); x2[runif(150) < 0.3] <- NA
imp <- gmm_impute(cbind(x1, x2), N = 1L, m = 10L, seed = 1L)
proxy_fmi(imp, "x2")
#>        x2 
#> 0.1765858 
```
