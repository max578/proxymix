# Convert imputations to a mice multiply-imputed dataset

Packages a
[gmm_imputation](https://max578.github.io/proxymix/reference/gmm_imputation.md)
as a [`mice::mids`](https://amices.org/mice/reference/mids.html) object
so that an arbitrary model estimand can be fitted and pooled with the
established mice workflow. The joint Gaussian-mixture imputations –
including the multimodal and heteroscedastic shapes a univariate imputer
cannot produce – flow through unchanged; mice supplies
[`with()`](https://rdrr.io/r/base/with.html),
[`mice::pool()`](https://amices.org/mice/reference/pool.html), and the
pooled diagnostics.

## Usage

``` r
as_mids(object)
```

## Arguments

- object:

  A
  [gmm_imputation](https://max578.github.io/proxymix/reference/gmm_imputation.md).

## Value

A [`mice::mids`](https://amices.org/mice/reference/mids.html) object
with `m` imputations.

## See also

[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)
for the closed-form column-mean pooling.

Other imputation:
[`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md),
[`gmm_imputation()`](https://max578.github.io/proxymix/reference/gmm_imputation.md),
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`mechanism`](https://max578.github.io/proxymix/reference/mechanism.md),
[`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md),
[`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md),
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)

## Examples

``` r
set.seed(1)
x1 <- rnorm(150); x2 <- x1 + rnorm(150)
x2[runif(150) < 0.3] <- NA
imp <- gmm_impute(cbind(x1, x2), N = 1L, m = 10L, seed = 1L)
if (requireNamespace("mice", quietly = TRUE)) {
  fit <- with(as_mids(imp), lm(x2 ~ x1))
  summary(mice::pool(fit))
}
#>          term   estimate  std.error statistic        df      p.value
#> 1 (Intercept) 0.01909702 0.09597906 0.1989707  71.05456 8.428539e-01
#> 2          x1 0.89798425 0.09967653 9.0089838 111.29150 6.668483e-15
```
