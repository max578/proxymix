# Extract completed datasets from a `gmm_imputation`

Extract completed datasets from a `gmm_imputation`

## Usage

``` r
gmm_complete(object, which = 1L)
```

## Arguments

- object:

  A
  [gmm_imputation](https://max578.github.io/proxymix/reference/gmm_imputation.md).

- which:

  Either an integer vector of imputation indices, or `"all"` for every
  completion. Default `1L`.

## Value

When `which` selects one completion, a single completed dataset (matrix,
or data frame if the input was one); otherwise a list of them.

## See also

Other imputation:
[`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md),
[`gmm_imputation()`](https://max578.github.io/proxymix/reference/gmm_imputation.md),
[`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md),
[`mechanism`](https://max578.github.io/proxymix/reference/mechanism.md),
[`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md),
[`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md),
[`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)

## Examples

``` r
x1 <- rnorm(60)
x2 <- x1 + rnorm(60, sd = 0.4)
x2[sample(60, 10)] <- NA
imp <- gmm_impute(cbind(x1, x2), N = 1L, m = 10L, seed = 1L)
first <- gmm_complete(imp, which = 1L)   # one completed dataset
anyNA(first)
#> [1] FALSE
all_completions <- gmm_complete(imp, which = "all")   # a list of ten
length(all_completions)
#> [1] 10
```
