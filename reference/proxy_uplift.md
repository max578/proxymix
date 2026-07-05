# Uplift (alias of [`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md) for a binary treatment)

For a binary treatment, the uplift is exactly the conditional average
treatment effect. This is a thin alias of
[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md)
kept for the next-best-action vocabulary.

## Usage

``` r
proxy_uplift(model, newdata, ...)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

- newdata:

  A data frame carrying the covariate columns.

- ...:

  Forwarded to
  [`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
  inside the `"mc"` refits.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
– see
[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md).

## See also

Other decision:
[`fit_uplift()`](https://max578.github.io/proxymix/reference/fit_uplift.md),
[`gmm_cf_mean()`](https://max578.github.io/proxymix/reference/gmm_cf_mean.md),
[`gmm_cf_tail_prob()`](https://max578.github.io/proxymix/reference/gmm_cf_tail_prob.md),
[`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md),
[`proxy_confounding_gap()`](https://max578.github.io/proxymix/reference/proxy_confounding_gap.md),
[`proxy_decide()`](https://max578.github.io/proxymix/reference/proxy_decide.md),
[`proxy_identification_report()`](https://max578.github.io/proxymix/reference/proxy_identification_report.md),
[`proxy_overlap()`](https://max578.github.io/proxymix/reference/proxy_overlap.md),
[`proxy_policy_value()`](https://max578.github.io/proxymix/reference/proxy_policy_value.md),
[`proxy_predict()`](https://max578.github.io/proxymix/reference/proxy_predict.md),
[`proxy_regime_segments()`](https://max578.github.io/proxymix/reference/proxy_regime_segments.md),
[`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/reference/proxy_retrospective_uplift.md),
[`uplift_identification()`](https://max578.github.io/proxymix/reference/uplift_identification.md),
[`uplift_model()`](https://max578.github.io/proxymix/reference/uplift_model.md)

## Examples

``` r
set.seed(1)
dat <- data.frame(y = stats::rnorm(200), t = stats::rbinom(200, 1L, 0.5),
                  x = stats::rnorm(200))
m <- fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment")
proxy_uplift(m, newdata = data.frame(x = 0))
#>       id        tau        se      ci_lo     ci_hi overlap_flag
#>    <int>      <num>     <num>      <num>     <num>       <lgcl>
#> 1:     1 0.02079754 0.1311432 -0.2362383 0.2778334        FALSE
```
