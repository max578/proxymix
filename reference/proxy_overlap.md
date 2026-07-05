# Per-unit overlap / positivity diagnostic

Flags units whose `(treatment, covariate)` configuration is poorly
covered by the fitted joint – the proxy's mass coverage is the
positivity diagnostic. For each treatment arm the squared Mahalanobis
distance to the nearest regime centre is converted to an upper-tail
chi-square coverage probability; the reported `coverage` is the minimum
across arms, since the treatment effect needs both arms supported. Units
below `floor` are flagged and excluded from
[`proxy_policy_value()`](https://max578.github.io/proxymix/reference/proxy_policy_value.md)
by default.

## Usage

``` r
proxy_overlap(model, newdata, t1 = 1, t0 = 0, floor = 0.01)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

- newdata:

  A data frame carrying the covariate columns.

- t1, t0:

  The treated and control treatment values. Default `1` and `0`.

- floor:

  Coverage probability below which a unit is flagged. Default `0.01`.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `id`, `coverage`, `overlap_flag`.

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
[`proxy_policy_value()`](https://max578.github.io/proxymix/reference/proxy_policy_value.md),
[`proxy_predict()`](https://max578.github.io/proxymix/reference/proxy_predict.md),
[`proxy_regime_segments()`](https://max578.github.io/proxymix/reference/proxy_regime_segments.md),
[`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/reference/proxy_retrospective_uplift.md),
[`proxy_uplift()`](https://max578.github.io/proxymix/reference/proxy_uplift.md),
[`uplift_identification()`](https://max578.github.io/proxymix/reference/uplift_identification.md),
[`uplift_model()`](https://max578.github.io/proxymix/reference/uplift_model.md)

## Examples

``` r
set.seed(1)
dat <- data.frame(y = stats::rnorm(200), t = stats::rbinom(200, 1L, 0.5),
                  x = stats::rnorm(200))
m <- fit_uplift(dat, "y", "t", "x", N = 1L, regime = "moment")
proxy_overlap(m, newdata = data.frame(x = c(0, 8)))
#>       id     coverage overlap_flag
#>    <int>        <num>       <lgcl>
#> 1:     1 5.887126e-01        FALSE
#> 2:     2 4.893549e-13         TRUE
```
