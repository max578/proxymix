# Retrospective (counterfactual-mean) uplift for observed units

For each observed unit `(y, t, x)`, the counterfactual-mean uplift of
moving from `t0` to `t1`, \\E\[Y\_{t_1} \mid y, t, x\] - E\[Y\_{t_0}
\mid y, t, x\]\\, computed by
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md).
Unlike
[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md),
the abduction gate uses the observed outcome `y` as well, sharpening the
per-unit estimate. Only the counterfactual *mean* is identified; the
spread is not (see
[`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md)).

## Usage

``` r
proxy_retrospective_uplift(model, observed, t1 = 1, t0 = 0)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

- observed:

  A data frame carrying the outcome, treatment, and covariate columns of
  the observed units.

- t1, t0:

  The treated and control treatment values. Default `1` and `0`.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `id`, `y_obs`, `t_obs`, `cf_mean_t1`, `retro_uplift`.

## See also

[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md)

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
[`proxy_uplift()`](https://max578.github.io/proxymix/reference/proxy_uplift.md),
[`uplift_identification()`](https://max578.github.io/proxymix/reference/uplift_identification.md),
[`uplift_model()`](https://max578.github.io/proxymix/reference/uplift_model.md)

## Examples

``` r
set.seed(1)
n <- 400L
x <- stats::rnorm(n)
t <- stats::rbinom(n, 1L, 0.5)
y <- 1 + (0.5 + x) * t + stats::rnorm(n, sd = 0.5)
dat <- data.frame(y = y, t = t, x = x)
m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                max_iter = 80L, seed = 1L)
proxy_retrospective_uplift(m, observed = dat[1:5, ])
#>       id     y_obs t_obs cf_mean_t1 retro_uplift
#>    <int>     <num> <num>      <num>        <num>
#> 1:     1 0.7030127     1  0.7030127            0
#> 2:     2 1.7512123     0  1.7512123            0
#> 3:     3 0.9285252     1  0.9285252            0
#> 4:     4 3.3663765     1  3.3663765            0
#> 5:     5 0.9316633     0  0.9316633            0
```
