# The fitted regimes as an interpretable segment table

Exposes the `K` mixture components as decision segments: each regime's
prevalence (`weight`), its within-segment treatment effect (the
within-class treatment slope), its residual standard deviation, and its
covariate centre. This is the interpretable by-product the closed-form
reading gives for free.

## Usage

``` r
proxy_regime_segments(model, t1 = NULL, t0 = NULL)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

- t1, t0:

  The treated and control treatment values used to scale the
  within-segment effect. Default the treatment levels observed at fit
  time (`model@treatment_levels`).

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `regime`, `weight`, `effect`, `sigma`, and one column per
covariate centre.

## Details

When each regime holds a single treatment arm, as often happens with a
binary treatment, every within-segment effect is zero and a warning of
class `proxymix_constant_treatment` is raised.

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
[`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/reference/proxy_retrospective_uplift.md),
[`proxy_uplift()`](https://max578.github.io/proxymix/reference/proxy_uplift.md),
[`uplift_identification()`](https://max578.github.io/proxymix/reference/uplift_identification.md),
[`uplift_model()`](https://max578.github.io/proxymix/reference/uplift_model.md)

## Examples

``` r
set.seed(1)
n <- 600L
u <- stats::rbinom(n, 1L, 0.5) # unobserved group
x <- stats::rnorm(n)
t <- stats::rbinom(n, 1L, 0.3 + 0.4 * u)
y <- 1 + 5 * u + 0.5 * t + x + stats::rnorm(n, sd = 0.5)
dat <- data.frame(y = y, t = t, x = x)
m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                max_iter = 80L, seed = 1L)
proxy_regime_segments(m)
#>    regime    weight    effect     sigma           x
#>     <int>     <num>     <num>     <num>       <num>
#> 1:      1 0.5233333 0.5685978 0.5184802  0.01701770
#> 2:      2 0.4766667 0.4458673 0.5168527 -0.09606112
```
