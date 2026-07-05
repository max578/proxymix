# Optimal action and expected incremental value per unit

Turns the per-unit treatment effect into a next-best action under a
linear value model: treat when the value of the effect exceeds the cost,
i.e. \\d^\*(x) = 1\\\\\text{value} \cdot \tau(x) \> \text{cost}\\\\\\,
with expected incremental value \\\text{value} \cdot \tau(x) -
\text{cost}\\. The standard error of `tau` is propagated to an
action-flip probability – the chance the recommended action would
reverse under sampling noise.

## Usage

``` r
proxy_decide(
  model,
  newdata,
  value,
  cost = 0,
  t1 = 1,
  t0 = 0,
  se_method = c("delta", "mc"),
  ...
)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

- newdata:

  A data frame carrying the covariate columns.

- value:

  Numeric scalar – the value of one unit of outcome.

- cost:

  Numeric scalar – the cost of treating one unit. Default `0`.

- t1, t0:

  The treated and control treatment values. Default `1` and `0`.

- se_method:

  One of `"delta"` (default) or `"mc"`.

- ...:

  Forwarded to
  [`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md).

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `id`, `action`, `expected_value`, `tau`, `se`, `flip_prob`,
`overlap_flag`.

## See also

[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md),
[`proxy_policy_value()`](https://max578.github.io/proxymix/reference/proxy_policy_value.md)

Other decision:
[`fit_uplift()`](https://max578.github.io/proxymix/reference/fit_uplift.md),
[`gmm_cf_mean()`](https://max578.github.io/proxymix/reference/gmm_cf_mean.md),
[`gmm_cf_tail_prob()`](https://max578.github.io/proxymix/reference/gmm_cf_tail_prob.md),
[`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md),
[`proxy_confounding_gap()`](https://max578.github.io/proxymix/reference/proxy_confounding_gap.md),
[`proxy_identification_report()`](https://max578.github.io/proxymix/reference/proxy_identification_report.md),
[`proxy_overlap()`](https://max578.github.io/proxymix/reference/proxy_overlap.md),
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
n <- 400L
x <- stats::rnorm(n)
t <- stats::rbinom(n, 1L, 0.5)
y <- 1 + (x > 0) * t + stats::rnorm(n, sd = 0.5)
dat <- data.frame(y = y, t = t, x = x)
m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                max_iter = 50L, seed = 1L)
proxy_decide(m, data.frame(x = c(-1, 1)), value = 1, cost = 0.2)
#>       id action expected_value        tau         se    flip_prob overlap_flag
#>    <int>  <int>          <num>      <num>      <num>        <num>       <lgcl>
#> 1:     1      0     -0.1090610 0.09093897 0.08741388 1.060815e-01        FALSE
#> 2:     2      1      0.7190012 0.91900122 0.08274963 1.830228e-18        FALSE
```
