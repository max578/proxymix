# Off-line value of a targeting policy

Estimates the expected value of deploying a per-unit targeting policy,
\\V(d) = E_X\[\\\text{value}\cdot E\[Y \mid do(T = d(X)), X\] -
\text{cost} \cdot d(X)\\\]\\, from the fitted model alone – no live A/B
test. Units that fail the overlap diagnostic are excluded by default and
their count is reported, never silently dropped.

## Usage

``` r
proxy_policy_value(
  model,
  newdata,
  policy,
  value,
  cost = 0,
  t1 = 1,
  t0 = 0,
  exclude_low_overlap = TRUE
)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

- newdata:

  A data frame carrying the covariate columns – the population the
  policy would be deployed on.

- policy:

  A per-unit action specification: a 0/1 vector of length
  `nrow(newdata)`, a function of the
  [`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md)
  table returning actions, or one of the strings `"all"`, `"none"`,
  `"optimal"`.

- value:

  Numeric scalar – the value of one unit of outcome.

- cost:

  Numeric scalar – the cost of treating one unit. Default `0`.

- t1, t0:

  The treated and control treatment values. Default `1` and `0`.

- exclude_low_overlap:

  Logical – drop overlap-flagged units from the average (and report the
  count). Default `TRUE`.

## Value

A one-row
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `policy_value`, `n_used`, `n_excluded`, `n_treated`.

## See also

[`proxy_decide()`](https://max578.github.io/proxymix/reference/proxy_decide.md),
[`proxy_overlap()`](https://max578.github.io/proxymix/reference/proxy_overlap.md)

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
[`proxy_predict()`](https://max578.github.io/proxymix/reference/proxy_predict.md),
[`proxy_regime_segments()`](https://max578.github.io/proxymix/reference/proxy_regime_segments.md),
[`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/reference/proxy_retrospective_uplift.md),
[`proxy_uplift()`](https://max578.github.io/proxymix/reference/proxy_uplift.md),
[`uplift_identification()`](https://max578.github.io/proxymix/reference/uplift_identification.md),
[`uplift_model()`](https://max578.github.io/proxymix/reference/uplift_model.md)

## Examples

``` r
set.seed(1)
n <- 600L
x <- stats::rnorm(n)
t <- stats::rbinom(n, 1L, 0.5)
y <- 1 + (0.4 + x) * t + stats::rnorm(n, sd = 0.5)
dat <- data.frame(y = y, t = t, x = x)
m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                max_iter = 80L, seed = 1L)
nd <- data.frame(x = stats::rnorm(200))
proxy_policy_value(m, nd, policy = "optimal", value = 1, cost = 0.3)
#> Excluded 1 low-overlap unit from the policy value.
#>    policy_value n_used n_excluded n_treated
#>           <num>  <int>      <int>     <int>
#> 1:     1.448247    199          1       105
```
