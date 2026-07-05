# Predicted outcome under a treatment (the seeing rung)

Per-unit predicted outcome \\E\[Y \mid do(T = t), X = x\]\\ – the first
rung of the ladder, risk / response scoring. Under `"ignorability"` this
is the component-gated conditional mean; under `"latent_confounder"` it
is the regime-gated interventional mean. For a binary outcome with
`scale = "response"` the prediction is the discretised predictive
probability `P(Y > threshold)`.

## Usage

``` r
proxy_predict(
  model,
  newdata,
  t,
  scale = c("link", "response"),
  threshold = 0.5
)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

- newdata:

  A data frame carrying the covariate columns.

- t:

  The treatment value to predict the outcome under.

- scale:

  One of `"link"` (default) or `"response"`.

- threshold:

  Decision threshold for the binary discretised predictive. Default
  `0.5`.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `id` and `prediction`.

## See also

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
proxy_predict(m, data.frame(x = c(-1, 0, 1)), t = 1)
#>       id   prediction
#>    <int>        <num>
#> 1:     1 -0.002176397
#> 2:     2  0.045330305
#> 3:     3  0.092837007
```
