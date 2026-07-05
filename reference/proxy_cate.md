# Heterogeneous treatment effects (CATE / uplift)

Per-unit conditional average treatment effect \\\tau(x) = E\[Y \mid do(T
= t_1), X = x\] - E\[Y \mid do(T = t_0), X = x\]\\, read in closed form
off the fitted joint mixture. Under the model's default `"ignorability"`
assumption this is the contrast of two component-gated conditional
means; under `"latent_confounder"` it is the regime-gated within-class
slope (the do-operator). The two coincide when treatment carries no
information about the regime beyond `X`; their difference is
[`proxy_confounding_gap()`](https://max578.github.io/proxymix/reference/proxy_confounding_gap.md).

## Usage

``` r
proxy_cate(
  model,
  newdata,
  t1 = 1,
  t0 = 0,
  se = TRUE,
  se_method = c("delta", "mc"),
  level = 0.95,
  B = 200L,
  scale = c("link", "response"),
  threshold = 0.5,
  ...
)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

- newdata:

  A data frame carrying the covariate columns.

- t1, t0:

  The treated and control treatment values. Default `1` and `0`.

- se:

  Logical – compute standard errors and confidence intervals.

- se_method:

  One of `"delta"` (closed form, the default) or `"mc"` (resampling).

- level:

  Confidence level for the interval. Default `0.95`.

- B:

  Number of bootstrap refits when `se_method = "mc"`. Default `200`.

- scale:

  One of `"link"` (the latent / continuous scale, the default) or
  `"response"`. For a binary outcome the response scale reports the
  effect on the discretised predictive probability `P(Y > threshold)`;
  for continuous and count outcomes the two scales coincide.

- threshold:

  Decision threshold for the binary discretised predictive. Default
  `0.5`.

- ...:

  Forwarded to
  [`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
  inside the `"mc"` refits.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `id`, `tau`, `se`, `ci_lo`, `ci_hi`, `overlap_flag`.

## Details

The default delta-method standard error is the within-component
prediction variance, holding the regime gate fixed; it reduces to the
ordinary least-squares standard error of the treatment effect at K = 1.
Set `se_method = "mc"` for a resampling standard error that also
reflects gate uncertainty.

## See also

[`proxy_decide()`](https://max578.github.io/proxymix/reference/proxy_decide.md),
[`proxy_confounding_gap()`](https://max578.github.io/proxymix/reference/proxy_confounding_gap.md),
[`proxy_overlap()`](https://max578.github.io/proxymix/reference/proxy_overlap.md)

Other decision:
[`fit_uplift()`](https://max578.github.io/proxymix/reference/fit_uplift.md),
[`gmm_cf_mean()`](https://max578.github.io/proxymix/reference/gmm_cf_mean.md),
[`gmm_cf_tail_prob()`](https://max578.github.io/proxymix/reference/gmm_cf_tail_prob.md),
[`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`proxy_confounding_gap()`](https://max578.github.io/proxymix/reference/proxy_confounding_gap.md),
[`proxy_decide()`](https://max578.github.io/proxymix/reference/proxy_decide.md),
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
y <- 1 + 0.5 * x + (1 + x) * t + stats::rnorm(n, sd = 0.5)
dat <- data.frame(y = y, t = t, x = x)
m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                max_iter = 50L, seed = 1L)
proxy_cate(m, newdata = data.frame(x = c(-1, 0, 1)))
#>       id        tau         se      ci_lo     ci_hi overlap_flag
#>    <int>      <num>      <num>      <num>     <num>       <lgcl>
#> 1:     1 0.02760485 0.08010995 -0.1294078 0.1846175        FALSE
#> 2:     2 1.01893392 0.05416744  0.9127677 1.1251002        FALSE
#> 3:     3 2.01026300 0.07624457  1.8608264 2.1596996        FALSE
```
