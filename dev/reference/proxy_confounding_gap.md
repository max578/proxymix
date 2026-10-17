# Confounding gap: the sensitivity of the effect to the latent regime

Per-unit difference between the ignorability-mode and do-mode effects,
\\\Delta(x) = \tau\_{\mathrm{obs}}(x) - \tau\_{\mathrm{do}}(x)\\. Under
ignorability the two coincide and \\\Delta \equiv 0\\; a non-zero gap is
a **sensitivity signal** – how much the estimated effect would move if a
fitted regime confounded treatment and outcome beyond `X` – not a
correction the data licenses.

## Usage

``` r
proxy_confounding_gap(model, newdata, t1 = NULL, t0 = NULL)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/dev/reference/uplift_model.md).

- newdata:

  A data frame carrying the covariate columns.

- t1, t0:

  The treated and control treatment values. Default the treatment levels
  observed at fit time (`model@treatment_levels`); a value matching
  neither observed level aborts.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with columns `id`, `tau_obs`, `tau_do`, `gap`, `overlap_flag`.

## Details

When each regime holds a single treatment arm, as often happens with a
binary treatment, the within-regime effect is zero, the gap equals the
whole estimated effect, and a warning of class
`proxymix_constant_treatment` is raised.

## See also

[`proxy_cate()`](https://max578.github.io/proxymix/dev/reference/proxy_cate.md),
[`proxy_identification_report()`](https://max578.github.io/proxymix/dev/reference/proxy_identification_report.md)

Other decision:
[`fit_uplift()`](https://max578.github.io/proxymix/dev/reference/fit_uplift.md),
[`gmm_cf_mean()`](https://max578.github.io/proxymix/dev/reference/gmm_cf_mean.md),
[`gmm_cf_tail_prob()`](https://max578.github.io/proxymix/dev/reference/gmm_cf_tail_prob.md),
[`gmm_cf_variance()`](https://max578.github.io/proxymix/dev/reference/gmm_cf_variance.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/dev/reference/gmm_counterfactual.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/dev/reference/gmm_intervene.md),
[`proxy_cate()`](https://max578.github.io/proxymix/dev/reference/proxy_cate.md),
[`proxy_decide()`](https://max578.github.io/proxymix/dev/reference/proxy_decide.md),
[`proxy_identification_report()`](https://max578.github.io/proxymix/dev/reference/proxy_identification_report.md),
[`proxy_overlap()`](https://max578.github.io/proxymix/dev/reference/proxy_overlap.md),
[`proxy_policy_value()`](https://max578.github.io/proxymix/dev/reference/proxy_policy_value.md),
[`proxy_predict()`](https://max578.github.io/proxymix/dev/reference/proxy_predict.md),
[`proxy_regime_segments()`](https://max578.github.io/proxymix/dev/reference/proxy_regime_segments.md),
[`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/dev/reference/proxy_retrospective_uplift.md),
[`proxy_uplift()`](https://max578.github.io/proxymix/dev/reference/proxy_uplift.md),
[`uplift_identification()`](https://max578.github.io/proxymix/dev/reference/uplift_identification.md),
[`uplift_model()`](https://max578.github.io/proxymix/dev/reference/uplift_model.md)

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
proxy_confounding_gap(m, data.frame(x = c(-1, 0, 1)))
#>       id  tau_obs    tau_do      gap overlap_flag
#>    <int>    <num>     <num>    <num>       <lgcl>
#> 1:     1 2.908393 0.5072191 2.401174        FALSE
#> 2:     2 2.868516 0.5105320 2.357984        FALSE
#> 3:     3 2.801447 0.5132839 2.288164        FALSE
```
