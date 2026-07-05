# The identification report (an executive one-pager)

The differentiator: a structured audit of what the decision model
identifies, what it assumes, and what it cannot answer. Carries the
estimand, the identification regime and its requirement, the overlap
rate on the supplied population, the confounding-gap magnitude (the
value at risk from unobserved confounding), and the explicit
non-identification of the individual counterfactual law.

## Usage

``` r
proxy_identification_report(model, newdata, t1 = 1, t0 = 0)
```

## Arguments

- model:

  An
  [uplift_model](https://max578.github.io/proxymix/reference/uplift_model.md).

- newdata:

  A data frame carrying the covariate columns – the population the
  report is computed over.

- t1, t0:

  The treated and control treatment values. Default `1` and `0`.

## Value

An S7 object of class `uplift_identification` with a `print` method.

## See also

[`proxy_confounding_gap()`](https://max578.github.io/proxymix/reference/proxy_confounding_gap.md),
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
n <- 600L
x <- stats::rnorm(n)
t <- stats::rbinom(n, 1L, 0.5)
y <- 0.5 * t + x + stats::rnorm(n, sd = 0.5)
dat <- data.frame(y = y, t = t, x = x)
m <- fit_uplift(dat, "y", "t", "x", N = 2L, regime = "sample",
                max_iter = 80L, seed = 1L)
proxy_identification_report(m, data.frame(x = stats::rnorm(100)))
#> == Identification report ==================================
#>   Estimand   : CATE / uplift: E[Y | do(T=t1), X] - E[Y | do(T=t0), X]
#>   Assumption : ignorability
#>                requires (Y(0), Y(1)) independent of T given X.
#>   Regimes    : K = 2   Outcome scale: continuous
#>   Units      : 100
#>   Overlap    : 99.0% of units adequately supported
#>   Confounding gap (value at risk if a latent regime confounds):
#>                mean |Delta| = 0.4842, max |Delta| = 0.572
#>   NOT identified : the individual counterfactual law
#>                    (its variance and tail probabilities).
#> ===========================================================
```
