# The counterfactual variance

Returns the variance of a
[gmm_counterfactual_law](https://max578.github.io/proxymix/dev/reference/gmm_counterfactual_law.md),
the weighted spread of its atoms around the counterfactual mean. Like
the law itself, it assumes the unit keeps its structural residual under
the intervention. Under that assumption the only uncertainty left is
which component generated the unit; if the residual were redrawn, the
same mixture would give a different variance.

## Usage

``` r
gmm_cf_variance(x)
```

## Arguments

- x:

  A `gmm_counterfactual_law`.

## Value

Numeric scalar. Zero when the mixture has one component.

## See also

Other decision:
[`fit_uplift()`](https://max578.github.io/proxymix/dev/reference/fit_uplift.md),
[`gmm_cf_mean()`](https://max578.github.io/proxymix/dev/reference/gmm_cf_mean.md),
[`gmm_cf_tail_prob()`](https://max578.github.io/proxymix/dev/reference/gmm_cf_tail_prob.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/dev/reference/gmm_counterfactual.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/dev/reference/gmm_intervene.md),
[`proxy_cate()`](https://max578.github.io/proxymix/dev/reference/proxy_cate.md),
[`proxy_confounding_gap()`](https://max578.github.io/proxymix/dev/reference/proxy_confounding_gap.md),
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
g <- gmm(weights = c(0.6, 0.4),
         means = list(c(0, 0, 0), c(2, 1, 1)),
         covariances = list(diag(3), diag(3)))
cf <- gmm_counterfactual(g, evidence = c(1.2, 0, 0.5),
                         do = c(NA, 1, NA), query = 1L)
gmm_cf_variance(cf)
#> [1] 0
```
