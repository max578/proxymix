# Refused: a tail probability of an individual counterfactual law

Like
[`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md),
a per-unit counterfactual tail probability \\P(Y\_{t'} \> c \mid y, t,
x)\\ is not identified – it is a functional of the unidentified
counterfactual law, not of its identified mean. This accessor refuses
rather than mislead.

## Usage

``` r
gmm_cf_tail_prob(x, threshold)
```

## Arguments

- x:

  A
  [gmm_counterfactual_law](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md).

- threshold:

  Numeric scalar `c`.

## Value

Never returns; always raises a `proxymix_not_identified` error.

## See also

Other decision:
[`fit_uplift()`](https://max578.github.io/proxymix/reference/fit_uplift.md),
[`gmm_cf_mean()`](https://max578.github.io/proxymix/reference/gmm_cf_mean.md),
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
[`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/reference/proxy_retrospective_uplift.md),
[`proxy_uplift()`](https://max578.github.io/proxymix/reference/proxy_uplift.md),
[`uplift_identification()`](https://max578.github.io/proxymix/reference/uplift_identification.md),
[`uplift_model()`](https://max578.github.io/proxymix/reference/uplift_model.md)

## Examples

``` r
g <- gmm(weights = 1, means = list(c(0, 0, 0)),
         covariances = list(diag(3)))
cf <- gmm_counterfactual(g, evidence = c(1, 0, 0.2),
                         do = c(NA, 1, NA), query = 1L)
try(gmm_cf_tail_prob(cf, threshold = 2))
#> Error in gmm_cf_tail_prob(cf, threshold = 2) : 
#>   Individual counterfactual tail probabilities are not identified.
#> ℹ They are functionals of the unidentified counterfactual law, not its mean.
#> ℹ Only the counterfactual mean is identified; see `gmm_cf_mean()`.
```
