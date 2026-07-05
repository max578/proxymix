# Identification-report object

The structured return type of
[`proxy_identification_report()`](https://max578.github.io/proxymix/reference/proxy_identification_report.md).
Print it for the executive one-pager.

## Usage

``` r
uplift_identification(
  estimand = character(0),
  assume = character(0),
  n_units = integer(0),
  overlap_pct = numeric(0),
  confounding_gap_mean = numeric(0),
  confounding_gap_max = numeric(0),
  K = integer(0),
  outcome_type = "continuous"
)
```

## Arguments

- estimand:

  Character – the target estimand.

- assume:

  Character – the identification regime.

- n_units:

  Integer – units the report covers.

- overlap_pct:

  Numeric – percentage of units with adequate overlap.

- confounding_gap_mean, confounding_gap_max:

  Numeric – mean and maximum absolute confounding gap over the
  population.

- K:

  Integer – the number of fitted regimes.

- outcome_type:

  Character – the outcome scale.

## Value

An S7 object of class `uplift_identification`.

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
[`proxy_regime_segments()`](https://max578.github.io/proxymix/reference/proxy_regime_segments.md),
[`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/reference/proxy_retrospective_uplift.md),
[`proxy_uplift()`](https://max578.github.io/proxymix/reference/proxy_uplift.md),
[`uplift_model()`](https://max578.github.io/proxymix/reference/uplift_model.md)
