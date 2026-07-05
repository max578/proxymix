# A fitted uplift / next-best-action model

The object returned by
[`fit_uplift()`](https://max578.github.io/proxymix/reference/fit_uplift.md):
a joint
[gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) over
`(outcome, treatment, covariates)` together with the column roles, the
identification assumption it will be read under, the outcome type, and
the training sample (retained so that resampling standard errors and
overlap diagnostics are available). The decision verbs dispatch on this
class.

## Usage

``` r
uplift_model(
  fit = NULL,
  roles = list(),
  assume = "ignorability",
  outcome_type = "continuous",
  data = NULL,
  n_train = integer(0),
  treatment_levels = c(0, 1),
  name = "uplift_model",
  metadata = list()
)
```

## Arguments

- fit:

  The joint
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) over
  the stacked `(outcome, treatment, covariates)` coordinates, in that
  column order.

- roles:

  A list with integer indices `outcome`, `treatment`, `covariate` and
  the matching `outcome_name`, `treatment_name`, `covariate_names`.

- assume:

  One of `"ignorability"` (the default uplift assumption) or
  `"latent_confounder"` (the do-operator reading; always flagged).

- outcome_type:

  One of `"continuous"`, `"binary"`, `"count"`.

- data:

  The `n` by `p` numeric training matrix
  `cbind(outcome, treatment, covariates)`.

- n_train:

  Integer scalar – the training sample size.

- treatment_levels:

  Numeric length-2 vector `c(t0, t1)` – the control and treated values
  used at fit time (the observed treatment arms).

- name:

  Human-readable name.

- metadata:

  Optional list of descriptors (e.g. the K-selection trace).

## Value

An S7 object of class `uplift_model`.

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
[`uplift_identification()`](https://max578.github.io/proxymix/reference/uplift_identification.md)
