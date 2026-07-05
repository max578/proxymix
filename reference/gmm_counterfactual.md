# Counterfactual law of one unit (abduction, action, prediction)

Computes the per-unit counterfactual of a `query` coordinate (the
outcome) for an observed unit, under a `do()` intervention – Pearl's
third rung on the latent-class structural causal model read off a fitted
Gaussian mixture.

## Usage

``` r
gmm_counterfactual(g, evidence, do, query, ridge_eps = 1e-06)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)) in
  `R^p`.

- evidence:

  A length-`p` numeric vector of the observed unit. Unobserved
  coordinates are `NA`; the `query` coordinate must be observed.

- do:

  A length-`p` numeric vector. Intervened coordinates take their
  counterfactual value; the rest are `NA`. Intervened coordinates must
  be observed in `evidence`.

- query:

  A single integer coordinate index in `seq_len(p)` – the outcome whose
  counterfactual is sought.

- ridge_eps:

  Tiny ridge added to the conditioning covariances for numerical
  hygiene.

## Value

A `gmm_counterfactual_law` object carrying the `K` atoms, their
abduction weights, and the identified counterfactual `mean`.

## Details

The three steps are closed form. *Abduction* recovers the regime
posterior \\\pi_k(\text{evidence})\\ of the observed unit and, within
each component, its structural residual. *Action* sets the `do`
coordinates. *Prediction* re-evaluates the query coordinate. For a
binary treatment the result is a discrete law on `K` atoms, \$\$Y\_{t'}
\mid (y, t, x) = \sum_k \pi_k(y, t, x)\\ \delta\\\big(y +
\beta_k^{T}(t' - t)\big),\$\$ where \\\beta_k^T\\ is component `k`'s
within-class treatment slope. Only the **mean** of this law is
identified; its spread reflects regime uncertainty, not the
(unidentified) cross-world coupling, so the variance and tail accessors
refuse to answer (see
[`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md)).

## See also

[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md),
[`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/reference/proxy_retrospective_uplift.md)

Other operators:
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md),
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`gmm_missing()`](https://max578.github.io/proxymix/reference/gmm_missing.md),
[`gmm_mix()`](https://max578.github.io/proxymix/reference/gmm_mix.md),
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md),
[`gmm_product()`](https://max578.github.io/proxymix/reference/gmm_product.md),
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)

Other decision:
[`fit_uplift()`](https://max578.github.io/proxymix/reference/fit_uplift.md),
[`gmm_cf_mean()`](https://max578.github.io/proxymix/reference/gmm_cf_mean.md),
[`gmm_cf_tail_prob()`](https://max578.github.io/proxymix/reference/gmm_cf_tail_prob.md),
[`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md),
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
## Observed unit (y = 1.2, t = 0, x = 0.5); imagine t = 1.
g <- gmm(weights = c(0.6, 0.4),
         means = list(c(0, 0, 0), c(2, 1, 1)),
         covariances = list(diag(3), diag(3)))
cf <- gmm_counterfactual(g, evidence = c(1.2, 0, 0.5),
                         do = c(NA, 1, NA), query = 1L)
cf@mean
#> [1] 1.2
```
