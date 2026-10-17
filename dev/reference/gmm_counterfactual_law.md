# A per-unit counterfactual law

The return type of
[`gmm_counterfactual()`](https://max578.github.io/proxymix/dev/reference/gmm_counterfactual.md):
a discrete law on `K` atoms (the per-component counterfactual outcomes)
with abduction weights, plus the counterfactual `mean`. All three assume
the unit keeps its structural residual under the intervention, which the
data cannot confirm. Under that assumption the atom spread is
uncertainty about which component generated the unit, reported by
[`gmm_cf_variance()`](https://max578.github.io/proxymix/dev/reference/gmm_cf_variance.md)
and
[`gmm_cf_tail_prob()`](https://max578.github.io/proxymix/dev/reference/gmm_cf_tail_prob.md).

## Usage

``` r
gmm_counterfactual_law(
  atoms = numeric(0),
  weights = numeric(0),
  mean = numeric(0),
  query = integer(0),
  evidence = numeric(0),
  do = numeric(0),
  name = "gmm_counterfactual_law",
  metadata = list()
)
```

## Arguments

- atoms:

  Numeric vector of length `K` – the per-component counterfactual
  outcomes.

- weights:

  Numeric vector of length `K` – the abduction weights
  \\\pi_k(\text{evidence})\\.

- mean:

  Numeric scalar – the counterfactual mean.

- query:

  Integer scalar – the queried coordinate index.

- evidence:

  Numeric vector – the observed unit.

- do:

  Numeric vector – the intervention.

- name:

  Human-readable name.

- metadata:

  Optional list of descriptors.

## Value

An S7 object of class `gmm_counterfactual_law`.

## See also

Other classes:
[`autoplot.gmm_fit`](https://max578.github.io/proxymix/dev/reference/autoplot.gmm_fit.md),
[`glance.gmm_fit`](https://max578.github.io/proxymix/dev/reference/glance.gmm_fit.md),
[`gmm()`](https://max578.github.io/proxymix/dev/reference/gmm.md),
[`gmm_dim()`](https://max578.github.io/proxymix/dev/reference/gmm_dim.md),
[`gmm_fit()`](https://max578.github.io/proxymix/dev/reference/gmm_fit.md),
[`gmm_n_components()`](https://max578.github.io/proxymix/dev/reference/gmm_n_components.md),
[`gmm_target()`](https://max578.github.io/proxymix/dev/reference/gmm_target.md),
[`gmm_weights()`](https://max578.github.io/proxymix/dev/reference/gmm_weights.md),
[`is_proposal()`](https://max578.github.io/proxymix/dev/reference/is_proposal.md),
[`tidy.gmm`](https://max578.github.io/proxymix/dev/reference/tidy.gmm.md)
