# Interventional law of a Gaussian mixture (the do-operator)

Reads a fitted Gaussian mixture as a latent-class structural causal
model and returns the interventional distribution of the free
coordinates under `do()` of some coordinates, optionally conditioning on
others.

## Usage

``` r
gmm_intervene(g, do, given = NULL, ridge_eps = 1e-06)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)) in
  `R^p`.

- do:

  A length-`p` numeric vector. Coordinates to intervene on take their
  `do`-value; coordinates not intervened on are `NA`.

- given:

  A length-`p` numeric vector, or `NULL` (the default, meaning no
  conditioning). Coordinates to condition on take their value;
  coordinates not conditioned on are `NA`. A coordinate may not appear
  in both `do` and `given`.

- ridge_eps:

  Tiny ridge added to the conditional covariances for numerical hygiene.
  Set to `0` to disable.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) over the
free coordinates (those `NA` in both `do` and `given`), with weights
re-gated by the `given` evidence only.

## Details

Intervened (`do`) coordinates are *set* inside every component but do
not re-weight the regime gate – this is the graph surgery that
distinguishes \\p(\cdot \mid do(T = t))\\ from \\p(\cdot \mid T = t)\\.
Conditioned (`given`) coordinates re-weight the gate in the usual
Bayesian way. Writing the component prior as \\\pi_k\\, the
within-component conditional mean as \\\mu_k\\, and the given-coordinate
evidence as \\e_k\\, the returned mixture has weights \\\pi_k(given)
\propto \pi_k\\ e_k\\ and per-component parameters from the Schur
conditional on the union of the `do` and `given` coordinates.

For a joint fit over \\(Y, T, X)\\,
`gmm_intervene(fit, do = T = 1, given = X = x)` is the do-response \\p(Y
\mid do(T = 1), X = x)\\; its mean is \\\sum_k \pi_k(x)\\ \mu_k^{y \mid
1, x}\\, the latent-confounder-mode interventional mean of
[proxy_cate](https://max578.github.io/proxymix/reference/proxy_cate.md).

## See also

[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md)

Other operators:
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md),
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
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
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
## Joint (Y, T, X): set T = 1 while conditioning on X = 0.3.
g <- gmm(weights = c(0.5, 0.5),
         means = list(c(0, 0, 0), c(2, 1, 1)),
         covariances = list(diag(3), diag(3)))
gmm_intervene(g, do = c(NA, 1, NA), given = c(NA, NA, 0.3))
#> <intervene(gmm)>: K = 2 components in p = 1 dimensions
#>   [1] w = 0.5498, |mu| = 0.0000, tr(Sigma) = 1.0000
#>   [2] w = 0.4502, |mu| = 2.0000, tr(Sigma) = 1.0000
```
