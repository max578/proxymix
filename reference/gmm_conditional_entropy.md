# Conditional predictive entropy of a Gaussian mixture

Returns the differential entropy of the conditional mixture \\g\_{Y \mid
X = x}\\ obtained from
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md)
– the predictive uncertainty of the target coordinates given the
conditioned ones. The order-2 Renyi entropy is closed-form;
`order = "shannon"` falls back to Monte Carlo. Multiple conditioning
configurations are evaluated row-by-row.

## Usage

``` r
gmm_conditional_entropy(
  g,
  given,
  order = c("renyi2", "shannon"),
  n_mc = 5000L,
  seed = NULL
)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  joint mixture.

- given:

  Either a numeric vector with one entry per coordinate, or a matrix
  whose rows are such vectors. `NA` marks a target (kept) coordinate; a
  numeric value conditions on that coordinate (the
  [`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md)
  convention).

- order:

  `"renyi2"` (closed-form, the default) or `"shannon"`.

- n_mc, seed:

  Passed to
  [`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md)
  for `order = "shannon"`.

## Value

A numeric scalar for a single configuration, or a numeric vector with
one entropy per row of `given`.

## See also

[`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md)

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md),
[`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md),
[`ess_trace()`](https://max578.github.io/proxymix/reference/ess_trace.md),
[`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md),
[`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md),
[`gmm_evidence()`](https://max578.github.io/proxymix/reference/gmm_evidence.md),
[`gmm_fit_ensemble()`](https://max578.github.io/proxymix/reference/gmm_fit_ensemble.md),
[`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md),
[`gmm_independence_graph()`](https://max578.github.io/proxymix/reference/gmm_independence_graph.md),
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
## Joint over (Y, X); predictive entropy of Y at several X values.
s <- matrix(c(2, 0.8, 0.8, 1), 2, 2)
g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(s))
gmm_conditional_entropy(g, given = rbind(c(NA, 0), c(NA, 1)))
#> [1] 1.419254 1.419254
```
