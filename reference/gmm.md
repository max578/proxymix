# A Gaussian mixture

Lightweight S7 class representing an `N`-component multivariate Gaussian
mixture on \\\mathbb{R}^p\\. Use `gmm()` to construct,
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md) /
[`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md) to
evaluate or sample, and
[`gmm_marginalise()`](https://max578.github.io/proxymix/reference/gmm_marginalise.md)
/
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md)
for closed-form operations.

## Usage

``` r
gmm(
  weights = numeric(0),
  means = list(),
  covariances = list(),
  name = "gmm",
  metadata = list()
)
```

## Arguments

- weights:

  Numeric vector of length `K`, non-negative, summing to one.

- means:

  List of length `K`, each element a length-`p` numeric vector.

- covariances:

  List of length `K`, each element a `p`-by-`p` symmetric
  positive-definite numeric matrix.

- name:

  Optional human-readable name.

- metadata:

  Optional list of arbitrary metadata (regime tags, diagnostic
  snapshots, etc.).

## Value

An S7 object inheriting from `gmm`.

## See also

Other classes:
[`autoplot.gmm_fit`](https://max578.github.io/proxymix/reference/autoplot.gmm_fit.md),
[`glance.gmm_fit`](https://max578.github.io/proxymix/reference/glance.gmm_fit.md),
[`gmm_counterfactual_law()`](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md),
[`gmm_dim()`](https://max578.github.io/proxymix/reference/gmm_dim.md),
[`gmm_fit()`](https://max578.github.io/proxymix/reference/gmm_fit.md),
[`gmm_n_components()`](https://max578.github.io/proxymix/reference/gmm_n_components.md),
[`gmm_target()`](https://max578.github.io/proxymix/reference/gmm_target.md),
[`gmm_weights()`](https://max578.github.io/proxymix/reference/gmm_weights.md),
[`is_proposal()`](https://max578.github.io/proxymix/reference/is_proposal.md),
[`tidy.gmm`](https://max578.github.io/proxymix/reference/tidy.gmm.md)

## Examples

``` r
g <- gmm(
  weights = c(0.4, 0.6),
  means = list(c(-1, 0), c(1, 0)),
  covariances = list(diag(2), diag(2))
)
g
#> <gmm>: K = 2 components in p = 2 dimensions
#>   [1] w = 0.4000, |mu| = 1.0000, tr(Sigma) = 2.0000
#>   [2] w = 0.6000, |mu| = 1.0000, tr(Sigma) = 2.0000
```
