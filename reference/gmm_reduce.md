# Reduce a Gaussian mixture to fewer components

Collapses a Gaussian mixture to at most `k_max` components. The default
`method = "merge"` is a greedy, moment-preserving pairwise merge: at
each step the cheapest pair of components is replaced by the single
Gaussian that preserves their combined weight, mean and covariance,
until the component budget is met. Because every merge is
moment-preserving, the reduced mixture has the **same global mean and
covariance** as the original, and reducing all the way to a single
component returns the moment-matched Gaussian.

## Usage

``` r
gmm_reduce(
  g,
  k_max,
  method = c("merge", "anneal"),
  cost = c("kl", "cs"),
  draws = 5000L,
  seed = NULL,
  ridge_eps = 0
)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)).

- k_max:

  Maximum number of components to keep (a positive integer). When
  `k_max` is at least the current component count the mixture is
  returned unchanged.

- method:

  `"merge"` (the default, moment-preserving greedy merge) or `"anneal"`
  (the merge refined by an annealed re-fit; never worse than the merge).

- cost:

  The pairwise merge cost: `"kl"` (the Runnalls Kullback-Leibler bound,
  the default) or `"cs"` (the Cauchy-Schwarz divergence).

- draws:

  Number of draws used by the `"anneal"` re-fit. Ignored when
  `method = "merge"`.

- seed:

  Optional integer seed for the `"anneal"` re-fit (the result is
  deterministic given a seed). Ignored when `method = "merge"`.

- ridge_eps:

  Ridge added to each merged covariance. The moment-preserving merge is
  positive-definite by construction, so the default `0` keeps the global
  moments exact; set a small positive value for extra numerical headroom
  in very long reduction chains (at the cost of a tiny moment drift).

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) with at most
`k_max` components.

## Details

The merge cost decides which pair is collapsed first. `cost = "kl"` is
the Kullback-Leibler upper bound of Runnalls (2007), the standard
mixture- reduction criterion; `cost = "cs"` is the closed-form
Cauchy-Schwarz divergence between the two-component sub-mixture and its
merge (the same Gaussian-product identity that underlies
[`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md)),
scaled by the merged mass. Both are non-negative and zero for identical
components.

`method = "anneal"` additionally draws a sample from the mixture and
refits a `k_max`-component proxy by annealed EM (see
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md)),
then returns whichever of the merge and the re-fit has the smaller
Cauchy-Schwarz divergence from the original. The re-fit can improve on
the greedy merge for smooth, over-parameterised mixtures, where a
globally fitted proxy beats any sequence of pairwise merges; the merge
is returned when it is at least as good, so the result is never worse
than `method = "merge"`. Unlike the merge, the re-fit is a Monte Carlo
fit and does not preserve the global moments exactly; raise `draws` for
a closer re-fit.

Reduction is the closing operation of a Gaussian-sum filter: repeated
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md)
/
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md)
steps under a Gaussian-mixture noise or dynamics model multiply the
component count, and `gmm_reduce()` returns it to a fixed budget.

## References

Runnalls, A. R. (2007) Kullback-Leibler approach to Gaussian mixture
reduction. *IEEE Transactions on Aerospace and Electronic Systems*
43(3), 989–999.
[doi:10.1109/TAES.2007.4383588](https://doi.org/10.1109/TAES.2007.4383588)

## See also

Other operators:
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`gmm_missing()`](https://max578.github.io/proxymix/reference/gmm_missing.md),
[`gmm_mix()`](https://max578.github.io/proxymix/reference/gmm_mix.md),
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md),
[`gmm_product()`](https://max578.github.io/proxymix/reference/gmm_product.md)

## Examples

``` r
## A six-component mixture with three near-duplicate pairs.
g <- gmm(
  weights = rep(1 / 6, 6),
  means = list(c(-4, 0), c(-4, 0.1), c(4, 0), c(4.1, 0), c(0, 5), c(0, 5.1)),
  covariances = rep(list(diag(2)), 6)
)
gmm_reduce(g, k_max = 3L)
#> <reduce[3](gmm)>: K = 3 components in p = 2 dimensions
#>   [1] w = 0.3333, |mu| = 4.0003, tr(Sigma) = 2.0025
#>   [2] w = 0.3333, |mu| = 4.0500, tr(Sigma) = 2.0025
#>   [3] w = 0.3333, |mu| = 5.0500, tr(Sigma) = 2.0025
```
