# Affine pushforward of a Gaussian mixture

Returns the (closed-form) distribution of \\Y = A X + b + \epsilon\\
when \\X \sim g\\ is a Gaussian mixture and \\\epsilon \sim
\mathcal{N}(0, R)\\ is independent additive Gaussian noise.

## Usage

``` r
gmm_affine(g, A, b = 0, noise_cov = NULL, ridge_eps = 1e-06)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)) in
  `R^p`.

- A:

  An `m` by `p` numeric matrix.

- b:

  Numeric scalar or length-`m` vector. Default `0`.

- noise_cov:

  `m` by `m` SPD numeric matrix, or `NULL` (treated as the zero matrix —
  a deterministic channel).

- ridge_eps:

  Tiny ridge added to the output covariances for numerical hygiene. Set
  to `0` to disable.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) in `R^m`
with the same number of components and the same weights as `g`.

## Details

For each component `k`, the pushed-forward parameters are \$\$\mu'\_k =
A \mu_k + b, \qquad \Sigma'\_k = A \Sigma_k A^\top + R,\$\$ and the
mixture weights are unchanged. This is the finite-mixture analogue of a
Kalman-style predict step.

The channel is required to be **affine** in `x` and the noise is
required to be Gaussian. Non-linear channels are not closed form and are
not silently approximated: push samples through the map instead
([`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md) then
the transform) and refit with
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md)
when a mixture of the image is needed.

## See also

Other operators:
[`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md),
[`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md),
[`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md),
[`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md),
[`gmm_missing()`](https://max578.github.io/proxymix/reference/gmm_missing.md),
[`gmm_mix()`](https://max578.github.io/proxymix/reference/gmm_mix.md),
[`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md),
[`gmm_product()`](https://max578.github.io/proxymix/reference/gmm_product.md),
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)

## Examples

``` r
g <- gmm(weights = c(0.5, 0.5),
         means = list(c(-1, 0), c(1, 0)),
         covariances = list(diag(2), diag(2)))
A <- matrix(c(1, 0, 0, 1, 1, -1), nrow = 3L, byrow = TRUE)
gmm_affine(g, A, b = c(0, 0, 0), noise_cov = 0.01 * diag(3))
#> <affine(gmm)>: K = 2 components in p = 3 dimensions
#>   [1] w = 0.5000, |mu| = 1.4142, tr(Sigma) = 4.0300
#>   [2] w = 0.5000, |mu| = 1.4142, tr(Sigma) = 4.0300
```
