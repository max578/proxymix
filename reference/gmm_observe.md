# Bayesian update of a Gaussian mixture on a noisy linear observation

Conditions a Gaussian mixture `g` on a single noisy linear observation
\\y = A x + b + \epsilon\\, \\\epsilon \sim \mathcal{N}(0, R)\\. Per
component, applies the Kalman gain \$\$K_k = \Sigma_k A^\top S_k^{-1},
\qquad S_k = A \Sigma_k A^\top + R,\$\$ and updates \$\$\mu'\_k =
\mu_k + K_k (y - A \mu_k - b), \qquad \Sigma'\_k = (I - K_k A)
\Sigma_k.\$\$ Component weights are multiplied by the marginal evidence
\\\pi_k \mathcal{N}(y; A \mu_k + b, S_k)\\ and renormalised. This is the
finite-mixture analogue of a Kalman update step.

## Usage

``` r
gmm_observe(g, A, y, noise_cov, b = 0, ridge_eps = 1e-06)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)) in
  `R^p`.

- A:

  An `m` by `p` numeric matrix.

- y:

  A length-`m` numeric vector (the observation).

- noise_cov:

  An `m` by `m` SPD numeric matrix (the observation noise covariance
  `R`). Required.

- b:

  Numeric scalar or length-`m` vector. Default `0`.

- ridge_eps:

  Tiny ridge added to updated covariances for numerical hygiene. Set to
  `0` to disable.

## Value

A [gmm](https://max578.github.io/proxymix/reference/gmm.md) in `R^p`
with the same number of components and the reweighted component weights.

## Details

If the marginal evidence vanishes at every component (e.g. `y` is many
standard deviations from every component), the function issues a warning
and returns `g` unchanged with `metadata$gmm_observe_no_update = TRUE`.

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
[`gmm_product()`](https://max578.github.io/proxymix/reference/gmm_product.md),
[`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)

## Examples

``` r
g <- gmm(weights = c(0.5, 0.5),
         means = list(c(-1, 0), c(1, 0)),
         covariances = list(diag(2), diag(2)))
A <- matrix(c(1, 0), nrow = 1L)
gmm_observe(g, A = A, y = 0.8, noise_cov = matrix(0.25, 1, 1))
#> <observe(gmm)>: K = 2 components in p = 2 dimensions
#>   [1] w = 0.2176, |mu| = 0.4400, tr(Sigma) = 1.2000
#>   [2] w = 0.7824, |mu| = 0.8400, tr(Sigma) = 1.2000
```
