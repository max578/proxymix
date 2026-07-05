# Modes of a Gaussian mixture

Returns the distinct local modes of a Gaussian-mixture density by
Gaussian mean-shift (the fixed-point hill-climb of Carreira-Perpinan
2000) started from each component mean, together with the mixture
density at each mode. Nearby converged points are merged so that each
genuine mode is reported once.

## Usage

``` r
gmm_modes(object, starts = NULL, tol = 1e-05, dedup = NULL, max_iter = 200L)
```

## Arguments

- object:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md).

- starts:

  Optional `m`-by-`p` matrix of starting points for the mean-shift.
  Default `NULL` uses the component means.

- tol:

  Convergence tolerance on the mean-shift step length. Default `1e-5`.

- dedup:

  Distance below which two converged points are treated as the same
  mode. Default `NULL` derives a small fraction of the spread of the
  component means.

- max_iter:

  Maximum mean-shift iterations per start. Default `200L`.

## Value

A list with `modes` (an `n`-by-`p` matrix of distinct modes, ordered by
descending mixture density), `density` (the mixture density at each
mode), and `n` (the number of modes).

## Details

A mixture of `K` components has at most `K` modes, and fewer when
components overlap; mean-shift from every component mean finds them
robustly without a grid. The companion of
[`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md):
applied to the fitted map it returns the recovered optima (ordered by
density, which for a Gibbs proxy ranks the deepest optima first).

## See also

[`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md),
whose fitted map this resolves into optima.

Other ops:
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md),
[`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md),
[`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md),
[`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md),
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md),
[`gmm_marginalise()`](https://max578.github.io/proxymix/reference/gmm_marginalise.md),
[`gmm_mean()`](https://max578.github.io/proxymix/reference/gmm_mean.md),
[`pgmm()`](https://max578.github.io/proxymix/reference/pgmm.md),
[`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md)

## Examples

``` r
g <- gmm(
  weights = rep(1 / 3, 3),
  means = list(c(-3, 0), c(3, 0), c(0, 4)),
  covariances = rep(list(0.3 * diag(2)), 3)
)
gmm_modes(g)$modes
#>      [,1]         [,2]
#> [1,]   -3 3.209642e-18
#> [2,]    3 3.209642e-18
#> [3,]    0 4.000000e+00
```
