# Fit concentric normals to an elliptical characteristic function

Fits a Gaussian mixture to a symmetric distribution whose characteristic
function has the form \\\varphi(t) = \exp(i t^\top \mu)\\ \rho(t^\top
\Sigma t)\\, where \\\rho(u)\\ is a mixture of \\\exp(-s u / 2)\\ over
\\s \> 0\\ (a Gaussian scale mixture). Elliptical stable laws, the
multivariate Cauchy, Student and Laplace distributions are of this form.
Every component has mean \\\mu\\ and a covariance proportional to
\\\Sigma\\, and the fit reaches further into heavy tails than
[`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md)
does.

## Usage

``` r
from_cf_elliptical(
  cf,
  n_dim,
  N = 20L,
  scale_range = NULL,
  seed = NULL,
  bracket = TRUE
)
```

## Arguments

- cf, n_dim:

  As in
  [`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md).

- N:

  Number of values on the grid of component standard deviations. Default
  `20L`.

- scale_range:

  Smallest and largest multiplier \\s\\. A component has covariance
  \\s^2 \Sigma\\, and \\s = 1\\ is the normal whose half-height
  frequency matches the target's. The default `NULL` reads the tail
  index from `cf` and runs from 0.05 to 30 times the 0.999 quantile of a
  stable law with that index (to 20 when the index is 1.9 or more), both
  in units of that law's scale.

- seed:

  Optional integer seed for the random directions. When supplied, the
  random number generator state is restored on exit.

- bracket:

  Logical. Compute Kolmogorov bounds against `cf` along eight
  directions, as in
  [`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md).
  Default `TRUE`.

## Value

A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) with
`regime` set to `"cf"`. Its `diagnostics` hold `shape` (the matrix
\\\Sigma\\), `shape_ratio`, `asymmetry` (the largest imaginary part of
the centred characteristic function, relative to its modulus), `scales`
(all `N` grid values, including those dropped), `tail_index`,
`radial_error` (the largest absolute difference between \\\rho\\ and its
fit on the frequencies used), `heldout_rms` and `heldout_max`, and
`bracket` and `bracket_directions`, as in
[`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md).

## Details

The centre \\\mu\\ is read from the phase of `cf`. The matrix \\\Sigma\\
is found by least squares from the frequency at which \\\|\varphi\|\\
falls to one half along each of a set of random directions. The standard
deviations of the components are `N` values on a geometric grid. Their
weights are found by non-negative least squares on the relative error of
\\1 - \rho\\, over frequencies from very small to well past the
half-height. The smallest frequencies fix the tails. Grid values that
receive no weight are dropped, so the result can have fewer than `N`
components.

`diagnostics$shape_ratio` should be 1 for an elliptical target. It is
the smallest and largest value, over 200 random directions, of fitted
scale divided by target scale. A warning of class
`proxymix_cf_not_elliptical` is given if either value is outside 0.9 to
1.1, or if the imaginary part of the centred characteristic function
exceeds 0.001 of its modulus.

## References

Lawson, C. L. and Hanson, R. J. (1974). *Solving Least Squares
Problems*. Prentice-Hall.

## See also

[`cf_target()`](https://max578.github.io/proxymix/reference/cf_target.md),
[`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md),
[`from_cf_sources()`](https://max578.github.io/proxymix/reference/from_cf_sources.md),
[`gmm_cf_distance()`](https://max578.github.io/proxymix/reference/gmm_cf_distance.md).

Other fitting:
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md),
[`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md),
[`fit_mppca_samples()`](https://max578.github.io/proxymix/reference/fit_mppca_samples.md),
[`from_cf()`](https://max578.github.io/proxymix/reference/from_cf.md),
[`from_cf_sources()`](https://max578.github.io/proxymix/reference/from_cf_sources.md),
[`from_kde()`](https://max578.github.io/proxymix/reference/from_kde.md),
[`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md),
[`select_N()`](https://max578.github.io/proxymix/reference/select_N.md)

## Examples

``` r
## A bivariate Cauchy distribution.
S <- matrix(c(2, 0.5, 0.5, 1), 2L)
cf <- function(t) exp(1i * t[, 1] - sqrt(rowSums((t %*% S) * t)))
fit <- from_cf_elliptical(cf, n_dim = 2L, seed = 1L)
## the print shows the Kolmogorov bounds along eight directions
fit
#> <gmm_fit>: regime = "cf", K = 14, p = 2
#>   target     : <none>
#>   iterations : NA
#>   converged  : TRUE
#>   cf held-out RMS error : 0.000206
#>   Kolmogorov distance   : at least 9.83e-05, at most 0.005 (8 directions)
#>   [1] w = 0.3014, |mu| = 1.0000, tr(Sigma) = 4.5137
#>   [2] w = 0.2342, |mu| = 1.0000, tr(Sigma) = 1.2550
#>   [3] w = 0.1986, |mu| = 1.0000, tr(Sigma) = 16.2342
#>   [4] w = 0.1137, |mu| = 1.0000, tr(Sigma) = 58.3883
#>   [5] w = 0.0599, |mu| = 1.0000, tr(Sigma) = 210.0006
#>   ... 9 more components
fit@diagnostics$shape_ratio
#> [1] 0.9999999 1.0000002
## the 0.99 quantile of the first coordinate against the exact value
qgmm(0.99, gmm_marginalise(fit, 1L))
#> [1] 45.98639
1 + sqrt(2) * qcauchy(0.99)
#> [1] 46.00101
```
