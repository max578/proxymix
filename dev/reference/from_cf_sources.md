# Fit a characteristic function that is a linear mix of independent sources

Fits a distribution of the form \\X = A S + G\\, where the coordinates
of \\S\\ are independent and not Gaussian, \\A\\ is an invertible matrix
and \\G\\ is Gaussian noise or absent. Only the characteristic function
of \\X\\ is needed. The mixing matrix is read from `cf`, and each source
is then fitted in one dimension with
[`from_cf()`](https://max578.github.io/proxymix/dev/reference/from_cf.md).
A joint fit in `n_dim` dimensions would be less accurate.

## Usage

``` r
from_cf_sources(
  cf,
  n_dim,
  N,
  n_freq = 8192L,
  bandwidths = c(0.3, 0.6, 1.2, 2.5, 5),
  max_iter = 300L,
  seed = NULL
)
```

## Arguments

- cf, n_dim:

  As in
  [`from_cf()`](https://max578.github.io/proxymix/dev/reference/from_cf.md).
  `n_dim` must be at least 2.

- N:

  Number of mixture components for each source.

- n_freq, bandwidths, max_iter:

  Passed to
  [`from_cf()`](https://max578.github.io/proxymix/dev/reference/from_cf.md)
  for each source.

- seed:

  Optional integer seed. When supplied, the random number generator
  state is restored on exit.

## Value

A list of class `proxymix_cf_sources` with

- `mixing`, the matrix \\\hat A\\. Its columns are the mixing
  directions, scaled so that each coordinate of \\Z = \hat A^{-1} X\\
  has unit variance. Column order and sign are not identified.

- `sources`, a list of `n_dim` one-dimensional
  [gmm_fit](https://max578.github.io/proxymix/dev/reference/gmm_fit.md)
  objects, the fits to the coordinates of \\Z\\.

- `noise_cov`, the off-diagonal covariance of \\Z\\, with the diagonal
  set to zero. It is zero when there is no noise. The noise variance on
  the diagonal stays inside `sources`.

- `cf(t)`, the characteristic function of the fit at the rows of `t`.

- `draw(n)`, an `n` by `n_dim` matrix of draws from the fit.

- `project(direction, max_components = 1e5)`, the
  [gmm](https://max578.github.io/proxymix/dev/reference/gmm.md) of
  \\a^\top X\\ under the fit for a direction \\a\\. Sources whose
  coefficient is below `1e-3` of the largest are replaced by a single
  normal with the same mean and variance. It gives an error if the
  projection needs more than `max_components` components.

- `diagnostics`, with `off_diagonal` (the largest off-diagonal entry
  divided by the largest diagonal entry of the Hessian difference at the
  third frequency, transformed by the inverse of `mixing`, 0 if the
  sources are separated exactly, with a warning above 0.05),
  `dependence` (the largest over pairs of coordinates of \\Z\\ of the
  root of the integral, against the standard normal density in two
  dimensions, of the squared modulus of the joint characteristic
  function minus the product of the marginal ones, with the Gaussian
  noise term allowed for; 0 for independent sources, with a warning
  above 0.01), `imaginary` (the largest imaginary part of the directions
  before it is dropped) and `heldout_rms` (the root mean square of the
  difference between `cf` and the fit's characteristic function on fresh
  frequencies).

## Details

The fitted family, independent sources that are each a one-dimensional
Gaussian mixture, is the independent factor model of Attias (1999), who
fits it to a sample by maximum likelihood.

The result is a list, not a
[gmm](https://max578.github.io/proxymix/dev/reference/gmm.md). As a
single mixture it would need `N^n_dim` components. It has three
functions: `cf()`, `draw()` (sampling from the fit) and `project()` (the
one-dimensional mixture along a direction). `draw()` gives an error when
the noise covariance estimated from `cf` is too far from positive
semi-definite to be sampled.

The mixing matrix is found from second derivatives of the logarithm of
`cf` at two random frequencies (Yeredor, 2000). The coordinates it
separates are then tested for pairwise independence with the criterion
of Eriksson and Koivunen (2003), which compares the joint characteristic
function of each pair with the product of the two marginal ones at all
frequencies. If the criterion is above 0.01 the frequencies are drawn
again, up to four times, and the directions with the smallest value are
kept.

The target must have a finite, positive definite covariance. A warning
of class `proxymix_cf_not_sources` is raised when the criterion is still
above 0.01, or when the directions do not separate the sources at a
third frequency. Either is the sign that the target is not of the form
above.

## References

Attias, H. (1999). Independent factor analysis. *Neural Computation*,
11(4), 803–851.
[doi:10.1162/089976699300016458](https://doi.org/10.1162/089976699300016458)

Eriksson, J., and Koivunen, V. (2003). Characteristic-function-based
independent component analysis. *Signal Processing*, 83(10), 2195–2208.
[doi:10.1016/S0165-1684(03)00162-2](https://doi.org/10.1016/S0165-1684%2803%2900162-2)

Yeredor, A. (2000). Blind source separation via the second
characteristic function. *Signal Processing*, 80(5), 897–902.
[doi:10.1016/S0165-1684(00)00062-1](https://doi.org/10.1016/S0165-1684%2800%2900062-1)

## See also

[`from_cf()`](https://max578.github.io/proxymix/dev/reference/from_cf.md),
[`from_cf_elliptical()`](https://max578.github.io/proxymix/dev/reference/from_cf_elliptical.md).

Other fitting:
[`fit_em_samples()`](https://max578.github.io/proxymix/dev/reference/fit_em_samples.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/dev/reference/fit_kld_em.md),
[`fit_moment_match()`](https://max578.github.io/proxymix/dev/reference/fit_moment_match.md),
[`fit_mppca_samples()`](https://max578.github.io/proxymix/dev/reference/fit_mppca_samples.md),
[`from_cf()`](https://max578.github.io/proxymix/dev/reference/from_cf.md),
[`from_cf_elliptical()`](https://max578.github.io/proxymix/dev/reference/from_cf_elliptical.md),
[`from_kde()`](https://max578.github.io/proxymix/dev/reference/from_kde.md),
[`from_objective()`](https://max578.github.io/proxymix/dev/reference/from_objective.md),
[`select_N()`](https://max578.github.io/proxymix/dev/reference/select_N.md)

## Examples

``` r
## Two independent Gamma(3) sources, mixed, with N(0, 0.25 I) noise.
A <- matrix(c(0.7, 0.3, 0, -0.7), 2L)
cf <- function(t) {
  u <- t %*% A
  (1 - 1i * u[, 1])^(-3) * (1 - 1i * u[, 2])^(-3) *
    exp(-0.125 * rowSums(t^2))
}
fit <- from_cf_sources(cf, n_dim = 2L, N = 3L, n_freq = 1024L, seed = 1L)
fit$diagnostics$off_diagonal
#> [1] 2.906809e-06
first <- fit$project(c(1, 0))
gmm_cf_distance(first, function(s) cf(cbind(s[, 1], 0)))$distance
#> [1] 0.002356499
```
