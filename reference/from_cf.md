# Fit a Gaussian mixture to a characteristic function

Fits a Gaussian mixture to a distribution that is known through its
characteristic function \\\varphi(t) = E \exp(i t^\top X)\\. No density
and no sample of the distribution are needed. The result is an ordinary
fitted mixture, so
[`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md),
[`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md),
[`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md),
[`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md)
and the other mixture operations apply to it.

## Usage

``` r
from_cf(
  cf,
  n_dim,
  N,
  n_freq = 8192L,
  bandwidths = c(0.3, 0.6, 1.2, 2.5, 5),
  n_start = 1L,
  max_iter = 300L,
  init = NULL,
  seed = NULL,
  bracket = TRUE
)
```

## Arguments

- cf:

  A function that takes an `n` by `n_dim` matrix of frequencies, one per
  row, and returns the characteristic function at each row as a complex
  vector of length `n`. A real vector is accepted for a symmetric
  target.

- n_dim:

  Dimension of the target.

- N:

  Number of mixture components.

- n_freq:

  Number of frequencies in the fit, at least 64 times
  `length(bandwidths)`. Default `8192L`.

- bandwidths:

  Increasing positive numbers, the typical length of a frequency vector
  in units of the target's half-height frequency (each coordinate has
  standard deviation `bandwidth / sqrt(n_dim)` of that). The grown start
  and `init` use all bandwidths together. Random starts are fitted from
  the lowest bandwidth upwards. Extend the range upwards for a target
  with fine structure and downwards for a heavy-tailed target.

- n_start:

  Number of starts. The first start is `init` if given, otherwise the
  grown start described above; the others are random. The fit with the
  lowest loss is returned. Default `1L`. With ten or more components in
  five or more dimensions a single start can stop at a poor fit, so four
  starts are a safer choice. For ten components in ten dimensions the
  grown start stopped at a poor fit and a random start recovered the
  target.

- max_iter:

  Maximum number of L-BFGS iterations per frequency stage. A
  Levenberg-Marquardt step with at most 200 iterations follows when the
  fit has at most 700 parameters. Default `300L`.

- init:

  Optional [gmm](https://max578.github.io/proxymix/reference/gmm.md)
  with `N` components to start the first fit from, in place of the grown
  start. Default `NULL`.

- seed:

  Optional integer seed. When supplied, the random number generator
  state is restored on exit.

- bracket:

  Logical. Compute the Kolmogorov bounds described under *Checking a
  fit*. Default `TRUE`.

## Value

A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) with
`regime` set to `"cf"`. Its `diagnostics` hold `loss` (the final loss),
`loss_by_start`, `loss_path` (the loss as the grown start went from 1 to
`N` components, or `NULL` when `init` is given), `heldout_rms` and
`heldout_max`, `relocations` (the number of moves kept), `tail_index`
(one value for each coordinate), `cf_mean` and `cf_cov` (the target's
mean and covariance read from `cf`, or `NULL` when the target has no
finite, positive definite covariance), `n_freq`, `bandwidths`, and, with
`bracket = TRUE`, `bracket` and `bracket_directions`. `converged` is
`FALSE` when the final optimisation of the returned fit stopped at its
iteration limit.

## Details

The fit minimises the mean of \\\|\varphi(t_j) - \varphi_g(t_j)\|^2\\
over a fixed sample of frequencies \\t_j\\, where \\\varphi_g\\ is the
characteristic function of the mixture
([`gmm_cf()`](https://max578.github.io/proxymix/reference/gmm_cf.md)).
The frequencies are Gaussian with several bandwidths, scaled to the
target so that `bandwidths = 1` is the frequency at which
\\\|\varphi\|\\ has fallen to about one half. The same criterion is used
by Brück (2025) to train a neural sampler.

The starting mixture is built one component at a time. Each new
component is placed at the mean where the fit falls furthest short of
the target, after both are smoothed by a Gaussian. That mean is computed
from the characteristic-function residual, so no density of the target
is needed. All components are refitted after each addition. The fit then
tries moving each of the two least useful components to a better mean,
keeping a move if it cuts the loss by at least 10%, up to `N` moves.

## Checking a fit

The returned object reports the root mean square and the maximum of
\\\|\varphi - \varphi_g\|\\ on `n_freq` fresh frequencies that the fit
did not use (`heldout_rms` and `heldout_max` in `diagnostics`).

With `bracket = TRUE` (the default) the fit also carries a lower and an
upper bound for the Kolmogorov distance between the target and the fit,
computed by
[`gmm_cf_distance()`](https://max578.github.io/proxymix/reference/gmm_cf_distance.md)
along eight directions: the first `min(n_dim, 8)` coordinate axes and,
when `n_dim < 8`, random unit directions to make up eight (one direction
when `n_dim = 1`). The bounds for each direction are in
`diagnostics$bracket`, the directions in
`diagnostics$bracket_directions`, and the largest of each over the
directions in the quality record as `cf_distance` and `cf_upper`
([`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md)).
These bounds compare distribution functions; they are not a
Kullback-Leibler divergence. Drawing the random directions does not
change the state of the random number generator; with `seed` they follow
the seed. Set `bracket = FALSE` to skip the bounds.

## Heavy tails

A Gaussian mixture has finite moments of every order, so it cannot
reproduce infinite-variance tails and fits the centre only. The tail
index is read from \\-\log\|\varphi\|\\ near zero. It is 2 when the
variance is finite and equals the stability index for a stable law. When
it is below 1.9 a warning of class `proxymix_cf_heavy_tail` is raised.

## Scope

The target needs a characteristic function that decays, which excludes
lattice distributions and distributions with a point mass.

## References

Brück, F. (2025). Generative neural networks for characteristic
functions. *Journal of Computational and Graphical Statistics*, 34(4),
1456–1465.
[doi:10.1080/10618600.2025.2455135](https://doi.org/10.1080/10618600.2025.2455135)

## See also

[`cf_target()`](https://max578.github.io/proxymix/reference/cf_target.md)
and
[`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
for the same fit through the unified verb;
[`from_cf_elliptical()`](https://max578.github.io/proxymix/reference/from_cf_elliptical.md)
and
[`from_cf_sources()`](https://max578.github.io/proxymix/reference/from_cf_sources.md)
for structured targets;
[`gmm_cf()`](https://max578.github.io/proxymix/reference/gmm_cf.md),
[`gmm_cf_distance()`](https://max578.github.io/proxymix/reference/gmm_cf_distance.md).

Other fitting:
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md),
[`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md),
[`fit_mppca_samples()`](https://max578.github.io/proxymix/reference/fit_mppca_samples.md),
[`from_cf_elliptical()`](https://max578.github.io/proxymix/reference/from_cf_elliptical.md),
[`from_cf_sources()`](https://max578.github.io/proxymix/reference/from_cf_sources.md),
[`from_kde()`](https://max578.github.io/proxymix/reference/from_kde.md),
[`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md),
[`select_N()`](https://max578.github.io/proxymix/reference/select_N.md)

## Examples

``` r
## The sum of a Gamma(3) variable and independent N(0, 0.25) noise.
cf <- function(t) (1 - 1i * t[, 1])^(-3) * exp(-0.125 * t[, 1]^2)
fit <- from_cf(cf, n_dim = 1L, N = 3L, n_freq = 1024L, seed = 1L)
## the print shows the held-out error and the Kolmogorov bounds
fit
#> <gmm_fit>: regime = "cf", K = 3, p = 1
#>   target     : <none>
#>   iterations : NA
#>   converged  : TRUE
#>   cf held-out RMS error : 0.00477
#>   Kolmogorov distance   : at least 0.00398, at most 0.00859 (1 direction)
#>   [1] w = 0.4671, |mu| = 2.9662, tr(Sigma) = 1.4310
#>   [2] w = 0.3043, |mu| = 1.5676, tr(Sigma) = 0.6571
#>   [3] w = 0.2286, |mu| = 4.8617, tr(Sigma) = 3.6102
fit@diagnostics$bracket
#>   direction    distance      upper upper_inside   truncation   cutoff_gap
#> 1         1 0.003982315 0.00859195  0.006278226 8.263034e-14 7.939929e-12
```
