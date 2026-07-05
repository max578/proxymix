# Classical EM fit on samples

Implements regime (ii) of Hoek and Elliott (2024). Runs the textbook
expectation-maximisation algorithm for Gaussian mixtures on the supplied
samples, with diagonal ridge regularisation for numerical stability,
optional multi-start, and monotone-log-likelihood checking.

## Usage

``` r
fit_em_samples(
  target,
  N = 2L,
  init = NULL,
  max_iter = 100L,
  tol = 1e-06,
  ridge_eps = 1e-06,
  n_starts = 5L,
  anneal = FALSE,
  temp_schedule = NULL,
  seed = NULL,
  canonicalise = TRUE
)
```

## Arguments

- target:

  A
  [gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
  carrying an `n` by `p` `samples` matrix.

- N:

  Number of mixture components.

- init:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md)
  initialisation, or `NULL` to use
  [`init_kmeans()`](https://max578.github.io/proxymix/reference/init_kmeans.md).

- max_iter:

  Maximum number of EM iterations.

- tol:

  Relative-log-likelihood convergence tolerance.

- ridge_eps:

  Ridge added to each component covariance at every M-step.

- n_starts:

  Number of multi-start initialisations (only when `init` is `NULL`).
  The best fit by final log-likelihood is returned.

- anneal:

  Logical. If `TRUE`, a deterministic-annealing warm-start (see
  [`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md))
  replaces the multi-start: the components are annealed from a high
  temperature down to one, and the resulting parameters seed a single
  final (cold) EM polish. This attacks the local-optima sensitivity of
  cold EM at the cost of the schedule length. Defaults to `FALSE` (cold
  best-of-`n_starts`).

- temp_schedule:

  Optional numeric vector of descending temperatures for the annealing
  warm-start. `NULL` (the default) uses a geometric schedule from `10`
  down to `1` in covariance-whitened units. Ignored when
  `anneal = FALSE`.

- seed:

  Optional integer seed for the annealing perturbations (the warm-start
  is deterministic given a seed). Ignored when `anneal = FALSE`.

- canonicalise:

  Logical. If `TRUE` (the default), the fitted mixture is post-processed
  by
  [`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md)
  so that components are sorted by descending weight and (as a
  tiebreaker) by descending `||mu||`.

## Value

A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) with
`regime = "sample"`. When `anneal = TRUE` the diagnostics list also
carries `annealed = TRUE` and the `temp_schedule` used.

## See also

Other fitting:
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md),
[`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md),
[`from_kde()`](https://max578.github.io/proxymix/reference/from_kde.md),
[`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md),
[`select_N()`](https://max578.github.io/proxymix/reference/select_N.md)

## Examples

``` r
x <- matrix(stats::rnorm(200), ncol = 2)
tgt <- gmm_target_from_samples(x)
fit <- fit_em_samples(tgt, N = 2L, max_iter = 30L, n_starts = 2L)
fit@diagnostics$loglik_final
#> [1] -284.1687
```
