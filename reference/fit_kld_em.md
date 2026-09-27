# Importance-sampled KLD-EM fit (regime iii)

Implements regime (iii) of Hoek and Elliott (2024). Minimises
`KL(f || g_theta)` where `f` is supplied as an evaluable log-density on
the target, via expectation-maximisation against importance-sampled
draws from a user-chosen proposal `q`.

## Usage

``` r
fit_kld_em(
  target,
  N = 3L,
  proposal = NULL,
  is_size = 5000L,
  init = NULL,
  max_iter = 100L,
  tol = 1e-05,
  ridge_eps = 1e-06,
  min_ess = 50,
  on_low_ess = c("warn", "abort"),
  seed = NULL,
  validation_size = NULL,
  validation_proposal = NULL,
  validation_seed = NULL,
  support_warn = TRUE,
  adapt = c("none", "pmc"),
  refresh_every = 5L,
  defensive_gamma = 0.15,
  inflate = 1.5,
  anneal = FALSE,
  temp_schedule = NULL,
  canonicalise = TRUE
)
```

## Arguments

- target:

  A
  [gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md)
  with a non-NULL `log_density`.

- N:

  Number of mixture components.

- proposal:

  An
  [is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md).
  When `NULL` (the default) the proposal is chosen automatically: a
  support-matched
  [`is_uniform()`](https://max578.github.io/proxymix/reference/is_uniform.md)
  when the target declares a bounded or one-sided `support`, otherwise a
  multivariate-t with `df = 5` in `target@n_dim` dimensions. The
  automatic choice is announced with a one-line message so it is never
  silent.

- is_size:

  Number of importance-sampling draws used for fitting.

- init:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md)
  initialisation, or `NULL` to use a kmeans pass on the
  importance-resampled draws.

- max_iter:

  Maximum number of EM iterations.

- tol:

  Convergence tolerance on the relative change in the
  importance-weighted EM objective `Q(theta) = sum_n W_n log g(x_n)`.
  `Q` is invariant to the target's normalising constant, so the stopping
  rule behaves identically for normalised and unnormalised targets (the
  importance-sampled KLD estimate carries an additive `log Z(f)` offset
  and is therefore never used for stopping).

- ridge_eps:

  Ridge added to each component covariance at every M-step.

- min_ess:

  Minimum effective sample size below which the fit is flagged as
  degenerate: a classed warning (`proxymix_low_ess`) is issued (or, with
  `on_low_ess = "abort"`, a classed error `proxymix_degenerate_fit`),
  the fit's `converged` flag is forced to `FALSE`, and
  `degenerate = TRUE` is recorded in the diagnostics and the quality
  certificate.

- on_low_ess:

  What to do when the effective sample size falls below `min_ess`:
  `"warn"` (the default) flags and continues, `"abort"` refuses to
  return a degenerate fit.

- seed:

  Optional integer seed. When supplied, the fit is reproducible
  end-to-end: the fitting IS draw, the initialisation resample and
  kmeans pass, and any empty-component reseed draws are all derived from
  it. When `NULL`, those draws consume the ambient random-number stream.

- validation_size:

  Number of independent importance-sampling draws for the held-out
  estimate `validation_kld` and its standard error `validation_mc_se`.
  The standard error is roughly `sqrt(2 * KL / ESS)`, where ESS is the
  effective sample size of the validation draws, so a close fit needs
  many draws before the estimate is precise. The default `NULL` draws
  batches of `is_size` until the standard error is at most a tenth of
  the estimate (at most 0.001 when the estimate is below 0.01), stopping
  at `10 * is_size` draws. When the target's normalising constant is
  unknown, the goal is a standard error of 0.001. A number draws exactly
  that many; `0L` skips validation.

- validation_proposal:

  Optional
  [is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md)
  for the validation sample. Defaults to the same proposal used for
  fitting.

- validation_seed:

  Optional integer seed used when drawing the validation sample.
  Defaults to `seed + 1L` when `seed` is supplied, `NULL` otherwise.

- support_warn:

  Logical. If `TRUE` (the default), issue a warning when more than 5% of
  IS draws receive non-finite weights (typically because the proposal
  does not dominate the target's support).

- adapt:

  Proposal adaptation: `"none"` (the default; one fixed IS draw, the
  historical behaviour) or `"pmc"` (population-Monte-Carlo refresh of
  the proposal from the current iterate; see Details).

- refresh_every:

  With `adapt = "pmc"`, refresh the proposal after this many EM
  iterations on the current batch. Default `5L`.

- defensive_gamma:

  With `adapt = "pmc"`, the mass kept on the original proposal as a
  heavy-tailed defensive anchor at every refresh (bounds the
  importance-weight variance). Default `0.15`.

- inflate:

  With `adapt = "pmc"`, the factor inflating the current iterate's
  covariances inside the refreshed proposal. Default `1.5`.

- anneal:

  Logical. If `TRUE`, a deterministic-annealing warm-start (see
  [`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md))
  replaces the kmeans initialisation: components are annealed from a
  high temperature down to one on the importance-weighted draws, and the
  resulting parameters seed the (unchanged) cold KLD-EM loop. This
  attacks the local-optima sensitivity of cold EM. Defaults to `FALSE`.

- temp_schedule:

  Optional numeric vector of descending temperatures for the annealing
  warm-start. `NULL` (the default) uses a geometric schedule from `10`
  down to `1` in covariance-whitened units. Ignored when
  `anneal = FALSE`.

- canonicalise:

  Logical. If `TRUE` (the default), the fitted mixture is post-processed
  by
  [`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md).

## Value

A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md) with
`regime = "kld"`. The diagnostics list contains, among others,
`kld_trace`, `kld_final` (in-sample), `kld_is_shifted`,
`kld_final_absolute` (when computable), `ess`, `ess_relative`
(`ess / is_size`), `max_weight`, `support_fraction`, `mc_se_kld`,
`validation_kld` (held-out), `validation_ess`, `validation_max_weight`,
and `heldout_kld`, the held-out estimate of the KL divergence that
[`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md)
describes.

## Details

With `adapt = "none"` (the default) the Monte Carlo draws from `q` are
computed once at the start and the resulting self-normalised
importance-sampling weights are reused at every EM iteration. With
`adapt = "pmc"` the proposal is refreshed every `refresh_every`
iterations with a defensive mixture built from the current iterate – the
population-Monte-Carlo scheme: the fitted mixture (covariances inflated
by `inflate`) carries `1 - defensive_gamma` of the proposal mass and the
original proposal `q` keeps `defensive_gamma` as a heavy-tailed anchor,
a fresh IS batch is drawn, and EM continues on the refreshed weights.
Because the refreshed proposal tracks the target, the effective sample
size recovers from a poor initial proposal and the usable dimension
range extends well beyond what a fixed proposal reaches; the per-batch
ESS trace is reported as `diagnostics$ess_history`. While a batch is
degenerate (its effective sample size is below `min_ess`), the refresh
fires every iteration with an escalating covariance inflation floored at
a growing fraction of the batch's sample covariance, so a collapsed
iterate walks back out toward the target instead of freezing; and
convergence is only accepted on an adapted batch, so a run that
stabilises on the original proposal's draw is refreshed at least once
before it is allowed to stop. The scheme is the mixture
population-Monte-Carlo idea of Cappé et al. (2008) with the
defensive-mixture safeguard of Owen and Zhou (2000); it re-draws rather
than recycles batches (compare the adaptive multiple importance sampling
of Cornuet et al., 2012).

`kld_final` is estimated on the same importance draws the fit was tuned
to, so it reads low and can come out negative, more so for large `N` and
small `is_size`. When `validation_size > 0` the function also draws an
*independent* validation IS sample and reports its KLD estimate
(`validation_kld`), effective sample size, and largest weight share;
`validation_kld` is the estimate to quote as the fit's accuracy.

When the target's `normalised` property is `FALSE` or `NA`, the
importance-sampled `kld_final`, `kld_trace` and `validation_kld`
estimate \\KL(f \Vert g) + \log Z(f)\\, where \\Z(f)\\ is the integral
of `exp(log_density)`. The fit's diagnostics list records this via
`kld_is_shifted = TRUE` and a `kld_shift_explanation` string. When the
target also supplies a finite `log_normalizer`, the fit subtracts it and
reports the KL divergence as `kld_final_absolute` and
`validation_kld_absolute`.

## References

Cappé, O., Douc, R., Guillin, A., Marin, J.-M. and Robert, C. P. (2008)
Adaptive importance sampling in general mixture classes. *Statistics and
Computing* 18, 447–459.
[doi:10.1007/s11222-008-9059-x](https://doi.org/10.1007/s11222-008-9059-x)

Cornuet, J.-M., Marin, J.-M., Mira, A. and Robert, C. P. (2012) Adaptive
multiple importance sampling. *Scandinavian Journal of Statistics* 39,
798–812.
[doi:10.1111/j.1467-9469.2011.00756.x](https://doi.org/10.1111/j.1467-9469.2011.00756.x)

Owen, A. and Zhou, Y. (2000) Safe and effective importance sampling.
*Journal of the American Statistical Association* 95(449), 135–143.
[doi:10.1080/01621459.2000.10473909](https://doi.org/10.1080/01621459.2000.10473909)

## See also

Other fitting:
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md),
[`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md),
[`from_kde()`](https://max578.github.io/proxymix/reference/from_kde.md),
[`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md),
[`select_N()`](https://max578.github.io/proxymix/reference/select_N.md)

## Examples

``` r
tgt <- banana_target()
q <- is_mvt(n_dim = 2L, mean = c(0, 0),
            sigma = 4 * diag(2), df = 5)
fit <- fit_kld_em(tgt, N = 3L, proposal = q,
                  is_size = 1500L, max_iter = 25L, seed = 1L,
                  validation_size = 1500L)
fit@diagnostics$kld_final
#> [1] 0.006967304
fit@diagnostics$validation_kld
#> [1] 0.00926966
```
