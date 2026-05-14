# Affine-Gaussian operator calculus — v0.3 design note

> Pre-implementation design pass for proxymix v0.3. Fixes the maths, the
> closed-form updates, and the test obligations before any R code lands.
> Locked before v0.3.0 implementation.

## Why an operator calculus, not yet-more-fitters

Once a Gaussian-mixture proxy
\[
g(x) = \sum_{k=1}^{K} \pi_k \, \mathcal{N}(x; \mu_k, \Sigma_k)
\quad x \in \mathbb{R}^p
\]
has been fitted (regime (i), (ii) or (iii)), the user invariably wants
to do *more*: push through a sensor, condition on a noisy linear
observation, aggregate into a coarser variable, marginalise over
unobserved coordinates. All four are **closed-form** when the
transformation is affine-Gaussian — and that closed form is the
finite-mixture analogue of the Kalman update. v0.3 surfaces that
calculus as four exported operators.

Without v0.3 the user has only `gmm_marginalise()` (coordinate-subset)
and `gmm_conditionalise()` (Schur-complement, exact conditioning on
fully observed coordinates). Both are special cases of `gmm_observe()`
with `R = 0` and a particular `A`. v0.3 generalises the algebra.

## Mathematical claims

For a single Gaussian component `N(mu, Sigma)` in `R^p` and an
affine-Gaussian channel `y = A x + b + epsilon`, `epsilon ~ N(0, R)` in
`R^m`:

1. **Pushforward (predict step).** `y ~ N(A mu + b, A Sigma A' + R)`.
2. **Joint.** `(x, y) ~ N([mu; A mu + b], [[Sigma, Sigma A']; [A Sigma, A Sigma A' + R]])`.
3. **Conditional `x | y` (update step).** Apply standard Gaussian
   conditioning to the joint:
   \[
   \mathbb{E}[x | y] = \mu + \Sigma A^\top S^{-1} (y - A \mu - b),
   \qquad
   \mathrm{Cov}(x | y) = \Sigma - \Sigma A^\top S^{-1} A \Sigma
   \]
   where \(S = A \Sigma A^\top + R\) is the predictive covariance of
   `y` at the component.
4. **Marginal evidence.** \(p_k(y) = \mathcal{N}(y; A \mu_k + b, A \Sigma_k A^\top + R)\).

For a Gaussian mixture, these are applied **component-wise**, and the
component weights are re-normalised by the marginal evidence:

\[
\pi'_k = \frac{\pi_k\, p_k(y)}{\sum_j \pi_j\, p_j(y)}.
\]

That is the entire operator calculus. The exactness rests on (a)
linearity of the channel and (b) Gaussianity of `epsilon`. Both are
*non-negotiable* — see "Honesty constraints" below.

## Operator surface (v0.3.0)

### `gmm_affine(g, A, b = 0, noise_cov = NULL)`

Pushforward through `y = A x + b + epsilon`. Returns a `gmm` in `R^m`.

* `A`: `m`-by-`p` numeric matrix.
* `b`: numeric vector of length `m` (default zero).
* `noise_cov`: `m`-by-`m` SPD matrix, or `NULL` (deterministic channel,
  treated as `R = 0`).

For component `k`:
\[
\mu'_k = A \mu_k + b, \qquad
\Sigma'_k = A \Sigma_k A^\top + R_{\text{effective}},
\]
with weights unchanged. `R_effective = noise_cov` if supplied, else the
zero matrix.

**Edge case.** If `noise_cov` is `NULL` *and* `A` is rank-deficient, the
output covariance is rank-deficient. We accept that — the resulting
`gmm` is mathematically a degenerate Gaussian; subsequent operators
that need invertibility raise an explicit error.

### `gmm_observe(g, A, y, noise_cov)`

Condition `g` on a noisy linear observation `y = A x + b + epsilon`,
`epsilon ~ N(0, noise_cov)`. Returns a `gmm` in `R^p`.

For each component, apply the Kalman update to `(mu_k, Sigma_k)` and
multiply the weight by the marginal evidence `dnorm(y; A mu_k, S_k)`;
renormalise the weights.

**Edge case — vanishing weights.** If every per-component evidence is
numerically zero (e.g. the observation is many standard deviations from
*every* component), we issue `cli_warn()` and return the prior mixture
unchanged with a `metadata$gmm_observe_no_update = TRUE` flag rather
than producing `NaN` weights.

### `gmm_aggregate(g, A, noise_cov = NULL)`

A named alias of `gmm_affine()` for the case where `A` is a row-wise
aggregation matrix (e.g. block-sum, group-mean). Provides a clearer
public API for the downscaling use case and a hook for aggregation-
specific diagnostics in a later release.

### `gmm_missing(g, observed, values)`

Condition `g` on the *full* observation of a subset of coordinates.
Closed-form special case of `gmm_observe(g, A = selection_matrix,
y = values, noise_cov = 0)`. Re-uses `gmm_conditionalise()`'s Schur
machinery for efficiency.

* `observed`: integer indices in `seq_len(gmm_dim(g))`.
* `values`: numeric vector of length `length(observed)`.

## Honesty constraints (audit-mandated)

1. **Linearity is required.** Any non-linear `A`, `b`, or non-Gaussian
   `epsilon` is **not closed form** and must not be silently approximated.
   The pure operators error on non-numeric `A` / `b`. A separate Monte
   Carlo helper `gmm_pushforward_mc(g, channel_fn, n)` may be added
   later, *clearly labelled* and *not* part of the closed-form surface.
2. **No "approximate closed form".** The audit's stop criterion: if the
   maths is not exact, the function does not ship in this calculus.
3. **Mediator covariate shift is a research question.** Not encoded in
   operator signatures. Belongs in the v0.3 vignette and forthcoming
   manuscript draft.
4. **Comparison to GP / deconditional framework.** A vignette section
   contrasts the finite-mixture latent approach with the GP latent
   approach. GPs win on flexibility and uncertainty richness; mixtures
   win on closed-form aggregation, conditioning, and scenario
   generation.

## Numerical hygiene

* All M-step / update covariance matrices are symmetrised
  (`symmetrise(S)`) and ridged when supplied as inputs are nearly
  singular. We expose `ridge_eps` as a kwarg on every operator with the
  same default (1e-6) as the rest of the package.
* The Kalman gain is computed via `chol(S)` and `chol2inv()` rather
  than `solve(S)`. If `chol()` fails on a component, that component's
  update returns a near-zero gain and a warning is issued.
* Per-component log-evidence uses `mvnfast::dmvn(..., log = TRUE)`
  and a max-shifted `logsumexp` for renormalisation, so vanishing
  evidence never produces `0 / 0`.

## Public API stability

| Function | Signature (frozen at v0.3.0 entry) |
|---|---|
| `gmm_affine` | `gmm_affine(g, A, b = 0, noise_cov = NULL, ridge_eps = 1e-6)` |
| `gmm_observe` | `gmm_observe(g, A, y, noise_cov, b = 0, ridge_eps = 1e-6)` |
| `gmm_aggregate` | `gmm_aggregate(g, A, noise_cov = NULL, ridge_eps = 1e-6)` |
| `gmm_missing` | `gmm_missing(g, observed, values)` |

Return type: every operator returns a `gmm` (not a `gmm_fit`) so that the
operators compose freely. If `g` was a `gmm_fit`, the returned `gmm`
carries `metadata$source_fit = TRUE` and `metadata$source_fit_call`.

## Test obligations (v0.3.0 entry gate)

* **A0 — identity channel.** `gmm_affine(g, diag(p), b = 0)` is
  equal to `g` up to numerical tolerance.
* **A1 — affine of moments.** For each component, the resulting `mu'_k`
  matches `A %*% mu_k + b` and `Sigma'_k` matches `A Sigma A' + R`.
* **A2 — sum-of-Gaussians sanity.** For a single-component `gmm` with
  `A = I` and `R = R_0`, the result is the convolution
  `N(mu, Sigma + R_0)`.
* **O0 — Kalman parity.** For a single-component `gmm`, `gmm_observe()`
  matches a hand-coded Kalman update on `(mu, Sigma)`.
* **O1 — vanishing-evidence guard.** If `y` is 50 sigma away from every
  component, the function warns and returns the prior unchanged.
* **O2 — Bayes consistency.** Two sequential `gmm_observe()` calls on
  independent observations equal a single call on the stacked
  observation.
* **G0 — `gmm_aggregate()` is `gmm_affine()` for an aggregation matrix.**
  Round-trip parity test.
* **M0 — `gmm_missing()` matches `gmm_conditionalise()`** on a `gmm` with
  identity covariance.
* **C0 — Composition with `gmm_marginalise()`.** Pushforward, then
  marginalise, equals marginalise of the pushforward (i.e. operators
  commute where the maths says they should).

All tests live in `tests/testthat/test-operator-calculus.R`.

## Performance budget

For a `gmm` with `K = 10`, `p = 10`, `m = 10`:

* `gmm_affine`: O(K * (m * p + p * p * m)) per call; effectively
  K matrix multiplies of size up to (m, p) by (p, p). Target < 1 ms.
* `gmm_observe`: O(K * (m^3 + m * p)) per call; dominated by the
  per-component `chol(S)`. Target < 1 ms.

Both well within the regime-(iii) IS loop's per-iteration budget.

## What is explicitly *not* in v0.3.0

* `gmm_pushforward_mc()` for non-affine channels (out-of-scope; Tier 3).
* Mixture-of-experts conditional channels (i.e. `epsilon` whose
  distribution depends on `x`). Tier 3.
* Sequential filtering / smoothing wrappers on the operators. Tier 3.
* Numerical-only `gmm_apsim_aggregate()` (lives in a separate APSIM-bridge package).

## Pre-release gate

Kill criteria, before declaring v0.3.0 done:

1. **A0–C0 all green** on CI matrix.
2. **Vignette `operator_calculus.Rmd` builds clean**, shows the Kalman
   parity, the GP vs mixture contrast, and one downscaling example.
3. **NEWS entry articulates the *finite-mixture analogue of the Kalman
   update*** without overclaiming novelty.
4. **No silent non-affine fallbacks** — verify by code-review on the
   `R/operator_calculus.R` source.
5. **Companion manuscript skeleton** (`manuscripts/operator-calculus/`)
   created with section headings and method-statement boxes.

If any of (1)–(5) fail: hold the release, fix the failing gate, retest.

## References

* Hoek, J. van der and Elliott, R. J. (2024). *Mixtures of multivariate
  Gaussians.* Stochastic Analysis and Applications.
  <doi:10.1080/07362994.2024.2372605>.
* Kalman, R. E. (1960). *A new approach to linear filtering and
  prediction problems.* Transactions of the ASME — J. Basic Engineering.
* Murphy, K. P. (2012). *Machine Learning: A Probabilistic
  Perspective.* MIT Press. (Ch. 4: Gaussian models; closed-form
  conditioning.)
* Sejdinovic, D. et al. *Deconditional embedding for downscaling*
  (oral, Adelaide 2024-08; transcript at
  `references/Sejdinovic_talk1_unmatched_causal_and_downscaling_transcript.txt`).
