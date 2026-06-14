# proxymix 0.4.0

### New features

* **New: the C4 consumer seam -- compress a flexyBayes posterior directly.**
  `from_fb_posterior()` is the consumer entry point for the constellation's
  C4 contract (`fb_log_posterior` / `posterior_proxy`, owned by flexyBayes):
  it takes a flexyBayes posterior addressed only through its (unnormalised)
  log-density and returns a closed-form Gaussian-mixture proxy via
  importance-sampled KLD-EM, generalising the input source from a kernel-
  density estimate (`from_kde()`) or a PESTO ensemble to any Bayesian
  posterior. The producer interface the seam expects is materialised and
  documented by `fb_log_posterior_spec()`; `fb_producer_available()` is a
  capability probe (degrades to `FALSE`, never errors, when no real producer
  is installed); and `mock_fb_posterior()` is a synthetic producer (known
  Gaussian or banana log-density) for testing the path with no sibling
  package present. **Activation:** the fitted-flexyBayes-object path switches
  on when flexyBayes lands its `fb_log_posterior` producer (the B4 follow-up).
  Until then, a fitted flexyBayes object raises an informative seam error,
  and the bare-callable / `mock_fb_posterior()` paths work today.
  proxymix never `Imports:` flexyBayes -- the seam is a soft contract and
  `R CMD check` is clean with no flexyBayes installed; the C4 return path
  (proxymix serving as a `posterior_proxy(type = "gmm")` backend) stays an
  adapter on the flexyBayes side, so the constellation's acyclic invariant
  holds.

* **New: `autoplot()` method for `gmm_fit`.** Render a fitted proxy with
  `ggplot2::autoplot(fit)` — a marginal density curve in one dimension, or a
  viridis density raster with per-component ellipses in two. Any ambient
  dimension is supported: the requested coordinates are reduced through the
  package's own closed-form `gmm_marginalise()` before plotting (e.g.
  `autoplot(fit, dims = c(1L, 3L))`). `ggplot2` stays an optional dependency
  — the method registers only when `ggplot2` is installed, so `R CMD check`
  remains clean with no sibling package present.

## proxymix 0.3.0 (2026-05-14)

Second methodological extension. Brings a complete affine-Gaussian
operator calculus to Gaussian-mixture proxies — pushforward,
Bayesian update on a noisy linear observation, aggregation,
missing-data conditioning — each closed-form and component-wise.

### User-visible changes

* **New: `gmm_affine(g, A, b, noise_cov)`** — closed-form pushforward
  of a Gaussian mixture through `y = A x + b + epsilon`,
  `epsilon ~ N(0, noise_cov)`. Returns the mixture in `R^m` with
  `mu'_k = A mu_k + b`, `Sigma'_k = A Sigma_k A' + noise_cov`,
  weights unchanged.
* **New: `gmm_observe(g, A, y, noise_cov)`** — Bayesian update on a
  noisy linear observation. Applies the Kalman gain per component and
  reweights component weights by per-component marginal evidence. The
  finite-mixture analogue of a Kalman update.
* **New: `gmm_aggregate(g, A, noise_cov)`** — named alias for
  `gmm_affine()` aimed at downscaling / aggregation pipelines.
* **New: `gmm_missing(g, observed, values)`** — Schur-complement
  conditioning routed through an integer-index API for missing-data
  pipelines.

### Design and validation

* `docs/design/operator_calculus_v0.3.md` — pre-implementation design
  note: maths, honesty constraints (no non-affine fallbacks, no
  approximate closed form), public API freeze, test obligations,
  performance budget, pre-release gate.
* `vignettes/operator_calculus.Rmd` — educational vignette with
  Kalman parity check, sequential vs stacked observations,
  aggregation through a coarsening matrix, and a comparison to a
  Gaussian-process latent.
* `inst/validation/operator_calculus_pinned.R` — three pinned
  reference pipelines (Kalman parity, sequential vs stacked,
  aggregate-then-observe) with hand-coded acceptance ranges.

### Tests

* `test-operator-calculus.R` (12 tests, 46 expectations): A0–A2
  (affine of moments), O0–O2 (Kalman parity, vanishing-evidence
  guard, Bayes consistency), G0 (aggregate alias), M0 (missing vs
  conditionalise), C0 (composition with marginalise), plus full
  input-validation coverage.

### Internal

* `R/operator_calculus.R` consolidates the four operators with
  shared validation helpers (`.validate_A`, `.validate_b`,
  `.validate_noise_cov`) and a single numerical-hygiene policy
  (ridge after each output covariance, symmetrisation, chol-based
  inverse with retry on near-singular matrices).
* `gmm_observe()` issues a `proxymix_observe_no_update` warning when
  the marginal evidence is numerically zero at every component and
  returns the prior unchanged with a metadata flag.

## proxymix 0.2.0 (2026-05-14)

First Tier-2 graduation. Two methodological extensions that compose
cleanly with the regime-(iii) wedge:

### User-visible changes

* **New: `from_kde()` (Tier-2 graduation).** Compiles a kernel density
  estimate over an `n` by `p` sample matrix into a closed-form
  Gaussian-mixture proxy via regime-(iii) KLD-EM. Supports scalar and
  diagonal bandwidths (`"silverman"`, `"scott"`, numeric scalar, or
  per-coordinate numeric vector). Dimensional guard: `p <= 5`
  recommended, `p <= 10` allowed with warning, `p > 10` rejected. The
  KDE-target is `normalised = TRUE` by construction, so downstream
  KLD and Hellinger diagnostics report absolute values. Companion
  vignette: `vignettes/from_kde.Rmd`.
* **New: `gmm_target_from_posterior()` (Contract A constructor).**
  S3 generic that compiles an (unnormalised) Bayesian posterior into a
  `gmm_target`. The `function` method accepts a bare vectorised callable
  with required `parameter_names`; the `default` method points users at
  either a registered Bayesian-package method (`flexyBayes`, `brms`,
  `Stan`, ...) or the function-based path. Vectorisation contract is
  enforced at construction by a probe call.
* **URL and BugReports.** `DESCRIPTION` now ships the canonical GitHub
  namespace at `github.com/max578/proxymix`.

### Tests

* `test-from-kde.R` (8 tests, 25 expectations): end-to-end recovery,
  bandwidth selection branches, dimensional guards, normalisation
  declaration, default proposal sanity, metadata pass-through.
* `test-from-posterior.R` (7 tests, 21 expectations): vectorisation
  contract enforcement, log-normalizer pass-through, default-method
  hinting, name validation, round-trip through `fit_proxymix(regime =
  "kld")`, attribute-based parameter-name support.
* `inst/validation/from_kde_pinned_fits.R`: pinned validation across
  three reference KDE -> GMM pipelines (bimodal, banana, mixture)
  with MC-SE-aware acceptance ranges.

### Documentation

* `vignettes/from_kde.Rmd`: educational walk-through covering scope,
  bandwidth sensitivity, recovery on a known mixture, and the
  contrast between KDE and proxy log-densities.

### Internal

* `gmm_target_from_posterior` registers an S3 generic, paving the way
  for `flexyBayes::gmm_target_from_posterior.flexybayes` (and analogous
  methods for `brms`, `Stan`, `pymc`-via-reticulate) without coupling
  proxymix to any specific Bayesian backend.
* From-KDE log-density evaluation uses chunked matrix builds so that
  peak memory stays bounded for large IS samples.

## proxymix 0.1.1 (2026-05-14)

Audit-driven scientific hardening pass. No new Tier-2 bodies; the wedge
is made harder to misuse.

### User-visible changes

* **Normalisation-aware targets.** `gmm_target` gains two new properties,
  `normalised` (logical or `NA`) and `log_normalizer` (numeric or `NA`),
  so that an unnormalised `log_density` can be supplied without making
  downstream KLD or Hellinger diagnostics misleading. All three built-in
  targets (`banana_target()`, `donut_target()`, `mixture_target()`)
  declare `normalised = TRUE`; the unnormalised case is now explicitly
  documented at the target level.
* **Canonical component ordering.** A new `gmm_canonicalise()` function
  reorders the components of a `gmm` (or `gmm_fit`) by descending
  weight, then by descending `||mu||` as a tiebreaker. `fit_proxymix()`
  and the regime-specific fitters now canonicalise their outputs by
  default (`canonicalise = TRUE`), making prints, snapshot tests, and
  cross-run comparisons reproducible. Set `canonicalise = FALSE` to
  retain the raw EM-order parameters.
* **Held-out importance-sample validation.** `fit_kld_em()` (and
  therefore `fit_proxymix(regime = "kld")`) accepts `validation_size`
  and `validation_proposal`. When `validation_size > 0`, a second
  independent IS sample is drawn and the fit's diagnostics list
  records `validation_kld`, `validation_ess`, and
  `validation_max_weight`. This lets users tell the difference between
  in-sample overfit and a fit that generalises across IS draws.
* **Richer IS diagnostics.** `fit_kld_em()` now records
  `ess_relative` (ESS / `is_size`), `max_weight` (largest
  self-normalised weight), `support_fraction` (fraction of IS draws
  with finite log-density under target *and* proposal), and a
  Monte-Carlo standard error for the final KLD estimate
  (`mc_se_kld`). A new `ess_summary()` helper returns the headline
  numbers as a small list.
* **Shifted-KLD labelling.** Diagnostics now record
  `kld_is_shifted` and `kld_shift_explanation` whenever the target is
  unnormalised or its normalisation is unknown, so users do not silently
  read a shifted MC integral as an absolute divergence.
* **Hellinger guard.** `hellinger_mc()` now warns when the target is
  not declared `normalised = TRUE` — the squared Hellinger distance is
  not meaningful against an unnormalised target.
* **Proposal-support warning.** `fit_kld_em()` issues a `cli` warning
  when more than 5% of importance-sample draws fall outside the
  proposal's support or carry non-finite weights. The most common
  trigger is an `is_uniform()` proposal whose box does not cover the
  target's mass.

### Validation corpus

* `inst/validation/regime_iii_pinned_fits.R` — a runnable validation
  script that fits the three built-in targets with pinned seeds and
  records final KLD, ESS, max weight, validation KLD, and runtime;
  intended as the seed of a growing `inst/validation/` corpus per the
  audit's recommendation.

### Tests

* New: `test-canonicalise.R`, `test-normalisation.R`,
  `test-validation-split.R`, `test-support-warning.R`, and
  `test-monotone-objective.R`. The last asserts monotonicity of the
  fixed IS-weighted objective \(\sum_n W_n \log g_\theta(x_n)\) under
  exact KLD-EM updates, which is a tighter check than the previous
  generic "trace decreases" test.

### Documentation

* `critical_review_20260514.md` — itemised response to the audit.
* `plan/proxymix_plan_v0.2_methodological.md` — forward methodological
  plan: v0.2 (`from_kde()` graduation guard-railed), v0.3 (affine-
  Gaussian operator calculus), and the audit-mandated five-phase
  protocol for the collider / DAG research branch.

### Internal

* `gmm_canonicalise()` is the single source of truth for component
  ordering — used by all three fitters and the dispatcher.

## proxymix 0.0.1 (2026-05-13)

Initial development release. Local-only; not yet on CRAN.

### Tier 1 — implemented

* `fit_proxymix()` top-level dispatcher with three fitting regimes:
  `"moment"` (closed-form moment matching), `"sample"` (classical EM on
  i.i.d. samples), and `"kld"` (importance-sampled KLD-EM against an
  evaluable-only target density). The `"auto"` regime picks the cheapest
  applicable regime from the structure of the supplied `gmm_target`.
* S7 class hierarchy: `gmm_target`, `gmm_fit`, `is_proposal`, with
  `print()` / `format()` methods and validators.
* Closed-form GMM operators in `gmm_ops.R`: `dgmm()`, `rgmm()`,
  `gmm_marginalise()`, `gmm_conditionalise()` (Schur complement),
  `gmm_kld()` (Monte Carlo estimator with variational upper / lower
  bounds for sanity).
* Importance-sampling proposals in `proposals.R`: `is_uniform()`,
  `is_mvn()`, `is_mvt()`; all wrap a `is_proposal` instance.
* Diagnostics: `kld_trace()`, `ess_trace()`, `hellinger_mc()`,
  `bic_aic()`.
* Multi-start best-of (Karlis & Xekalaki) initialisation in `init.R`,
  plus `init_random()`, `init_kmeans()`, `init_moment_seed()`.
* Built-in target factories used in the vignettes: `banana_target()`,
  `donut_target()`, `mixture_target()`, plus the from-samples and
  from-function constructors `gmm_target_from_samples()` and
  `gmm_target()`.
* Four vignettes: `quickstart`, `three_regimes`, `density_shapes`
  (the wedge demonstration), and `roadmap` (Tier-2 stubs).

### Tier 2 — provisioned stubs only

The following functions ship with stable signatures, full documentation,
and signature-stability tests; the body raises a "not yet implemented"
condition with a pointer to `vignettes/roadmap.Rmd`.

* `from_kde()` — KDE to GMM proxy via KLD-EM.
* `from_aggregate_likelihood()` — aggregate-likelihood downscaling
  (Sejdinovic et al. kernel-downsizing framework).
* `fit_kld_em_collider()` — KLD-EM under DAG-implied conditional
  independence constraints.
* `to_apsim_scenarios()` — Gaussian-mixture samples to APSIM scenario
  tables.
* `from_simulator()` — wrap an expensive simulator as a
  `gmm_target` via kernel-density or empirical-likelihood
  bridges.

### Tier 3 — deferred (not in scope)

Adaptive importance sampling, variational boosting, normalising-flow
proposals, Stan / INLA inter-operation.
