# proxymix 0.13.0

### New features

* **The algebra is completed.** Four operations whose mathematics already
  lived inside the package are now first-class operators:
  `gmm_product()` (the normalised pointwise product of two mixtures -- the
  conjugate Bayes update, with the marginal evidence returned as
  `metadata$log_integral` and an optional `reduce` to cap the `K1 * K2`
  component growth), `gmm_convolve()` (the exact distribution of the sum of
  independent mixture variables), `gmm_mix()` (model averaging /
  mixture-of-mixtures), and `gmm_mean()` / `gmm_cov()` (the exact first two
  moments). `pgmm()` / `qgmm()` complete the d/p/q/r quartet in one
  dimension. All follow the operator metadata policy: the quality
  certificate travels and the provenance chain grows.
* **`gmm_evidence()`: the normalising constant as a first-class output.**
  With the fitted proxy as the importance proposal, `log Z` of the target
  (the log marginal likelihood, for a posterior handed over as
  `likelihood x prior`) is estimated in the log domain with a delta-method
  standard error and a heavy-tail diagnostic (classed warning
  `proxymix_heavy_tail` when the proxy's tails look lighter than the
  target's). Validated against known constants for unimodal and multimodal
  targets.
* **Accessors and tidiers.** `gmm_weights()`, `gmm_means()`,
  `gmm_covariances()` read the component parameters without reaching into
  the `@` property layout; broom-style `tidy()` (component table) and
  `glance()` (one-row fit summary) register against the `generics` package
  when installed.
* `proposal_uniform()`, `proposal_mvn()`, `proposal_mvt()` are the
  preferred names of the importance-proposal constructors (the historical
  `is_*` prefix reads as a logical predicate); the `is_*` names remain as
  aliases.

### API changes

* The four planned-interface placeholders (`from_aggregate_likelihood()`,
  `fit_kld_em_collider()`, `to_apsim_scenarios()`, `from_simulator()`) and
  the posterior-producer seam helpers (`from_fb_posterior()`,
  `fb_log_posterior_spec()`, `fb_producer_available()`,
  `mock_fb_posterior()`) are internal as of this release: an exported
  function whose only behaviour is to error, or whose contract awaits an
  unreleased counterpart, does not belong on the public surface. Their
  signatures and behaviour are unchanged and remain under test; the
  general-purpose S3 generic `gmm_target_from_posterior()` remains the
  public route for external posteriors.
* The imputation mechanism layer is explicitly sealed: `.as_gate()` rejects
  unknown gate types with a clear error instead of letting a third-party
  gate object fail deep inside the engine.

### Behaviour changes

* `gmm_impute()` on complete data now does what its warning says: the
  completions are the data verbatim, with no bootstrap refits, no
  imputation draws, and no random-number consumption (previously it fitted
  and drew anyway).
* `proxy_pool(method = "analytic")` on a gated (MNAR / censored) mechanism
  announces the downgrade to Rubin's rules instead of switching silently.

### Performance

* The `gmm` validator takes a Cholesky-first fast path and pays for an
  eigendecomposition only when the Cholesky fails -- mixtures are
  constructed inside hot loops (every filter step builds several).

# proxymix 0.12.0

### New features

* **Fit-quality certificate.** Every fitter now stamps a quality certificate
  into the result's metadata (regime, convergence, degeneracy, the
  effective-sample-size profile including a new per-component ESS, the
  support fraction, and the held-out validation gap), readable with the new
  `gmm_fit_quality()`. Every closed-form operator carries the certificate
  through unchanged, together with a `provenance` vector recording the chain
  of operations, so the certificate can be read off a marginal, a filtered
  belief, or any other derived mixture. Downstream verbs (`gmm_entropy()`,
  `gmm_mutual_information()`, `gmm_independence_graph()`, `gmm_intervene()`,
  `gmm_counterfactual()`, `gmm_filter()`) read the same certificate and raise
  a one-shot advisory (class `proxymix_low_quality`) when the source fit is
  flagged.
* **Degeneracy is a state, not a footnote.** An importance-sampling collapse
  (`ESS < min_ess`) now flags the fit as degenerate: the warning is classed
  (`proxymix_low_ess`), `converged` is forced to `FALSE`, and the new
  `on_low_ess = "abort"` refuses to return the fit at all (classed error
  `proxymix_degenerate_fit`). Previously a fit sitting on one effective draw
  could report `converged = TRUE` with a warning as the only trace.
* **Held-out validation on by default.** `validation_size` now defaults to
  `ceiling(is_size / 4)` rather than `0`, so the overfit-vs-generalise
  diagnostic exists on every regime-(iii) fit. Pass `validation_size = 0L`
  to disable.
* `autoplot()` is now registered for plain `gmm` objects too, so
  operator-calculus results (marginals, conditionals, filtered beliefs) plot
  directly rather than only freshly fitted proxies.
* A dimension disclosure now lives at the core fitter: `fit_kld_em()` notes
  `p > 5` and warns (classed `proxymix_high_dimension`) beyond `p = 10`,
  instead of only the wrapper entry points knowing the scaling story.
* Classed conditions throughout the fitting path: `proxymix_low_ess`,
  `proxymix_degenerate_fit`, `proxymix_high_dimension`, `proxymix_support`,
  `proxymix_nonmonotone`, `proxymix_low_quality` -- pipelines can
  condition-handle instead of matching message text.

### Behaviour changes

* **Seeded fits are now reproducible end-to-end.** `fit_kld_em(seed =)`
  previously seeded only the importance-sampling draw; the initialisation
  resample, the kmeans pass, and empty-component reseeds consumed the ambient
  random-number stream, so two calls with the same seed could return
  different fits. All internal draws now derive from the seed. Fits under a
  given seed therefore differ from 0.11.x. `fit_em_samples(seed =)` likewise
  drives its multi-start (previously hard-coded restart seeds).
* `from_objective()` derives a distinct seed per cooling step: previously the
  importance-sampling stream was re-seeded identically at every temperature,
  so with the default `exploration = 0.5` half of the draws were
  byte-identical across the whole ladder and a basin missed by the first
  exploration draw was never probed again.
* **Data-scaled ridge in the fitters.** The EM fitters scale their
  covariance ridge by the data's covariance scale (previously an absolute
  constant), so the same default regularises identically at data scale
  `1e-8` and `1e+8`. The floor is constant within a fit -- a ridge relative
  to the component's own diagonal would shrink together with a collapsing
  component and stop flooring exactly when needed. At unit data scale the
  behaviour is unchanged to first order; the operators' tiny hygiene ridge
  is unchanged.
* An empty (dead) EM component now has its weight reset alongside its mean
  and covariance (previously the weight stayed at zero, the reseeded mean
  was unreachable, and the reseed re-fired every iteration); the reseed
  covariance is at data scale. Both EM fitters also warn (classed
  `proxymix_nonmonotone`) if their objective decreases beyond numerical
  tolerance, which the documentation previously claimed and the code did
  not do.

### Bug fixes

* Mixture-reduction merge costs are computed in the log domain: the raw
  Gaussian-product density underflows for large dimension times scale
  (reaching exactly zero near `p = 115` at scale `1e5`), which made every
  merge cost `0/Inf/NaN` and the merge order arbitrary.
* `hellinger_mc()` reports the standard error of the self-normalised
  weighted estimator; the previous `sd(ratio)/sqrt(n)` treated skewed
  importance weights as uniform and could understate the Monte Carlo error
  by orders of magnitude.

### Performance

* Greedy mixture reduction caches pairwise merge costs and re-costs only the
  pairs touching the merged component (previously all `O(K^2)` pairs were
  re-costed, each with its own Cholesky, after every single merge -- the
  bottleneck inside long Gaussian-sum-filter runs).
* The log-sum-exp row maximum is computed by a vectorised `pmax` reduce
  rather than a per-row loop, on every E-step's hot path.

### Tests

* The shipped conformance case set now executes during test runs (off CRAN)
  through its built-in driver, so the two-sided contract sweep gates checks
  rather than living only in an external harness; the two remaining stubs
  gained negative cases.
* New independent oracles: closed-form Gaussian KL and quadrature for
  `gmm_kld()` at `K = 1` (the package's namesake divergence previously had
  only positivity checks), closed-form Gaussian Hellinger for
  `hellinger_mc()`, quadrature for the logit-link gated moments (the shipped
  default link previously had no moment-level oracle), and a probit-link
  end-to-end MNAR recovery.
* New metamorphic invariants: Renyi-2 affine equivariance
  (`H(AX + b) = H(X) + log|det A|`), marginal-vs-joint quadrature
  consistency, `gmm_observe` with a near-noiseless selection row against
  `gmm_conditionalise`, post-operator mass preservation, and an
  `rgmm`/`dgmm` Kolmogorov-Smirnov cross-check.
* The uplift Monte Carlo studies are skipped on CRAN and grade monotonicity
  with slack plus a rank correlation, instead of strict positivity of every
  adjacent difference from a local-optimum EM.

# proxymix 0.11.6

### Bug fixes

* The regime-(iii) stopping rule is now invariant to the target's normalising
  constant. `fit_kld_em()` previously judged convergence on the relative change
  of the importance-sampled KLD trace, whose magnitude carries the arbitrary
  `-log Z(f)` offset of an unnormalised target -- so the same target, shifted
  by a constant, could stop after two iterations where the unshifted run used
  a hundred. Convergence is now judged on the importance-weighted EM objective
  `Q(theta) = sum_n W_n log g(x_n)`, which never touches the constant.
  Iteration counts (and therefore fits that previously stopped early or late
  for this reason) can change; a regression test pins the invariance.
* `gmm_eos_test()` now resolves its model specification through the same
  machinery as `gmm_filter()`: dynamics and measurement offsets (`b`, `d`) are
  honoured (they were previously ignored without notice, silently changing the
  innovations), function-valued (time-varying) specifications are accepted,
  and Gaussian-sum (mixture) noise is rejected with an informative error since
  both calibrations assume Gaussian innovations.
* The internal log-sum-exp kernel now excludes `NaN` entries from the sum
  (previously they were neutralised only in the row maximum and could still
  poison the result) and returns `+Inf` for a row containing `+Inf` (previously
  `-Inf`, the exact inversion).
* The smooth-gate moment computation in gated (MNAR) imputation floors its
  normaliser exactly as the censored path does, so a component conditional
  sitting deep in the never-missing region can no longer send `NaN` through
  the M-step and abort the EM with an opaque Cholesky error.

### Documentation

* Documentation, NEWS, and shipped validation scripts have been reworded to
  remove internal development vocabulary and references to unreleased
  companion packages; the posterior-producer probe `fb_producer_available()`
  now reads the producer package name from
  `getOption("proxymix.producer_package")` instead of hard-coding one.
* `README` gains a GitHub installation section, a "Why not MCMC?" positioning
  note with the dimensional limits of regime (iii) stated up front, a
  "Which function do I need?" routing table, and the full eleven-vignette
  list; the quick-start example now runs to convergence.
* The missing-data vignette's scope section now reflects the `mnar()` /
  `censored()` mechanisms available since 0.11.0 and cross-references the
  companion vignette.
* `inst/validation/contract-conformance.R` is now self-contained: a minimal
  serial sweep driver ships in the file, so the two-sided case set is
  executable anywhere (an external driver can still be injected via
  `PROXYMIX_CONFORMANCE_ENGINE`).
* The pkgdown reference index gains the missing interoperability section and
  the articles index gains the MNAR vignette (both previously broke the site
  build).

# proxymix 0.11.5

### Fixes

* `withr` is now a hard dependency (`Imports`, previously `Suggests`). It is
  called on the core seeded-sampling paths (annealing, entropy estimation,
  importance-sampled KLD-EM fitting, mixture reduction, the `gmm_target`
  sample constructors, initialisation, and the seeded diagnostics), so a clean
  installation without `withr` would have failed at runtime. Declaring it as the
  runtime dependency it already is removes that failure mode.

### Internal

* Metadata coherence: the `_pkgdown.yml` site `url:` now resolves to the
  package's GitHub Pages host (`max578.github.io/proxymix`), and `CITATION.cff`
  `version` and `date-released` track the released version.
* The censored-imputation recovery check now grades the estimator on its mean
  absolute error across several simulated draws rather than a single draw. With
  more than half of the target column left-censored, finite-sample recovery
  carries a small expected positive bias, so one draw could legitimately sit
  further from the truth than a tight per-draw bound allowed; the comparison
  against the LOD/2 substitute is unchanged.

# proxymix 0.11.4

### Internal

* Ships the durable two-sided conformance case set in `inst/validation/` (the
  authored validation cases an external conformance harness executes to measure
  coverage of the documented function surface). No change to package code; this
  gives the previously local-only validation material a tracked, travelling home.

# proxymix 0.11.3

### New features

* `gmm_independence_graph()` returns the undirected second-order
  conditional-independence (Gaussian graphical model) structure of a fitted
  mixture -- the partial-correlation skeleton of its overall covariance, in
  closed form. Composed with `fit_kld_em()`, it recovers the dependency
  structure of a target that can only be *evaluated* (an unnormalised energy /
  Gibbs density), where no sample exists to drive a sampling-based estimator.
  It is a graphical-model diagnostic, not a causal-discovery method: it recovers
  the undirected Markov skeleton, not edge directions.

# proxymix 0.11.2

### Bug fixes

* `gmm()` now rejects a covariance that is not positive-definite. A component
  covariance with a negative variance (or, more generally, a negative eigenvalue
  beyond numerical tolerance) was previously accepted silently, deferring the
  failure to a later evaluation. The validator now checks the smallest eigenvalue
  of each finite covariance against a scale-relative tolerance, so a clearly
  indefinite matrix is caught at construction while a numerically near-singular
  fit still constructs.

# proxymix 0.11.1

### Bug fixes

* `gmm_impute(mechanism = censored(...))` is now numerically stable when a mixture
  component sits well past the censoring bound. The fit could previously produce
  `NA` responsibilities and abort, because it initialised the missing entries from
  the observed mean -- the wrong side of a one-sided bound -- and the
  truncated-conditional moments could overshoot the interval. It now initialises
  inside the censored interval and clamps the truncated moments to it.

# proxymix 0.11.0

### New features

* **Missing-not-at-random and censored imputation.** `gmm_impute()` now takes a
  missingness mechanism. `mnar()` specifies a selection model in which an entry of
  a coordinate is missing with a probability that depends on its own unobserved
  value, fitted jointly with the mixture so the imputations are not biased by an
  ignorable model. The slope is a sensitivity parameter, supplied rather than
  estimated, and the intercept is calibrated to the observed missingness rate.
  `censored()` handles a known interval, such as a detection limit, by drawing
  from the mixture conditional truncated to that interval in closed form. `mar()`
  is the default and reproduces the earlier behaviour.
* **`proxy_mnar_sensitivity()`** sweeps the sensitivity slope over a grid and
  pools the coordinate mean at each value, the input to a tipping-point analysis.

# proxymix 0.10.0

### New features

* **End-of-sample instability testing.** `gmm_eos_test()` asks whether the last
  `m` observations of a series are consistent with a linear-Gaussian state-space
  model fitted on the rest, in the small-`m` regime (even `m = 1`) where ordinary
  structural-break tests are undefined because the post-break parameters cannot
  be estimated. It scores the last `m` standardised one-step filter innovations
  and calibrates the score either parametrically (a chi-square reference, exact
  under Gaussian innovations) or by the distribution-free subsampling test of
  Andrews (2003), which stays valid under heavy-tailed observation noise. A new
  vignette, *Testing the last observation for instability*, works through a
  local-level example.

# proxymix 0.9.0

### New features

* **Multiple imputation by Gaussian-mixture conditioning.** `gmm_impute()`
  fits a Gaussian mixture to a numeric dataset that contains missing values
  and draws `m` completed datasets from the mixture conditional
  `p(x_missing | x_observed)`, the same Schur-complement algebra as
  `gmm_conditionalise()`. Because the mixture can be multimodal and
  heteroscedastic, the imputations follow the shape of the joint
  distribution, which keeps downstream inference valid on data a
  single-Gaussian or linear-Gaussian imputer mis-specifies. The mixture is
  fitted to the incomplete data by an expectation-maximisation that uses each
  row's observed margin and restores the conditional covariance of the filled
  entries; each completion is drawn under a mixture fitted to a bootstrap
  resample of the rows, so the pooled inference reflects both imputation and
  parameter uncertainty. This release covers numeric data missing at random.

* **Pooling, diagnostics, and mice interoperability.** `gmm_complete()`
  extracts the completed datasets. `proxy_pool()` pools a column mean in
  closed form -- the exact large-sample limit of the between-imputation
  variance, with no Monte-Carlo noise and an imputation / parameter variance
  split -- and `proxy_fmi()` reports its fraction of missing information. For a
  regression or any other model estimand, `as_mids()` packages the completions
  as a `mice` object so the joint mixture imputations flow into `mice::pool()`
  unchanged: proxymix supplies the imputation model, mice the pooling. A new
  vignette, *Imputing missing data with a mixture*, works through a multimodal
  example.

# proxymix 0.8.0

### New features

* **A mixture proxy for the optima of an objective.** `from_objective()`
  fits a Gaussian-mixture proxy to the Gibbs measure `exp(-f(x) / T)` of a
  user-supplied objective `f` over a bounded box, by cooling a short
  temperature ladder through regime-(iii) KLD-EM (`fit_kld_em()`). The Gibbs
  measure can be evaluated point-wise but not directly sampled, so this is
  regime (iii) applied to an objective: the returned mixture is a closed-form
  map over the low regions of `f`, and a multimodal objective is recovered as
  a whole rather than one optimum at a time. A new vignette, *Mapping the
  optima of an objective*, works through the bimodal and Himmelblau cases.

* **Modes of a Gaussian mixture.** `gmm_modes()` returns the distinct local
  modes of a mixture density by Gaussian mean-shift from each component mean
  (Carreira-Perpinan 2000), with the mixture density at each mode. It
  resolves the map from `from_objective()` into the recovered optima and
  applies to any `gmm` or `gmm_fit`.

# proxymix 0.7.0

### New features

* **Closed-form entropy diagnostics.** `gmm_entropy()` returns the closed-form
  quadratic (order-2) Renyi entropy of a Gaussian mixture, or a Monte-Carlo
  Shannon estimate bracketed by an analytic upper bound. `gmm_divergence()`
  returns the closed-form, symmetric Cauchy-Schwarz divergence between two
  mixtures (`type = "cs"`), with `type = "kl"` delegating to `gmm_kld()`. Both
  rest on the Gaussian-product integral, so the closed-form quantities are exact
  finite sums of Gaussian-density evaluations. A new vignette,
  *Entropy diagnostics with proxymix*, demonstrates the layer.

* **Mutual information and conditional predictive entropy.**
  `gmm_mutual_information()` returns the closed-form, non-negative Cauchy-Schwarz
  mutual information between two coordinate blocks of a fitted joint (the
  divergence between the joint and the product of the marginals; zero exactly
  under independence). `gmm_conditional_entropy()` returns the order-2 Renyi
  entropy of the conditional mixture from `gmm_conditionalise()` -- the
  predictive uncertainty of the target coordinates given the conditioned ones,
  evaluated row-by-row.

* **Deterministic-annealing fitting and phase-transition component discovery.**
  `fit_em_samples()` and `fit_kld_em()` gain an opt-in `anneal = TRUE` argument
  that locates the mixture components by deterministic annealing -- cooling a
  temperature from a high value toward one -- before the unchanged cold EM loop
  polishes the fit. The annealed warm-start is markedly less sensitive to local
  optima than a cold multi-start. The companion diagnostic `gmm_anneal_path()`
  tracks the number of distinct centroids as the temperature falls, a
  physics-derived component-count read whose first phase transition has the
  closed-form critical temperature `lambda_max(Sigma^{-1} C)`.

* **Maximum-entropy targets and the ICL criterion.** `maxent_target()`
  constructs the least-committal `gmm_target` consistent with the supplied
  constraints: the Gaussian under second-moment constraints on full support, the
  uniform under a support constraint alone, and a truncated Gaussian under
  second-moment constraints on a box (fit via regime (iii) under the
  automatically selected support-matched proposal). `bic_aic()` now also returns
  the integrated completed likelihood `icl` and the `classification_entropy` of
  the fitted responsibilities (Biernacki, Celeux and Govaert 2000), which
  penalises overlapping components and equals the BIC for a single component.

* **Kernel regression in the unifying-primitive vignette.** The
  *One mixture, many methods* vignette now covers Nadaraya-Watson kernel
  regression as the nonparametric end of the conditioning axis: the conditional
  mean of a one-component-per-datum kernel density estimate equals the
  Nadaraya-Watson estimator exactly, so a single conditioning operation spans
  ordinary least squares (`K = 1`) through fully-local kernel smoothing
  (`K = n`).

* **Kalman filtering over time in the operator-calculus vignette.** The
  *Affine-Gaussian operator calculus* vignette now shows the predict
  (`gmm_affine`) and update (`gmm_observe`) operators run as a filter over a
  time series: at one component the recursion is exactly the classical Kalman
  filter (verified against a textbook implementation on a constant-velocity
  track), and at several components it is the Gaussian-sum filter.

* **Mixture reduction.** `gmm_reduce()` collapses a Gaussian mixture to a budget
  of at most `k_max` components. The default `method = "merge"` is a greedy,
  moment-preserving pairwise merge, using either the Runnalls (2007)
  Kullback-Leibler bound (`cost = "kl"`) or a closed-form Cauchy-Schwarz cost
  (`cost = "cs"`); every merge preserves the combined weight, mean and
  covariance, so the reduced mixture has the same global mean and covariance as
  the original, and reducing to one component returns the moment-matched
  Gaussian. `method = "anneal"` refines the merge with an annealed re-fit of a
  budget-sized proxy and keeps it when it improves on the merge (never worse),
  which helps for smooth, over-parameterised mixtures. Reduction bounds the
  component count of a Gaussian-sum filter built from `gmm_affine()` and
  `gmm_observe()`.

* **Bounded Gaussian-sum filtering.** `gmm_filter()` runs a filter over an
  observation series by alternating the predict operator (`gmm_affine()`), the
  update operator (`gmm_observe()`) and an optional reduction (`gmm_reduce()`).
  At one component it is the Kalman filter (verified against a textbook
  implementation); with a Gaussian-sum process or measurement noise -- a `gmm`
  supplied in place of a covariance matrix -- it is the Gaussian-sum filter of
  Alspach and Sorenson (1972), and the `k_max` cap holds the component count at
  budget over a long horizon. It returns the filtered mixture at each step, the
  filtered means and covariances, and a tidy per-step summary including the
  log marginal evidence. Constant and time-varying dynamics and measurements are
  both supported. A new section of the *Affine-Gaussian operator calculus*
  vignette demonstrates the verb.

# proxymix 0.6.0

### New features

* **New: a closed-form decision layer (uplift / next-best-action).** One
  joint Gaussian-mixture proxy over `(outcome, treatment, covariates)` is read
  -- in closed form, from that single fit -- as prediction, heterogeneous
  treatment effects, optimal per-unit actions, off-line policy value, and an
  identification audit. `fit_uplift()` assembles the joint fit and returns an
  `uplift_model`; the `proxy_*` verbs score it without re-fitting:
  `proxy_predict()` (the response / risk-scoring rung), `proxy_cate()` /
  `proxy_uplift()` (heterogeneous effects with a delta-method or resampling
  standard error), `proxy_decide()` (the revenue-maximising action plus an
  action-flip probability), `proxy_policy_value()` (off-line value of a
  targeting policy), `proxy_confounding_gap()` (the sensitivity to a latent
  confounder), `proxy_retrospective_uplift()` (counterfactual-mean uplift for
  observed units), `proxy_regime_segments()` (the fitted regimes as a segment
  table), `proxy_overlap()` (per-unit positivity / mass coverage), and
  `proxy_identification_report()` (the executive one-pager). The decision layer
  rides only identified quantities -- the conditional average treatment effect
  and counterfactual *means*.
* **New: the do-operator and the counterfactual as first-class operators.**
  `gmm_intervene()` returns the interventional law `p(. | do(.), .)` -- it sets
  the intervened coordinates inside every component without re-weighting the
  regime gate (the graph surgery that distinguishes `do(T = t)` from `T = t`).
  `gmm_counterfactual()` returns the `K`-atom counterfactual law of one observed
  unit by abduction, action, and prediction. Only the counterfactual mean is
  identified: the new `gmm_counterfactual_law` object exposes `gmm_cf_mean()`,
  while `gmm_cf_variance()` and `gmm_cf_tail_prob()` deliberately error
  (`proxymix_not_identified`) -- the individual counterfactual law depends on an
  unidentified cross-world coupling.
* **New: binary outcomes via latent-scale fitting with a discretised
  predictive.** `fit_uplift(outcome_type = "binary")` fits on the latent
  continuous scale; `proxy_cate(scale = "response")` and
  `proxy_predict(scale = "response")` report effects and predictions on the
  discretised predictive probability `P(Y > threshold)`. Count outcomes are
  supported on the same latent-scale reading.

### Documentation

* New vignette *One mixture, many methods* -- using one fitted Gaussian mixture
  in place of regression (`lm`), clustering (`kmeans` / `mclust`), principal
  components (`prcomp`), and ridge regression, with the trade-off of each
  substitution stated. At `N = 1` the conditional mean is exactly least squares,
  the covariance eigenvectors are exactly the principal components, and a
  covariance ridge is exactly the `L2` penalty.
* Vignette figures refreshed to publication quality -- proper mathematical
  subscripts and Greek symbols in titles and axis labels, shared contour levels,
  and a legend on the KDE-vs-proxy comparison.

### Internal and tests

* New `data.table` import for the decision-verb return tables.
* New test files lock the operators against hand-built linear structural causal
  models (`test-intervene.R`), the K = 1 reduction of `proxy_cate()` to the
  ordinary-least-squares treatment coefficient and the asymptotic agreement of
  the delta standard error with the `lm` coefficient standard error
  (`test-uplift-cate.R`), the audit verbs including recovery of a planted latent
  confounder (`test-uplift-audit.R`), and a six-process synthetic validation
  battery graded against known ground truth (`test-uplift-validation.R`).

# proxymix 0.5.0

### New features

* **New: support-aware importance proposals for bounded and one-sided
  targets.** `gmm_target()` gains an optional `support` argument
  (`list(lower = , upper = )`, with `-Inf` / `Inf` for unbounded
  coordinates). When a regime-(iii) fit is given such a target and no
  explicit `proposal`, `fit_kld_em()` now selects a support-matched
  `is_uniform()` proposal automatically -- inset inside a compact box, or
  data-derived for a one-sided coordinate -- instead of the default
  multivariate-t, which placed importance mass where the target log-density
  is `-Inf` and produced non-finite weights. The automatic choice is
  announced with a one-line message, never silently. Unbounded targets are
  unaffected: they keep the heavy-tailed default.
* **New: `epanechnikov_target()`.** A compact-support fixture (the
  Epanechnikov kernel `(3/4)(1 - u^2)` on a box) joining `banana_target()`
  / `donut_target()` / `mixture_target()`. It declares its `support`, so it
  fits via regime (iii) under the auto-selected uniform proposal with no NaN
  weights -- the canonical case where no mixture of full-support Gaussians
  can have compact support.

### Internal and tests

* `inst/validation/regime_iii_pinned_fits.R` gains a pinned Epanechnikov
  bounded-support fit (ESS, support fraction, no NaN weights).
* New regression tests lock the no-NaN-weight guarantee on compact and
  one-sided targets (`test-support-aware-proposal.R`), the exact `K = 1`
  conditional-mean and conditional-variance match against `lm`
  (`test-gmr-k1-lm.R`), and class / constructor / diagnostic contract
  branches. Line coverage raised to >= 90%.

# proxymix 0.4.0

### New features

* **New: a consumer seam for external Bayesian posteriors.**
  `from_fb_posterior()` takes a posterior addressed only through its
  (unnormalised) log-density and returns a closed-form Gaussian-mixture
  proxy via importance-sampled KLD-EM, generalising the input source from a
  kernel-density estimate (`from_kde()`) to any Bayesian posterior. The
  producer interface the seam expects is materialised and documented by
  `fb_log_posterior_spec()`; `fb_producer_available()` is a capability probe
  (degrades to `FALSE`, never errors, when no producer package is
  installed); and `mock_fb_posterior()` is a synthetic producer (known
  Gaussian or banana log-density) for testing the path with no producer
  package present. proxymix never `Imports:` a producer package -- the seam
  is a soft contract and `R CMD check` is clean with none installed.

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
  note: maths, scope constraints (no non-affine fallbacks, no
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

Two methodological extensions that compose cleanly with regime (iii):

### User-visible changes

* **New: `from_kde()`.** Compiles a kernel density
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
  either a registered Bayesian-package method (`brms`,
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
  for class-specific methods registered by Bayesian fitting packages
  (`brms`, `Stan`, `pymc`-via-reticulate) without coupling proxymix to
  any specific Bayesian backend.
* From-KDE log-density evaluation uses chunked matrix builds so that
  peak memory stays bounded for large IS samples.

## proxymix 0.1.1 (2026-05-14)

Scientific hardening pass: regime (iii) is made harder to misuse.

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
  intended as the seed of a growing `inst/validation/` corpus.

### Tests

* New: `test-canonicalise.R`, `test-normalisation.R`,
  `test-validation-split.R`, `test-support-warning.R`, and
  `test-monotone-objective.R`. The last asserts monotonicity of the
  fixed IS-weighted objective \(\sum_n W_n \log g_\theta(x_n)\) under
  exact KLD-EM updates, which is a tighter check than the previous
  generic "trace decreases" test.

### Internal

* `gmm_canonicalise()` is the single source of truth for component
  ordering — used by all three fitters and the dispatcher.

## proxymix 0.0.1 (2026-05-13)

Initial development release. Local-only; not yet on CRAN.

### Implemented

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
  (the regime-(iii) demonstration), and `roadmap` (planned interfaces).

### Provisioned stubs

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

### Deferred (not in scope)

Adaptive importance sampling, variational boosting, normalising-flow
proposals, Stan / INLA inter-operation.
