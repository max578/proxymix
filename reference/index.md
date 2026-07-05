# Package index

## Top-level interface

The unified fitting verb and its dispatcher.

- [`fit_proxymix()`](https://max578.github.io/proxymix/reference/fit_proxymix.md)
  : Fit a Gaussian-mixture proxy to a target density

## Classes

The S7 class hierarchy.

- [`autoplot.gmm_fit`](https://max578.github.io/proxymix/reference/autoplot.gmm_fit.md)
  : Plot a fitted Gaussian-mixture proxy
- [`glance.gmm_fit`](https://max578.github.io/proxymix/reference/glance.gmm_fit.md)
  : Glance at a fitted Gaussian-mixture proxy
- [`gmm()`](https://max578.github.io/proxymix/reference/gmm.md) : A
  Gaussian mixture
- [`gmm_counterfactual_law()`](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md)
  : A per-unit counterfactual law
- [`gmm_dim()`](https://max578.github.io/proxymix/reference/gmm_dim.md)
  : Dimension of a Gaussian mixture
- [`gmm_fit()`](https://max578.github.io/proxymix/reference/gmm_fit.md)
  : A fitted Gaussian-mixture proxy
- [`gmm_n_components()`](https://max578.github.io/proxymix/reference/gmm_n_components.md)
  : Number of components in a Gaussian mixture
- [`gmm_target()`](https://max578.github.io/proxymix/reference/gmm_target.md)
  : A target density on R^p
- [`gmm_weights()`](https://max578.github.io/proxymix/reference/gmm_weights.md)
  [`gmm_means()`](https://max578.github.io/proxymix/reference/gmm_weights.md)
  [`gmm_covariances()`](https://max578.github.io/proxymix/reference/gmm_weights.md)
  : Component parameters of a Gaussian mixture
- [`is_proposal()`](https://max578.github.io/proxymix/reference/is_proposal.md)
  : An importance-sampling proposal
- [`tidy.gmm`](https://max578.github.io/proxymix/reference/tidy.gmm.md)
  : Tidy a Gaussian mixture into a component table

## Target constructors

Ways to build a `gmm_target` — from samples, from a log-density, or from
a built-in.

- [`banana_target()`](https://max578.github.io/proxymix/reference/banana_target.md)
  : Banana-shaped 2-D target
- [`donut_target()`](https://max578.github.io/proxymix/reference/donut_target.md)
  : Donut-shaped 2-D target
- [`epanechnikov_target()`](https://max578.github.io/proxymix/reference/epanechnikov_target.md)
  : Compact-support Epanechnikov target
- [`gmm_target_from_samples()`](https://max578.github.io/proxymix/reference/gmm_target_from_samples.md)
  : Build a target from samples alone
- [`maxent_target()`](https://max578.github.io/proxymix/reference/maxent_target.md)
  : Maximum-entropy target under moment and support constraints
- [`mixture_target()`](https://max578.github.io/proxymix/reference/mixture_target.md)
  : Three-component Gaussian-mixture target

## Fitting regimes

The three KL-optimal fitting regimes from Hoek and Elliott (2024).

- [`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md)
  : Classical EM fit on samples
- [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md)
  : Importance-sampled KLD-EM fit (regime iii)
- [`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md)
  : Closed-form moment-matching fit
- [`from_kde()`](https://max578.github.io/proxymix/reference/from_kde.md)
  : Compile a kernel-density estimate into a Gaussian-mixture proxy
- [`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md)
  : Map the optima of an objective with a Gaussian-mixture proxy
- [`select_N()`](https://max578.github.io/proxymix/reference/select_N.md)
  : Select the number of mixture components

## Closed-form mixture operators

Density, sampling, marginalisation, conditioning, divergence.

- [`dgmm()`](https://max578.github.io/proxymix/reference/dgmm.md) :
  Density of a Gaussian mixture
- [`gmm_canonicalise()`](https://max578.github.io/proxymix/reference/gmm_canonicalise.md)
  : Canonicalise the component ordering of a Gaussian mixture
- [`gmm_conditionalise()`](https://max578.github.io/proxymix/reference/gmm_conditionalise.md)
  : Conditional of a Gaussian mixture
- [`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md)
  : Divergence between two Gaussian mixtures
- [`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md)
  : Kullback-Leibler divergence between two Gaussian mixtures
- [`gmm_marginalise()`](https://max578.github.io/proxymix/reference/gmm_marginalise.md)
  : Marginal of a Gaussian mixture
- [`gmm_mean()`](https://max578.github.io/proxymix/reference/gmm_mean.md)
  [`gmm_cov()`](https://max578.github.io/proxymix/reference/gmm_mean.md)
  : Mean and covariance of a Gaussian mixture
- [`gmm_modes()`](https://max578.github.io/proxymix/reference/gmm_modes.md)
  : Modes of a Gaussian mixture
- [`pgmm()`](https://max578.github.io/proxymix/reference/pgmm.md)
  [`qgmm()`](https://max578.github.io/proxymix/reference/pgmm.md) :
  Distribution and quantile functions of a one-dimensional mixture
- [`rgmm()`](https://max578.github.io/proxymix/reference/rgmm.md) :
  Sample from a Gaussian mixture

## Affine and causal operators

Pushforward, Kalman update, conditioning, the do-operator and the
counterfactual.

- [`gmm_affine()`](https://max578.github.io/proxymix/reference/gmm_affine.md)
  : Affine pushforward of a Gaussian mixture
- [`gmm_aggregate()`](https://max578.github.io/proxymix/reference/gmm_aggregate.md)
  : Aggregation pushforward of a Gaussian mixture
- [`gmm_convolve()`](https://max578.github.io/proxymix/reference/gmm_convolve.md)
  : Convolution of two independent Gaussian mixtures
- [`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md)
  : Counterfactual law of one unit (abduction, action, prediction)
- [`gmm_filter()`](https://max578.github.io/proxymix/reference/gmm_filter.md)
  : Bounded Gaussian-sum filtering over an observation series
- [`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md)
  : Interventional law of a Gaussian mixture (the do-operator)
- [`gmm_missing()`](https://max578.github.io/proxymix/reference/gmm_missing.md)
  : Condition a Gaussian mixture on the exact values of some coordinates
- [`gmm_mix()`](https://max578.github.io/proxymix/reference/gmm_mix.md)
  : Mix Gaussian mixtures into one mixture
- [`gmm_observe()`](https://max578.github.io/proxymix/reference/gmm_observe.md)
  : Bayesian update of a Gaussian mixture on a noisy linear observation
- [`gmm_product()`](https://max578.github.io/proxymix/reference/gmm_product.md)
  : Pointwise product of two Gaussian mixtures
- [`gmm_reduce()`](https://max578.github.io/proxymix/reference/gmm_reduce.md)
  : Reduce a Gaussian mixture to fewer components

## Decision layer (uplift / next-best-action)

One joint fit read as CATE, optimal action, off-line policy value, and
an identification audit.

- [`fit_uplift()`](https://max578.github.io/proxymix/reference/fit_uplift.md)
  : Fit an uplift / next-best-action model from a data frame

- [`gmm_cf_mean()`](https://max578.github.io/proxymix/reference/gmm_cf_mean.md)
  : The identified counterfactual mean

- [`gmm_cf_tail_prob()`](https://max578.github.io/proxymix/reference/gmm_cf_tail_prob.md)
  : Refused: a tail probability of an individual counterfactual law

- [`gmm_cf_variance()`](https://max578.github.io/proxymix/reference/gmm_cf_variance.md)
  : Refused: the variance of an individual counterfactual law

- [`gmm_counterfactual()`](https://max578.github.io/proxymix/reference/gmm_counterfactual.md)
  : Counterfactual law of one unit (abduction, action, prediction)

- [`gmm_intervene()`](https://max578.github.io/proxymix/reference/gmm_intervene.md)
  : Interventional law of a Gaussian mixture (the do-operator)

- [`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md)
  : Heterogeneous treatment effects (CATE / uplift)

- [`proxy_confounding_gap()`](https://max578.github.io/proxymix/reference/proxy_confounding_gap.md)
  : Confounding gap: the sensitivity of the effect to the latent regime

- [`proxy_decide()`](https://max578.github.io/proxymix/reference/proxy_decide.md)
  : Optimal action and expected incremental value per unit

- [`proxy_identification_report()`](https://max578.github.io/proxymix/reference/proxy_identification_report.md)
  : The identification report (an executive one-pager)

- [`proxy_overlap()`](https://max578.github.io/proxymix/reference/proxy_overlap.md)
  : Per-unit overlap / positivity diagnostic

- [`proxy_policy_value()`](https://max578.github.io/proxymix/reference/proxy_policy_value.md)
  : Off-line value of a targeting policy

- [`proxy_predict()`](https://max578.github.io/proxymix/reference/proxy_predict.md)
  : Predicted outcome under a treatment (the seeing rung)

- [`proxy_regime_segments()`](https://max578.github.io/proxymix/reference/proxy_regime_segments.md)
  : The fitted regimes as an interpretable segment table

- [`proxy_retrospective_uplift()`](https://max578.github.io/proxymix/reference/proxy_retrospective_uplift.md)
  : Retrospective (counterfactual-mean) uplift for observed units

- [`proxy_uplift()`](https://max578.github.io/proxymix/reference/proxy_uplift.md)
  :

  Uplift (alias of
  [`proxy_cate()`](https://max578.github.io/proxymix/reference/proxy_cate.md)
  for a binary treatment)

- [`uplift_identification()`](https://max578.github.io/proxymix/reference/uplift_identification.md)
  : Identification-report object

- [`uplift_model()`](https://max578.github.io/proxymix/reference/uplift_model.md)
  : A fitted uplift / next-best-action model

## Missing data

Multiple imputation by mixture conditioning, with pooling and the
fraction of missing information.

- [`as_mids()`](https://max578.github.io/proxymix/reference/as_mids.md)
  : Convert imputations to a mice multiply-imputed dataset

- [`gmm_complete()`](https://max578.github.io/proxymix/reference/gmm_complete.md)
  :

  Extract completed datasets from a `gmm_imputation`

- [`gmm_imputation()`](https://max578.github.io/proxymix/reference/gmm_imputation.md)
  : A Gaussian-mixture multiple-imputation result

- [`gmm_impute()`](https://max578.github.io/proxymix/reference/gmm_impute.md)
  : Multiple imputation by Gaussian-mixture conditioning

- [`mar()`](https://max578.github.io/proxymix/reference/mechanism.md)
  [`mnar()`](https://max578.github.io/proxymix/reference/mechanism.md)
  [`censored()`](https://max578.github.io/proxymix/reference/mechanism.md)
  : Missingness mechanisms for multiple imputation

- [`proxy_fmi()`](https://max578.github.io/proxymix/reference/proxy_fmi.md)
  : Fraction of missing information for a column mean

- [`proxy_mnar_sensitivity()`](https://max578.github.io/proxymix/reference/proxy_mnar_sensitivity.md)
  : Missing-not-at-random sensitivity analysis for a coordinate mean

- [`proxy_pool()`](https://max578.github.io/proxymix/reference/proxy_pool.md)
  : Pool a column mean across imputations

## State-space instability testing

End-of-sample structural-break test on the Gaussian-sum filter.

- [`gmm_eos_test()`](https://max578.github.io/proxymix/reference/gmm_eos_test.md)
  : End-of-sample instability test on a Gaussian state-space filter

## Importance-sampling proposals

Pluggable proposals for regime (iii).

- [`is_mvn()`](https://max578.github.io/proxymix/reference/is_mvn.md) :
  Multivariate-normal proposal
- [`is_mvt()`](https://max578.github.io/proxymix/reference/is_mvt.md) :
  Multivariate-t proposal
- [`is_uniform()`](https://max578.github.io/proxymix/reference/is_uniform.md)
  : Uniform-on-a-box proposal
- [`proposal_uniform()`](https://max578.github.io/proxymix/reference/proposal_uniform.md)
  [`proposal_mvn()`](https://max578.github.io/proxymix/reference/proposal_uniform.md)
  [`proposal_mvt()`](https://max578.github.io/proxymix/reference/proposal_uniform.md)
  : Preferred names for the importance-proposal constructors

## Interoperability

Compiling external Bayesian posteriors into proxies.

- [`gmm_target_from_posterior()`](https://max578.github.io/proxymix/reference/gmm_target_from_posterior.md)
  :

  Compile an unnormalised Bayesian posterior into a `gmm_target`

## Initialisation

Initial guesses + best-of-multistart wrapper.

- [`init_kmeans()`](https://max578.github.io/proxymix/reference/init_kmeans.md)
  : k-means initialisation
- [`init_moment_seed()`](https://max578.github.io/proxymix/reference/init_moment_seed.md)
  : Moment-seed initialisation
- [`init_random()`](https://max578.github.io/proxymix/reference/init_random.md)
  : Random initialisation
- [`init_warm_start()`](https://max578.github.io/proxymix/reference/init_warm_start.md)
  : Warm-start initialisation from an existing fit
- [`multi_start_best_of()`](https://max578.github.io/proxymix/reference/multi_start_best_of.md)
  : Multi-start best-of wrapper

## Diagnostics

KLD trace, ESS, Hellinger MC, information criteria.

- [`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md)
  : Information criteria: BIC, AIC, and ICL
- [`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md)
  : Summary of importance-sampling diagnostics
- [`ess_trace()`](https://max578.github.io/proxymix/reference/ess_trace.md)
  : Effective sample size of the importance-sampling weights
- [`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md)
  : Phase-transition component discovery by deterministic annealing
- [`gmm_conditional_entropy()`](https://max578.github.io/proxymix/reference/gmm_conditional_entropy.md)
  : Conditional predictive entropy of a Gaussian mixture
- [`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md)
  : Differential entropy of a Gaussian mixture
- [`gmm_evidence()`](https://max578.github.io/proxymix/reference/gmm_evidence.md)
  : Estimate the target's normalising constant from a fitted proxy
- [`gmm_fit_ensemble()`](https://max578.github.io/proxymix/reference/gmm_fit_ensemble.md)
  : Bootstrap ensemble of a fitted proxy
- [`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md)
  : The quality certificate of a fit or derived mixture
- [`gmm_independence_graph()`](https://max578.github.io/proxymix/reference/gmm_independence_graph.md)
  : Conditional-independence (Gaussian graphical model) structure of a
  mixture
- [`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md)
  : Cauchy-Schwarz mutual information between two coordinate blocks
- [`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md)
  : Monte-Carlo Hellinger distance between a fit and its target
- [`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md)
  : Per-iteration KLD trace of a fit
- [`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)
  : Percentile interval for any functional of a fitted proxy
