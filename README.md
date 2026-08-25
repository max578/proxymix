
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# proxymix

<!-- badges: start -->

<!-- badges: end -->

`proxymix` fits multivariate Gaussian-mixture proxies that are
**Kullback–Leibler optimal** to user-supplied target densities on
$\mathbb{R}^p$. Three regimes from van der Hoek and Elliott (2024) are
unified under one verb:

| Regime | When it applies | Method |
|----|----|----|
| **(i) moment** | $N = 1$ component | Closed-form moment matching |
| **(ii) sample** | i.i.d. samples from the target are available | Classical EM |
| **(iii) kld** | target density `f(x)` can be evaluated but **not** (cheaply) sampled | **KLD-EM with importance sampling** |

Regime (iii) is the reason the package exists. The sample-based mixture
packages (`mclust`, `mixtools`, `flexmix`) all assume i.i.d. draws from
the target; `proxymix` fits directly against an evaluable (possibly
unnormalised) log-density. The nearest CRAN neighbour is `AdMit`, which
adaptively fits a mixture of Student-t distributions to an evaluable
kernel as an importance/proposal density; `proxymix` differs in fitting
a **Gaussian** mixture that is Kullback–Leibler optimal, precisely so
that the fitted object then supports the closed-form operator calculus
(marginals, conditionals, Bayes updates, products, convolutions,
filtering) and carries a fit-quality certificate through every
operation.

**Why not MCMC?** If you can evaluate the unnormalised density you can
always run a sampler and then fit a mixture to the draws. The
regime-(iii) fit is the shortcut when what you want *is* the compact
closed-form object: no chain tuning or convergence diagnostics, a
deterministic pipeline given the seed, and a mixture whose marginals,
conditionals, moments and samples are then available in closed form
through the operator calculus. The trade-off is dimension: the
importance sampling that drives regime (iii) loses effective sample size
sharply beyond roughly $p = 5$–$10$, and every fit reports its effective
sample size so that limit is visible rather than silent. For
high-dimensional posteriors, sample with your favourite MCMC and use
regime (ii) on the draws.

## Installation

From GitHub:

    # install.packages("remotes")
    remotes::install_github("max578/proxymix")

or locally from the source tree:

    R CMD build proxymix
    R CMD INSTALL proxymix_*.tar.gz

Documentation site: <https://max578.github.io/proxymix/>.

## Quick start

``` r
library(proxymix)

## A target you can evaluate but not sample from -- a 2D "banana".
banana <- banana_target()

## Fit a 3-component Gaussian mixture proxy via KLD-EM with importance sampling.
fit <- fit_proxymix(banana, N = 3L, regime = "kld",
                    proposal = proposal_mvt(n_dim = 2L, df = 5),
                    is_size = 2000L, max_iter = 60L, seed = 1L)

print(fit)
#> <gmm_fit>: regime = "kld", K = 3, p = 2
#>   target     : banana
#>   iterations : 37
#>   converged  : TRUE
#>   [1] w = 0.6456, |mu| = 0.3299, tr(Sigma) = 1.3601
#>   [2] w = 0.2473, |mu| = 1.2996, tr(Sigma) = 2.2312
#>   [3] w = 0.1071, |mu| = 1.8882, tr(Sigma) = 4.9044

## Closed-form operations on the fitted mixture.
gmm_marginalise(fit, keep = 1L)
#> <marginalise(kld_em[N=3] on banana)>: K = 3 components in p = 1 dimensions
#>   [1] w = 0.6456, |mu| = 0.1818, tr(Sigma) = 0.4620
#>   [2] w = 0.2473, |mu| = 1.1518, tr(Sigma) = 0.4846
#>   [3] w = 0.1071, |mu| = 1.5954, tr(Sigma) = 0.6321
gmm_conditionalise(fit, given = c(NA, 0.5))
#> <conditionalise(kld_em[N=3] on banana)>: K = 3 components in p = 1 dimensions
#>   [1] w = 0.6730, |mu| = 0.2620, tr(Sigma) = 0.4524
#>   [2] w = 0.2576, |mu| = 1.1133, tr(Sigma) = 0.2352
#>   [3] w = 0.0694, |mu| = 1.4188, tr(Sigma) = 0.1199
```

## The unified fitting verb

    fit_proxymix(
      target,
      N        = 1L,
      regime   = c("auto", "moment", "sample", "kld"),
      ...
    )

- `target` is an S7 `gmm_target` produced by `gmm_target()` (from a
  log-density), `gmm_target_from_samples()` (from samples), or one of
  the built-in factories (`banana_target()`, `donut_target()`,
  `mixture_target()`).
- `regime = "auto"` picks the cheapest applicable regime from the
  target’s contents.

## Closed-form mixture operators

    dgmm(x, fit)              # density evaluation
    rgmm(n, fit)              # exact sampling
    gmm_marginalise(fit, keep)
    gmm_conditionalise(fit, given)
    gmm_kld(fit_p, fit_q)     # Monte Carlo + variational bounds

## Mapping the optima of an objective

    fit   <- from_objective(f, lower, upper, N = 10L)  # a mixture map of the optima
    gmm_modes(fit)$modes                               # the distinct optima

`from_objective()` treats an objective `f` as the Gibbs measure
`exp(-f / T)` – a regime-(iii) target you can evaluate but not sample –
and returns a closed-form mixture over its low regions, so a multimodal
`f` is recovered as a whole rather than one optimum at a time.
`gmm_modes()` resolves the fitted map into the recovered optima.

## Which function do I need?

| You have | You want | Reach for |
|----|----|----|
| an evaluable (unnormalised) log-density | a compact closed-form proxy | `gmm_target()` then `fit_proxymix(regime = "kld")` |
| i.i.d. samples | a mixture fit | `gmm_target_from_samples()` then `fit_proxymix()` |
| an objective function | a map of its optima | `from_objective()`, `gmm_modes()` |
| a kernel density estimate | a small closed-form surrogate | `from_kde()` |
| a fitted mixture | marginals, conditionals, updates | `gmm_marginalise()`, `gmm_conditionalise()`, `gmm_observe()`, `gmm_affine()` |
| data with holes | multiple imputation | `gmm_impute()` with `mar()` / `mnar()` / `censored()` |
| a time series + state-space model | filtering / stability testing | `gmm_filter()`, `gmm_eos_test()` |
| a fitted mixture | information-theoretic diagnostics | `gmm_entropy()`, `gmm_divergence()`, `gmm_mutual_information()` |
| a fitted mixture + treatment data | causal / decision quantities | `gmm_intervene()`, `fit_uplift()`, `proxy_cate()` |
| an unnormalised posterior | the marginal likelihood | `gmm_evidence()` |
| a fitted proxy | error bars on any functional | `gmm_fit_ensemble()`, `proxy_functional_ci()` |
| an evaluable-only target | the component count | `select_N()` |
| a fit collapsing in high dimension | an adaptive proposal | `fit_kld_em(adapt = "pmc")` |

## Vignettes

The package ships with twelve vignettes:

- **`quickstart`** – one-page tour.
- **`posterior_proxy`** – the flagship workflow: a real unnormalised
  Bayesian posterior compressed to a proxy, with the evidence,
  closed-form reads, and bootstrap error bars.
- **`three_regimes`** – a walk-through of regimes (i)–(iii) on toy 2-D
  targets, including the agreement of (i) and (iii) at $N = 1$.
- **`density_shapes`** – the regime-(iii) demonstration: banana, donut,
  three-mixture targets fit by importance-sampled KLD-EM.
- **`operator_calculus`** – closed-form pushforward, Bayesian update,
  aggregation and conditioning on a fitted mixture.
- **`from_kde`** – compressing a kernel density estimate into a
  Gaussian-mixture proxy.
- **`many_methods`** – one fitted mixture in place of regression,
  clustering, PCA and ridge regression.
- **`entropy`** – closed-form entropy, divergence and mutual-information
  diagnostics.
- **`calibration`** – mapping the optima of an objective via its Gibbs
  measure.
- **`missing_data`** – multiple imputation by conditioning the fitted
  mixture (missing at random).
- **`missing_data_mnar`** – imputation under value-dependent missingness
  and censoring, with sensitivity analysis.
- **`end_of_sample`** – testing whether the last few observations of a
  series are consistent with a fitted state-space model.

## Interactive tutorial

A standalone, single-file HTML primer for `proxymix` is hosted at

> **<https://max578.github.io/proxymix-tutorial/>**

The page is targeted at Adelaide University Mathematics and Statistics
undergraduates (Guided depth, default) and postgraduate / PhD readers
(Technical depth, via the in-page dial). It runs in the browser with no
R install required, and works the same way offline by double-click.
Source at
[`max578/proxymix-tutorial`](https://github.com/max578/proxymix-tutorial).

## Reference

Hoek, J. van der and Elliott, R. J. (2024). *Mixtures of multivariate
Gaussians.* Stochastic Analysis and Applications.
[doi:10.1080/07362994.2024.2372605](https://doi.org/10.1080/07362994.2024.2372605).

## Licence

MIT © Max Moldovan.
