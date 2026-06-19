
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# proxymix

<!-- badges: start -->

<!-- badges: end -->

`proxymix` fits multivariate Gaussian-mixture proxies that are
**Kullback–Leibler optimal** to user-supplied target densities on
$\mathbb{R}^p$. Three regimes from Hoek and Elliott (2024) are unified
under one verb:

| Regime | When it applies | Method |
|----|----|----|
| **(i) moment** | $N = 1$ component | Closed-form moment matching |
| **(ii) sample** | i.i.d. samples from the target are available | Classical EM |
| **(iii) kld** | target density `f(x)` can be evaluated but **not** (cheaply) sampled | **KLD-EM with importance sampling** — the wedge |

No other CRAN package implements regime (iii) for arbitrary parametric
$f$. `mclust`, `mixtools`, `flexmix` all assume samples.

## Installation

Local install from the source tree:

    R CMD build proxymix
    R CMD INSTALL proxymix_*.tar.gz

## Quick start

``` r
library(proxymix)

## A target you can evaluate but not sample from — a 2D "banana".
banana <- banana_target()

## Fit a 3-component Gaussian mixture proxy via KLD-EM with importance sampling.
fit <- fit_proxymix(banana, N = 3L, regime = "kld",
                    proposal = is_mvt(n_dim = 2L, df = 5),
                    is_size = 2000L, max_iter = 25L)

print(fit)
#> <gmm_fit>: regime = "kld", K = 3, p = 2
#>   target     : banana
#>   iterations : 25
#>   converged  : FALSE
#>   [1] w = 0.5185, |mu| = 0.3203, tr(Sigma) = 1.4429
#>   [2] w = 0.2664, |mu| = 1.0130, tr(Sigma) = 2.4879
#>   [3] w = 0.2151, |mu| = 1.1665, tr(Sigma) = 2.4686

## Closed-form operations on the fitted mixture.
gmm_marginalise(fit, keep = 1L)
#> <marginal[1] of kld_em[N=3] on banana>: K = 3 components in p = 1 dimensions
#>   [1] w = 0.5185, |mu| = 0.0075, tr(Sigma) = 0.4566
#>   [2] w = 0.2664, |mu| = 0.9980, tr(Sigma) = 0.5020
#>   [3] w = 0.2151, |mu| = 1.0880, tr(Sigma) = 0.5608
gmm_conditionalise(fit, given = c(NA, 0.5))
#> <conditional of kld_em[N=3] on banana>: K = 3 components in p = 1 dimensions
#>   [1] w = 0.5223, |mu| = 0.0245, tr(Sigma) = 0.4551
#>   [2] w = 0.2590, |mu| = 1.1247, tr(Sigma) = 0.2020
#>   [3] w = 0.2188, |mu| = 1.1205, tr(Sigma) = 0.2423
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
`exp(-f / T)` — a regime-(iii) target you can evaluate but not sample — and
returns a closed-form mixture over its low regions, so a multimodal `f` is
recovered as a whole rather than one optimum at a time. `gmm_modes()`
resolves the fitted map into the recovered optima.

## Vignettes

The package ships with eight vignettes:

- **`quickstart`** — one-page tour.
- **`three_regimes`** — a walk-through of regimes (i)–(iii) on toy 2-D
  targets, including the agreement of (i) and (iii) at $N = 1$.
- **`density_shapes`** — the wedge demonstration: banana, donut,
  three-mixture targets fit by importance-sampled KLD-EM.
- **`operator_calculus`** — closed-form pushforward, Bayesian update,
  aggregation and conditioning on a fitted mixture.
- **`from_kde`** — compressing a kernel density estimate into a
  Gaussian-mixture proxy.
- **`many_methods`** — one fitted mixture in place of regression,
  clustering, PCA and ridge regression.
- **`entropy`** — closed-form entropy, divergence and mutual-information
  diagnostics.
- **`calibration`** — mapping the optima of an objective via its Gibbs
  measure.

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
