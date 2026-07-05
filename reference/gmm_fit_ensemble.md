# Bootstrap ensemble of a fitted proxy

Quantifies the sampling variability of the fitted mixture itself by a
Bayesian (weighted) bootstrap: each replicate re-weights the fit's own
observations with Dirichlet(1, ..., 1) weights and refits by a
warm-started weighted EM. In regime `"kld"` the observations are the
fit's cached importance draws with their self-normalised weights, so a
replicate costs **zero new target evaluations**; in regimes `"sample"`
and `"moment"` the observations are the target's samples. Summaries of
any functional of the proxy are then read off the ensemble with
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)
– functional-space intervals sidestep the label-switching that makes
parameter-space intervals incoherent.

## Usage

``` r
gmm_fit_ensemble(fit, B = 200L, max_iter = 25L, tol = 1e-05, seed = NULL)
```

## Arguments

- fit:

  A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)
  whose fitting inputs are recoverable (a regime-`"kld"` fit carries its
  importance draws; regimes `"sample"` and `"moment"` require the target
  to carry its samples).

- B:

  Number of bootstrap replicates.

- max_iter, tol:

  Convergence controls for the per-replicate warm-started weighted EM.

- seed:

  Optional integer seed for the replicate weights.

## Value

A list of class `gmm_ensemble`: `fit` (the base fit), `members` (a
length-`B` list of
[gmm](https://max578.github.io/proxymix/reference/gmm.md)), `B`, and
`regime`.

## References

Rubin, D. B. (1981) The Bayesian bootstrap. *The Annals of Statistics*
9(1), 130–134.

## See also

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md),
[`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md),
[`ess_trace()`](https://max578.github.io/proxymix/reference/ess_trace.md),
[`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md),
[`gmm_conditional_entropy()`](https://max578.github.io/proxymix/reference/gmm_conditional_entropy.md),
[`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md),
[`gmm_evidence()`](https://max578.github.io/proxymix/reference/gmm_evidence.md),
[`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md),
[`gmm_independence_graph()`](https://max578.github.io/proxymix/reference/gmm_independence_graph.md),
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
                    is_size = 1500L, max_iter = 20L, seed = 1L)
ens <- gmm_fit_ensemble(fit, B = 30L, seed = 2L)
proxy_functional_ci(ens, gmm_mean)
#>   term     estimate     conf.low  conf.high
#> 1   f1 -0.004937457 -0.081337352 0.07347874
#> 2   f2  0.074421977 -0.003574058 0.19496631
```
