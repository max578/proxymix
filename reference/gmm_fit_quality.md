# The quality certificate of a fit or derived mixture

Returns the fit-quality certificate: a small list recording the fitting
regime, convergence, degeneracy, the effective-sample-size profile of
the importance weights, and (when a validation split was drawn) the
held-out validation gap. The certificate is stamped into the object's
metadata at fit time and carried unchanged through every closed-form
operator, so it can be read off a marginal, a conditional, a filtered
belief, or any other derived mixture – alongside the `provenance` vector
recording the chain of operations that produced it.

## Usage

``` r
gmm_fit_quality(g)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md).

## Value

A list with elements `regime`, `converged`, `degenerate`, `ess`,
`ess_relative`, `min_component_ess`, `max_weight`, `support_fraction`,
`kld_final`, and `validation_gap` (fields not applicable to the regime
are `NA`), or `NULL` for a mixture that was never fitted (e.g. built
directly with
[`gmm()`](https://max578.github.io/proxymix/reference/gmm.md)).

## Details

Downstream verbs read the same certificate and raise a one-shot advisory
(class `proxymix_low_quality`) when the source fit is flagged.

## See also

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md),
[`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md),
[`ess_trace()`](https://max578.github.io/proxymix/reference/ess_trace.md),
[`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md),
[`gmm_conditional_entropy()`](https://max578.github.io/proxymix/reference/gmm_conditional_entropy.md),
[`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md),
[`gmm_evidence()`](https://max578.github.io/proxymix/reference/gmm_evidence.md),
[`gmm_fit_ensemble()`](https://max578.github.io/proxymix/reference/gmm_fit_ensemble.md),
[`gmm_independence_graph()`](https://max578.github.io/proxymix/reference/gmm_independence_graph.md),
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
                    is_size = 1500L, max_iter = 15L, seed = 1L)
gmm_fit_quality(fit)
#> $regime
#> [1] "kld"
#> 
#> $converged
#> [1] FALSE
#> 
#> $degenerate
#> [1] FALSE
#> 
#> $ess
#> [1] 1048.118
#> 
#> $ess_relative
#> [1] 0.6987454
#> 
#> $min_component_ess
#> [1] 165.318
#> 
#> $max_weight
#> [1] 0.005472202
#> 
#> $support_fraction
#> [1] 1
#> 
#> $kld_final
#> [1] 0.144013
#> 
#> $validation_gap
#> [1] -0.05434317
#> 
## The certificate survives the operator calculus.
gmm_fit_quality(gmm_marginalise(fit, keep = 1L))
#> $regime
#> [1] "kld"
#> 
#> $converged
#> [1] FALSE
#> 
#> $degenerate
#> [1] FALSE
#> 
#> $ess
#> [1] 1048.118
#> 
#> $ess_relative
#> [1] 0.6987454
#> 
#> $min_component_ess
#> [1] 165.318
#> 
#> $max_weight
#> [1] 0.005472202
#> 
#> $support_fraction
#> [1] 1
#> 
#> $kld_final
#> [1] 0.144013
#> 
#> $validation_gap
#> [1] -0.05434317
#> 
```
