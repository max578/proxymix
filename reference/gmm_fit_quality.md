# The quality certificate of a fit or derived mixture

Returns the fit-quality certificate: a small list recording the fitting
regime, convergence, degeneracy, the effective-sample-size profile of
the importance weights, and (when a validation split was drawn) the
held-out validation gap. The certificate is stamped into the object's
metadata at fit time and carried through every closed-form operator: a
unary operator (marginal, conditional, affine map, Bayes update)
preserves it unchanged, and an operator with several operands (product,
convolution, mixing) combines their certificates conservatively –
worst-case of each field, with the individual per-operand certificates
retained under `quality_sources` and `regime` set to `"composite"`. It
can therefore be read off a marginal, a conditional, a product, a
filtered belief, or any other derived mixture – alongside the
`provenance` vector recording the chain of operations that produced it.

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
`kld_final`, `kld_approx`, and `validation_gap` (fields not applicable
to the regime are `NA`), or `NULL` for a mixture that was never fitted
(e.g. built directly with
[`gmm()`](https://max578.github.io/proxymix/reference/gmm.md)). A
certificate produced by an operator over several operands has
`regime = "composite"`, the conservative worst-case value of each
numeric field across the operands, and an additional `quality_sources`
element holding the operands' own certificates. `kld_final` is estimated
on the importance draws the fit was tuned to, so it reads low;
`validation_gap` is the held-out `validation_kld` (see
[`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md))
minus `kld_final`, and is `NA` when the fit was made with
`validation_size = 0`. `kld_approx` is half the importance-weighted
variance of `log f - log g` over the fitting draws. To second order it
equals the KL divergence from the target to the proxy, in nats. It does
not depend on the target's normalising constant, so it is available when
`kld_final` is shifted. It is `NA` outside regime `"kld"`.

## Details

Downstream verbs read the same certificate and raise a one-shot advisory
(class `proxymix_low_quality`) when the source fit is flagged. A fit is
flagged when it did not converge, when it is degenerate (fewer effective
importance draws than `min_ess`), or when `kld_approx` exceeds 0.3.
Relative ESS is not used: with a fixed proposal it describes the
proposal, not the fit, and stays the same however good the fit is.

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
#> $kld_approx
#> [1] 0.1029181
#> 
#> $validation_gap
#> [1] -0.008892188
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
#> $kld_approx
#> [1] 0.1029181
#> 
#> $validation_gap
#> [1] -0.008892188
#> 
```
