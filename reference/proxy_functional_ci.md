# Percentile interval for any functional of a fitted proxy

Applies a functional to every member of a bootstrap ensemble and returns
the base fit's point estimate with percentile confidence limits. The
functional may return a scalar or a fixed-length numeric vector – the
moments
([`gmm_mean()`](https://max578.github.io/proxymix/reference/gmm_mean.md),
[`gmm_cov()`](https://max578.github.io/proxymix/reference/gmm_mean.md)),
a tail probability via
[`pgmm()`](https://max578.github.io/proxymix/reference/pgmm.md) on a
marginal, an entropy, a conditional mean, or any composition of the
operator calculus.

## Usage

``` r
proxy_functional_ci(ensemble, fn, level = 0.9, ..., tail_warn = 0.001)
```

## Arguments

- ensemble:

  A `gmm_ensemble` from
  [`gmm_fit_ensemble()`](https://max578.github.io/proxymix/reference/gmm_fit_ensemble.md).

- fn:

  A function mapping a
  [gmm](https://max578.github.io/proxymix/reference/gmm.md) to a numeric
  scalar or vector.

- level:

  Confidence level. Default `0.9`.

- ...:

  Forwarded to `fn`.

- tail_warn:

  Threshold for the tail-probability warning described in Details.
  Default `1e-3`; `NULL` or `0` turns the warning off.

## Value

A data frame with one row per element of `fn`'s value: `term`,
`estimate` (the base fit's value), `conf.low`, `conf.high`.

## Details

The interval measures how much the functional moves when the fit is
repeated on re-weighted draws. It does not measure how far the
Gaussian-mixture family itself is from the target, so it can exclude the
target's true value when the mixture has too few components. The gap is
widest far in a tail: a mixture of Gaussians has Gaussian tails and can
misstate a small tail probability by orders of magnitude while the
interval stays narrow. A warning of class `proxymix_tail_functional` is
raised when every value of `fn` (at the base fit and at every member)
lies in `[0, 1]` and some element of the base fit's estimate, or its
complement, is below `tail_warn`. Check such a probability against
numerical integration of the target or a method built for tails.

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
[`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md),
[`gmm_independence_graph()`](https://max578.github.io/proxymix/reference/gmm_independence_graph.md),
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md)

## Examples

``` r
fit <- fit_proxymix(banana_target(), N = 2L, regime = "kld",
                    is_size = 1500L, max_iter = 20L, seed = 1L)
ens <- gmm_fit_ensemble(fit, B = 30L, seed = 2L)
proxy_functional_ci(ens, function(g) gmm_mean(g)[1L])
#>   term     estimate    conf.low  conf.high
#> 1   f1 -0.004937457 -0.08133735 0.07347874
```
