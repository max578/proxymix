# Monte-Carlo Hellinger distance between a fit and its target

Estimates the squared Hellinger distance
`H^2(f, g) = 1 - integral sqrt(f(x) g(x)) dx` as
`1 - mean(sqrt(f(x) / g(x)))` over `n_mc` fresh draws `x` from the fit
`g`, for every regime. The draws are independent of any sample used to
fit, so the estimate is unbiased, with standard error `se`. The target's
`log_density` must be supplied **and normalised**; otherwise the Monte
Carlo integral is biased by the missing \\\sqrt{Z(f)}\\. When the
target's `normalised` property is not `TRUE`, a warning is issued and
the returned value is flagged.

## Usage

``` r
hellinger_mc(fit, n_mc = 5000L, seed = NULL)
```

## Arguments

- fit:

  A [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)
  whose target carries a `log_density`.

- n_mc:

  Number of Monte Carlo draws from the fit.

- seed:

  Optional integer seed for the draws.

## Value

A list with components

- `h2` - estimate of `H^2(f, g)`,

- `se` - Monte Carlo standard error,

- `n_mc` - number of draws with a finite density ratio,

- `trustworthy` - `TRUE` when the target is declared normalised.

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
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
fit <- fit_proxymix(banana_target(), N = 3L, regime = "kld",
                    is_size = 2000L, max_iter = 25L, seed = 1L)
hellinger_mc(fit, n_mc = 1000L, seed = 1L)
#> $h2
#> [1] 0.00506444
#> 
#> $se
#> [1] 0.003431508
#> 
#> $n_mc
#> [1] 1000
#> 
#> $trustworthy
#> [1] TRUE
#> 
```
