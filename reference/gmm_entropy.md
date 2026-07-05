# Differential entropy of a Gaussian mixture

Computes the differential entropy of a Gaussian mixture. The quadratic
(order-2) Renyi entropy \\H_2(g) = -\log \int g(x)^2 \\ dx\\ is
available in closed form, because \\\int g^2\\ is a finite sum of
Gaussian-density evaluations. Shannon entropy has no closed form for a
mixture (the integrand carries the logarithm of a sum) and is estimated
by Monte Carlo, reported with its standard error and an analytic upper
bound that brackets it from above.

## Usage

``` r
gmm_entropy(g, order = c("renyi2", "shannon"), n_mc = 5000L, seed = NULL)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  object.

- order:

  `"renyi2"` (closed-form quadratic Renyi entropy, the default) or
  `"shannon"` (Monte-Carlo estimate with an analytic upper bound).

- n_mc:

  Number of Monte Carlo samples for `order = "shannon"`.

- seed:

  Optional integer seed for the Monte Carlo draw.

## Value

For `order = "renyi2"`, a numeric scalar. For `order = "shannon"`, a
list with components `mc` (the estimate), `mc_se` (its standard error),
`upper_bound` (the analytic upper bound), and `n_mc`.

## See also

[`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md),
[`gmm_kld()`](https://max578.github.io/proxymix/reference/gmm_kld.md)

Other diagnostics:
[`bic_aic()`](https://max578.github.io/proxymix/reference/bic_aic.md),
[`ess_summary()`](https://max578.github.io/proxymix/reference/ess_summary.md),
[`ess_trace()`](https://max578.github.io/proxymix/reference/ess_trace.md),
[`gmm_anneal_path()`](https://max578.github.io/proxymix/reference/gmm_anneal_path.md),
[`gmm_conditional_entropy()`](https://max578.github.io/proxymix/reference/gmm_conditional_entropy.md),
[`gmm_evidence()`](https://max578.github.io/proxymix/reference/gmm_evidence.md),
[`gmm_fit_ensemble()`](https://max578.github.io/proxymix/reference/gmm_fit_ensemble.md),
[`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md),
[`gmm_independence_graph()`](https://max578.github.io/proxymix/reference/gmm_independence_graph.md),
[`gmm_mutual_information()`](https://max578.github.io/proxymix/reference/gmm_mutual_information.md),
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
g <- gmm(weights = c(0.5, 0.5),
         means = list(c(-2, 0), c(2, 0)),
         covariances = list(diag(2), diag(2)))
gmm_entropy(g)
#> [1] 3.206021
gmm_entropy(g, order = "shannon", n_mc = 2000L, seed = 1L)
#> $mc
#> [1] 3.433183
#> 
#> $mc_se
#> [1] 0.02004944
#> 
#> $upper_bound
#> [1] 3.531024
#> 
#> $n_mc
#> [1] 2000
#> 
```
