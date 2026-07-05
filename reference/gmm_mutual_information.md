# Cauchy-Schwarz mutual information between two coordinate blocks

Measures the dependence between two disjoint coordinate blocks of a
fitted joint Gaussian mixture as the Cauchy-Schwarz divergence between
the joint over the two blocks and the product of their marginals,
\$\$I\_{\mathrm{CS}}(A; B) = D\_{\mathrm{CS}}(p\_{AB},\\ p_A\\ p_B).\$\$
The product of the marginals is itself a Gaussian mixture, so the
quantity is closed-form. It is non-negative and zero exactly when the
two blocks are independent. (The naive combination \\H_2(A) + H_2(B) -
H_2(A, B)\\ is **not** a valid mutual information: order-2 Renyi
entropies are not additive over independent blocks and that difference
can be negative.)

## Usage

``` r
gmm_mutual_information(g, block_a, block_b)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md))
  joint mixture.

- block_a, block_b:

  Disjoint integer vectors of coordinate indices (in `1..p`) naming the
  two blocks.

## Value

A non-negative numeric scalar.

## See also

[`gmm_entropy()`](https://max578.github.io/proxymix/reference/gmm_entropy.md),
[`gmm_divergence()`](https://max578.github.io/proxymix/reference/gmm_divergence.md)

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
[`hellinger_mc()`](https://max578.github.io/proxymix/reference/hellinger_mc.md),
[`kld_trace()`](https://max578.github.io/proxymix/reference/kld_trace.md),
[`proxy_functional_ci()`](https://max578.github.io/proxymix/reference/proxy_functional_ci.md)

## Examples

``` r
## A correlated bivariate Gaussian: mutual information grows with |rho|.
s <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
g <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(s))
gmm_mutual_information(g, 1L, 2L)
#> [1] 0.102997
```
