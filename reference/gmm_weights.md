# Component parameters of a Gaussian mixture

Read-only accessors for the component weights, means, and covariances.

## Usage

``` r
gmm_weights(g)

gmm_means(g)

gmm_covariances(g)
```

## Arguments

- g:

  A [gmm](https://max578.github.io/proxymix/reference/gmm.md) (or
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)).

## Value

`gmm_weights()` returns a numeric vector of length `K`; `gmm_means()`
and `gmm_covariances()` return length-`K` lists of length-`p` numeric
vectors and `p` by `p` matrices respectively.

## See also

[`gmm_mean()`](https://max578.github.io/proxymix/reference/gmm_mean.md)
/ [`gmm_cov()`](https://max578.github.io/proxymix/reference/gmm_mean.md)
for the moments of the mixture as a whole,
[`gmm_fit_quality()`](https://max578.github.io/proxymix/reference/gmm_fit_quality.md)
for the fit-quality certificate.

Other classes:
[`autoplot.gmm_fit`](https://max578.github.io/proxymix/reference/autoplot.gmm_fit.md),
[`glance.gmm_fit`](https://max578.github.io/proxymix/reference/glance.gmm_fit.md),
[`gmm()`](https://max578.github.io/proxymix/reference/gmm.md),
[`gmm_counterfactual_law()`](https://max578.github.io/proxymix/reference/gmm_counterfactual_law.md),
[`gmm_dim()`](https://max578.github.io/proxymix/reference/gmm_dim.md),
[`gmm_fit()`](https://max578.github.io/proxymix/reference/gmm_fit.md),
[`gmm_n_components()`](https://max578.github.io/proxymix/reference/gmm_n_components.md),
[`gmm_target()`](https://max578.github.io/proxymix/reference/gmm_target.md),
[`is_proposal()`](https://max578.github.io/proxymix/reference/is_proposal.md),
[`tidy.gmm`](https://max578.github.io/proxymix/reference/tidy.gmm.md)

## Examples

``` r
g <- gmm(weights = c(0.3, 0.7), means = list(-1, 2),
         covariances = list(matrix(1), matrix(0.5)))
gmm_weights(g)
#> [1] 0.3 0.7
gmm_means(g)
#> [[1]]
#> [1] -1
#> 
#> [[2]]
#> [1] 2
#> 
```
