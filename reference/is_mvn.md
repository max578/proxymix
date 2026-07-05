# Multivariate-normal proposal

Builds an
[is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md)
using a multivariate-normal `N(mean, cov)` density and sampler.

## Usage

``` r
is_mvn(n_dim, mean = rep(0, n_dim), cov = diag(n_dim))
```

## Arguments

- n_dim:

  Ambient dimension `p`.

- mean:

  Length-`p` numeric mean vector. Defaults to the zero vector.

- cov:

  A `p`-by-`p` symmetric positive-definite covariance matrix. Defaults
  to the identity.

## Value

An
[is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md)
object.

## See also

Other proposals:
[`is_mvt()`](https://max578.github.io/proxymix/reference/is_mvt.md),
[`is_uniform()`](https://max578.github.io/proxymix/reference/is_uniform.md),
[`proposal_uniform()`](https://max578.github.io/proxymix/reference/proposal_uniform.md)

## Examples

``` r
q <- is_mvn(n_dim = 2L, mean = c(0, 0), cov = 4 * diag(2))
q
#> <is_proposal>: "is_mvn" in p = 2 dimensions
```
