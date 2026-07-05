# Multivariate-t proposal

Builds an
[is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md)
using a multivariate-Student-t density and sampler with `df` degrees of
freedom, location `mean`, and scale matrix `sigma`. Heavier tails than
[`is_mvn()`](https://max578.github.io/proxymix/reference/is_mvn.md), so
often a safer importance proposal at moderate dimensions.

## Usage

``` r
is_mvt(n_dim, mean = rep(0, n_dim), sigma = diag(n_dim), df = 5)
```

## Arguments

- n_dim:

  Ambient dimension `p`.

- mean:

  Length-`p` numeric location vector. Defaults to the zero vector.

- sigma:

  A `p`-by-`p` symmetric positive-definite scale matrix. Defaults to the
  identity.

- df:

  Degrees of freedom (`df > 2` recommended for finite variance).

## Value

An
[is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md)
object.

## See also

Other proposals:
[`is_mvn()`](https://max578.github.io/proxymix/reference/is_mvn.md),
[`is_uniform()`](https://max578.github.io/proxymix/reference/is_uniform.md),
[`proposal_uniform()`](https://max578.github.io/proxymix/reference/proposal_uniform.md)

## Examples

``` r
q <- is_mvt(n_dim = 2L, df = 5)
q
#> <is_proposal>: "is_mvt" in p = 2 dimensions
```
