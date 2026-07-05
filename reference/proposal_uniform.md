# Preferred names for the importance-proposal constructors

`proposal_uniform()`, `proposal_mvn()`, and `proposal_mvt()` are the
preferred names of
[`is_uniform()`](https://max578.github.io/proxymix/reference/is_uniform.md),
[`is_mvn()`](https://max578.github.io/proxymix/reference/is_mvn.md), and
[`is_mvt()`](https://max578.github.io/proxymix/reference/is_mvt.md): the
historical `is_*` prefix reads as a logical predicate, which these
constructors are not. The `is_*` names remain available as aliases and
are not scheduled for removal.

## Usage

``` r
proposal_uniform(n_dim, lower = -1, upper = 1)

proposal_mvn(n_dim, mean = rep(0, n_dim), cov = diag(n_dim))

proposal_mvt(n_dim, mean = rep(0, n_dim), sigma = diag(n_dim), df = 5)
```

## Arguments

- n_dim:

  Ambient dimension `p`.

- lower:

  Length-`p` numeric vector of lower bounds (recycled from a single
  value).

- upper:

  Length-`p` numeric vector of upper bounds (recycled from a single
  value).

- mean, cov:

  Forwarded to
  [`is_mvn()`](https://max578.github.io/proxymix/reference/is_mvn.md).

- sigma, df:

  Forwarded to
  [`is_mvt()`](https://max578.github.io/proxymix/reference/is_mvt.md).

## Value

An
[is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md)
object.

## See also

Other proposals:
[`is_mvn()`](https://max578.github.io/proxymix/reference/is_mvn.md),
[`is_mvt()`](https://max578.github.io/proxymix/reference/is_mvt.md),
[`is_uniform()`](https://max578.github.io/proxymix/reference/is_uniform.md)
