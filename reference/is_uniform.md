# Uniform-on-a-box proposal

Builds an
[is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md)
that samples uniformly on the hyperrectangle \\\[\text{lower}\_1,
\text{upper}\_1\] \times \cdots \times \[\text{lower}\_p,
\text{upper}\_p\]\\ and reports the (constant) log-density on that box.
Outside the box the log-density is `-Inf`.

## Usage

``` r
is_uniform(n_dim, lower = -1, upper = 1)
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

## Value

An
[is_proposal](https://max578.github.io/proxymix/reference/is_proposal.md)
object.

## See also

Other proposals:
[`is_mvn()`](https://max578.github.io/proxymix/reference/is_mvn.md),
[`is_mvt()`](https://max578.github.io/proxymix/reference/is_mvt.md),
[`proposal_uniform()`](https://max578.github.io/proxymix/reference/proposal_uniform.md)

## Examples

``` r
q <- is_uniform(n_dim = 2L, lower = -5, upper = 5)
q
#> <is_proposal>: "is_uniform" in p = 2 dimensions
q@sample(3L)
#>           [,1]      [,2]
#> [1,] -3.506990 0.8660913
#> [2,] -1.363172 0.3978948
#> [3,]  2.871068 3.7388601
```
