# Select the number of mixture components

Fits every candidate component count and chooses one by a
regime-appropriate criterion. With samples (regime ii) the choice is the
smallest BIC. With an evaluable-only target (regime iii) each candidate
is scored by its held-out validation KLD – an independent importance
draw the fit never trained on – and the choice follows the
one-standard-error rule: the smallest `N` whose validation score is
within one Monte Carlo standard error of the best score. The scored
table is returned alongside the choice, and callers who prefer to choose
by eye can ignore the recommendation.

## Usage

``` r
select_N(
  target,
  candidates = 1:6,
  regime = c("auto", "sample", "kld"),
  seed = NULL,
  ...
)
```

## Arguments

- target:

  A
  [gmm_target](https://max578.github.io/proxymix/reference/gmm_target.md).

- candidates:

  Integer vector of component counts to try.

- regime:

  `"auto"` (default: `"sample"` when the target carries samples, `"kld"`
  otherwise), `"sample"`, or `"kld"`.

- seed:

  Optional integer seed shared across candidates (paired fits).

- ...:

  Forwarded to the regime's fitter
  ([`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md)
  or
  [`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md)).

## Value

A list of class `proxymix_selection`: `best_n` (the chosen count),
`best_fit` (its
[gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)), and
`table` (one row per candidate with the criterion values).

## See also

Other fitting:
[`fit_em_samples()`](https://max578.github.io/proxymix/reference/fit_em_samples.md),
[`fit_kld_em()`](https://max578.github.io/proxymix/reference/fit_kld_em.md),
[`fit_moment_match()`](https://max578.github.io/proxymix/reference/fit_moment_match.md),
[`from_kde()`](https://max578.github.io/proxymix/reference/from_kde.md),
[`from_objective()`](https://max578.github.io/proxymix/reference/from_objective.md)

## Examples

``` r
sel <- select_N(banana_target(), candidates = 1:3,
                is_size = 1500L, max_iter = 20L, seed = 1L)
sel$table
#>   N validation_kld validation_mc_se chosen
#> 1 1     0.17364347      0.023113220  FALSE
#> 2 2     0.09021690      0.023610948  FALSE
#> 3 3     0.04424349      0.007774635   TRUE
sel$best_n
#> [1] 3
```
