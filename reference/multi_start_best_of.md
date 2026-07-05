# Multi-start best-of wrapper

Runs the supplied fitter from each of several initialisations and
returns the fit with the best score, following Karlis and Xekalaki
(2003)'s recommendation.

## Usage

``` r
multi_start_best_of(fit_fn, inits, score_fn, ...)
```

## Arguments

- fit_fn:

  A function with signature `function(init, ...)` returning a
  [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md).

- inits:

  A list of [gmm](https://max578.github.io/proxymix/reference/gmm.md)
  initialisations.

- score_fn:

  A function `function(fit)` returning a numeric score — *larger is
  better* (typically the final log-target evidence).

- ...:

  Additional arguments forwarded to `fit_fn`.

## Value

The [gmm_fit](https://max578.github.io/proxymix/reference/gmm_fit.md)
with the largest `score_fn(fit)`.

## See also

Other init:
[`init_kmeans()`](https://max578.github.io/proxymix/reference/init_kmeans.md),
[`init_moment_seed()`](https://max578.github.io/proxymix/reference/init_moment_seed.md),
[`init_random()`](https://max578.github.io/proxymix/reference/init_random.md),
[`init_warm_start()`](https://max578.github.io/proxymix/reference/init_warm_start.md)

## Examples

``` r
x <- matrix(stats::rnorm(200), ncol = 2)
tgt <- gmm_target_from_samples(x)
inits <- list(init_random(2L, 2L, seed = 1L),
              init_moment_seed(x, N = 2L))
best <- multi_start_best_of(
  fit_fn   = function(init, ...) fit_em_samples(tgt, init = init, ...),
  inits    = inits,
  score_fn = function(fit) fit@diagnostics$loglik_final,
  max_iter = 25L
)
best@diagnostics$loglik_final
#> [1] -285.4343
```
