## Interoperability with mice.
##
## A `gmm_imputation` carries a joint Gaussian-mixture imputation -- which can
## be multimodal and heteroscedastic in a way mice's univariate methods are
## not -- and the established mice machinery (`with()`, `pool()`, diagnostics)
## is the natural place to take an arbitrary model estimand from there.
## `as_mids()` packages the completions as a mice multiply-imputed dataset so
## the two compose: impute the joint with proxymix, pool the model with mice.

#' Convert imputations to a mice multiply-imputed dataset
#'
#' Packages a [gmm_imputation] as a `mice::mids` object so that an arbitrary
#' model estimand can be fitted and pooled with the established mice workflow.
#' The joint Gaussian-mixture imputations -- including the multimodal and
#' heteroscedastic shapes a univariate imputer cannot produce -- flow through
#' unchanged; mice supplies `with()`, [mice::pool()], and the pooled
#' diagnostics.
#'
#' @param object A [gmm_imputation].
#'
#' @returns A `mice::mids` object with `m` imputations.
#' @family imputation
#' @seealso [gmm_impute()], [proxy_pool()] for the closed-form column-mean
#'   pooling.
#' @export
#' @examples
#' set.seed(1)
#' x1 <- rnorm(150); x2 <- x1 + rnorm(150)
#' x2[runif(150) < 0.3] <- NA
#' imp <- gmm_impute(cbind(x1, x2), N = 1L, m = 10L, seed = 1L)
#' if (requireNamespace("mice", quietly = TRUE)) {
#'   fit <- with(as_mids(imp), lm(x2 ~ x1))
#'   summary(mice::pool(fit))
#' }
as_mids <- function(object) {
  if (!S7::S7_inherits(object, gmm_imputation)) {
    cli::cli_abort("`object` must be a {.cls gmm_imputation}.")
  }
  if (!requireNamespace("mice", quietly = TRUE)) {
    cli::cli_abort("Install the {.pkg mice} package to use {.fn as_mids}.")
  }
  n <- nrow(object@data)
  base <- as.data.frame(object@data)
  names(base) <- object@var_names
  blocks <- vector("list", object@m + 1L)
  blocks[[1L]] <- cbind(base, .imp = 0L, .id = seq_len(n))
  for (i in seq_len(object@m)) {
    d <- as.data.frame(object@completions[[i]])
    names(d) <- object@var_names
    blocks[[i + 1L]] <- cbind(d, .imp = i, .id = seq_len(n))
  }
  long <- do.call(rbind, blocks)
  mice::as.mids(long, .imp = ".imp", .id = ".id")
}
