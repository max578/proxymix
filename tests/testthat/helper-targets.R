## Helpers shared across test files.
## Compact toy targets used to exercise the three regimes cheaply.

.helper_iid_normal <- function(n = 200, p = 2L, seed = 1L) {
  withr::with_seed(seed, {
    mu <- rep(0, p)
    sigma <- diag(p)
    matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  })
}

.helper_iid_two_component <- function(n = 400, seed = 1L) {
  withr::with_seed(seed, {
    z <- stats::rbinom(n, size = 1L, prob = 0.4)
    mu1 <- c(-2, 0)
    mu2 <- c(2, 0)
    s <- 0.7
    out <- matrix(0, nrow = n, ncol = 2L)
    n1 <- sum(z == 0L)
    n2 <- n - n1
    out[z == 0L, ] <- matrix(stats::rnorm(n1 * 2L, mean = 0, sd = s), ncol = 2L) +
      matrix(rep(mu1, each = n1), ncol = 2L)
    out[z == 1L, ] <- matrix(stats::rnorm(n2 * 2L, mean = 0, sd = s), ncol = 2L) +
      matrix(rep(mu2, each = n2), ncol = 2L)
    out
  })
}

.helper_banana_logf <- function() {
  function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    z1 <- x[, 1L]
    z2 <- x[, 2L] - 0.5 * (z1^2 - 1)
    -0.5 * (z1^2 + z2^2) - log(2 * pi)
  }
}
