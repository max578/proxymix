## A mixture proxy for a target given by its characteristic function.
##
## Some distributions have a characteristic function in closed form and no
## usable density: stable laws, sums and time aggregates of Levy increments,
## compound Poisson totals. `from_cf()` fits a Gaussian mixture to such a
## target by least squares between the two characteristic functions on a
## fixed sample of frequencies. `gmm_cf()` evaluates the characteristic
## function of a mixture, and `gmm_cf_distance()` measures a fit along one
## direction by the Kolmogorov distance, computed without sampling.
##
## Inside the fit a mixture is a list with weights `w`, a K-by-p matrix of
## means `m` and a list `L` of lower Cholesky factors, in coordinates in
## which the target has mean zero and unit scale.

# ---------------------------------------------------------------------------
# internal helpers: parameters, loss and its derivatives
# ---------------------------------------------------------------------------

## Frequency at which a standard normal characteristic function equals 1/2.
.cf_half_height <- sqrt(2 * log(2))

## Positions of the lower triangle of a p-by-p matrix (`idx`) and which of
## them lie on the diagonal (`dg`).
.cf_chol_index <- function(p) {
  idx <- which(lower.tri(diag(p), diag = TRUE))
  list(
    idx = idx,
    dg = which(row(diag(p))[idx] == col(diag(p))[idx]),
    n_l = length(idx)
  )
}

## Parameter vector: K - 1 weight logits, the K-by-p means, then the lower
## triangle of each Cholesky factor with its diagonal on the log scale.
.cf_pack <- function(w, m, L) {
  ix <- .cf_chol_index(ncol(m))
  lt <- vapply(L, function(Lk) {
    v <- Lk[ix$idx]
    v[ix$dg] <- log(v[ix$dg])
    v
  }, numeric(ix$n_l))
  c(log(w[-1L] / w[1L]), as.vector(m), as.vector(lt))
}

.cf_unpack <- function(par, K, p) {
  ix <- .cf_chol_index(p)
  a <- c(0, par[seq_len(K - 1L)])
  w <- exp(a - max(a))
  w <- w / sum(w)
  m <- matrix(par[K - 1L + seq_len(K * p)], K, p)
  lt <- matrix(par[K - 1L + K * p + seq_len(K * ix$n_l)], ix$n_l, K)
  L <- vector("list", K)
  for (k in seq_len(K)) {
    v <- lt[, k]
    v[ix$dg] <- exp(v[ix$dg])
    Lk <- matrix(0, p, p)
    Lk[ix$idx] <- v
    L[[k]] <- Lk
  }
  list(w = w, m = m, L = L)
}

## Component characteristic functions at the rows of `W`: an M-by-K complex
## matrix.
.cf_components <- function(W, g) {
  K <- length(g$w)
  C <- matrix(0i, nrow(W), K)
  for (k in seq_len(K)) {
    U <- W %*% g$L[[k]]
    C[, k] <- exp(1i * as.vector(W %*% g$m[k, ]) - 0.5 * rowSums(U^2))
  }
  C
}

.cf_mixture <- function(W, g) as.vector(.cf_components(W, g) %*% g$w)

## Loss (mean squared modulus of the difference between the target values
## `ph` and the mixture's characteristic function at the rows of `W`) and
## its gradient in the packed parameters.
.cf_loss_grad <- function(par, W, ph, K, p) {
  ix <- .cf_chol_index(p)
  g <- .cf_unpack(par, K, p)
  C <- .cf_components(W, g)
  pg <- as.vector(C %*% g$w)
  r <- ph - pg
  rc <- Conj(r) / nrow(W)
  A <- sweep(C, 2L, g$w, "*") * rc
  ga <- -2 * (colSums(Re(A)) - g$w * sum(Re(rc * pg)))
  gm <- 2 * crossprod(Im(A), W)
  gl <- matrix(0, ix$n_l, K)
  for (k in seq_len(K)) {
    Gk <- 2 * crossprod(W * Re(A[, k]), W %*% g$L[[k]])
    v <- Gk[ix$idx]
    v[ix$dg] <- v[ix$dg] * g$L[[k]][ix$idx][ix$dg]
    gl[, k] <- v
  }
  list(
    value = mean(Mod(r)^2),
    gradient = c(ga[-1L], as.vector(gm), as.vector(gl))
  )
}

## The loss as a sum of squares: residuals (real parts, then imaginary
## parts) and, when `jac = TRUE`, their Jacobian.
.cf_resid_jac <- function(par, W, ph, K, p, jac = TRUE) {
  ix <- .cf_chol_index(p)
  g <- .cf_unpack(par, K, p)
  C <- .cf_components(W, g)
  pg <- as.vector(C %*% g$w)
  M <- nrow(W)
  so <- 1 / sqrt(M)
  rho <- (ph - pg) * so
  out <- list(r = c(Re(rho), Im(rho)))
  if (!jac) return(out)
  J <- matrix(0, 2L * M, K - 1L + K * p + K * ix$n_l)
  rw <- row(diag(p))[ix$idx]
  cl <- col(diag(p))[ix$idx]
  for (k in seq_len(K)) {
    ck <- C[, k] * (g$w[k] * so)
    if (k > 1L) {
      d <- -(ck - g$w[k] * pg * so)
      J[, k - 1L] <- c(Re(d), Im(d))
    }
    J[, K - 1L + (seq_len(p) - 1L) * K + k] <- rbind(Im(ck) * W, -Re(ck) * W)
    U <- W %*% g$L[[k]]
    B <- W[, rw, drop = FALSE] * U[, cl, drop = FALSE]
    B[, ix$dg] <- sweep(
      B[, ix$dg, drop = FALSE], 2L, g$L[[k]][ix$idx][ix$dg], "*"
    )
    J[, K - 1L + K * p + (k - 1L) * ix$n_l + seq_len(ix$n_l)] <-
      rbind(Re(ck) * B, Im(ck) * B)
  }
  out$J <- J
  out
}

# ---------------------------------------------------------------------------
# internal helpers: optimisers
# ---------------------------------------------------------------------------

## Box bounds on the packed parameters. `sd_floor` is the smallest allowed
## Cholesky diagonal, in standardised units.
.cf_bounds <- function(K, p, sd_floor) {
  ix <- .cf_chol_index(p)
  lt_lb <- rep(-20, ix$n_l)
  lt_ub <- rep(20, ix$n_l)
  lt_lb[ix$dg] <- log(sd_floor)
  lt_ub[ix$dg] <- log(20)
  list(
    lb = c(rep(-10, K - 1L), rep(-25, K * p), rep(lt_lb, K)),
    ub = c(rep(10, K - 1L), rep(25, K * p), rep(lt_ub, K))
  )
}

.cf_lbfgs <- function(par, W, ph, K, p, lb, ub, maxit) {
  cache <- new.env()
  fn <- function(par) {
    o <- .cf_loss_grad(par, W, ph, K, p)
    cache$par <- par
    cache$g <- o$gradient
    o$value
  }
  gr <- function(par) {
    if (is.null(cache$par) || any(cache$par != par)) fn(par)
    cache$g
  }
  op <- stats::optim(
    pmin(pmax(par, lb), ub), fn, gr,
    method = "L-BFGS-B", lower = lb, upper = ub,
    control = list(maxit = maxit, lmm = 20L, factr = 10, pgtol = 0)
  )
  ## optim code 1 is the iteration limit
  list(par = op$par, value = op$value, hit_limit = op$convergence == 1L)
}

## Levenberg-Marquardt with box bounds, for the final polish of the
## sum-of-squares loss.
.cf_lm <- function(par, W, ph, K, p, lb, ub, maxit, tol = 1e-9) {
  par <- pmin(pmax(par, lb), ub)
  o <- .cf_resid_jac(par, W, ph, K, p)
  f <- sum(o$r^2)
  lam <- 1e-3
  it <- 0L
  hit_limit <- TRUE
  while (it < maxit) {
    it <- it + 1L
    gr <- as.vector(crossprod(o$J, o$r))
    free <- !((par <= lb + 1e-12 & gr > 0) | (par >= ub - 1e-12 & gr < 0))
    H <- crossprod(o$J[, free, drop = FALSE])
    dH <- diag(H)
    dH <- dH + 1e-10 * max(dH)
    ok <- FALSE
    for (tries in seq_len(25L)) {
      A <- H
      diag(A) <- diag(A) + lam * dH
      ## a singular system is retried with more damping
      st <- tryCatch(-solve(A, gr[free]), error = function(e) NULL)
      if (!is.null(st) && all(is.finite(st))) {
        cand <- par
        cand[free] <- pmin(pmax(par[free] + st, lb[free]), ub[free])
        f_new <- sum(.cf_resid_jac(cand, W, ph, K, p, jac = FALSE)$r^2)
        if (is.finite(f_new) && f_new < f) {
          ok <- TRUE
          break
        }
      }
      lam <- lam * 4
    }
    if (!ok) {
      hit_limit <- FALSE
      break
    }
    rel <- (f - f_new) / max(f, 1e-300)
    par <- cand
    f <- f_new
    lam <- max(lam / 3, 1e-12)
    if (rel < tol || f < 1e-28) {
      hit_limit <- FALSE
      break
    }
    o <- .cf_resid_jac(par, W, ph, K, p)
  }
  list(par = par, value = f, hit_limit = hit_limit)
}

## One local optimisation. Frequencies are used in stages, lowest
## bandwidths first; `lad` gives the bandwidth index of each row of `U`.
## Each stage bounds the component standard deviations below by what its
## highest frequency can resolve, which stops a covariance collapsing early.
.cf_local <- function(par, U, ph, lad, bandwidths, K, p, stages,
                      maxit_stage, maxit_lm) {
  for (s in stages) {
    use <- lad <= s
    bd <- .cf_bounds(K, p, 0.3 / (.cf_half_height * bandwidths[s]))
    op <- .cf_lbfgs(
      par, U[use, , drop = FALSE], ph[use], K, p, bd$lb, bd$ub, maxit_stage
    )
    par <- op$par
  }
  if (maxit_lm > 0L && length(par) <= 700L) {
    op <- .cf_lm(
      par, U[use, , drop = FALSE], ph[use], K, p, bd$lb, bd$ub, maxit_lm
    )
  }
  op
}

# ---------------------------------------------------------------------------
# internal helpers: placing components where the fit lacks mass
# ---------------------------------------------------------------------------

## Mixture weights by least squares for fixed components `C`, clipped to
## the simplex.
.cf_ls_weights <- function(C, ph) {
  G <- Re(crossprod(Conj(C), C))
  b <- Re(as.vector(crossprod(Conj(C), ph)))
  w <- tryCatch(
    solve(G + 1e-10 * diag(mean(diag(G)), ncol(C)), b),
    error = function(e) NULL
  )
  if (is.null(w)) w <- rep(1, ncol(C))
  w <- pmax(w, 1e-4)
  w / sum(w)
}

.cf_draw <- function(n, g) {
  p <- ncol(g$m)
  k <- sample.int(length(g$w), n, replace = TRUE, prob = g$w)
  X <- matrix(0, n, p)
  for (j in unique(k)) {
    i <- which(k == j)
    X[i, ] <- matrix(g$m[j, ], length(i), p, byrow = TRUE) +
      matrix(stats::rnorm(length(i) * p), length(i), p) %*% t(g$L[[j]])
  }
  X
}

## Candidate means for a new component: draws around the origin, the
## current means, and draws from the current mixture widened by a half.
.cf_candidates <- function(n, g, p) {
  X <- matrix(stats::rnorm(n * p, sd = 1.5), n, p)
  wide <- list(w = g$w, m = g$m, L = lapply(g$L, function(L) 1.5 * L))
  rbind(X, g$m, .cf_draw(n, wide))
}

## Best mean for a new component N(m, s_probe^2 I). The score of a mean is
## the rate at which the loss falls when mass moves to that component; it
## uses only target values already held. Returns the best of the rows of
## `X` after a short refinement.
.cf_place <- function(X, W, ph, g, s_probe) {
  rho <- ph - .cf_mixture(W, g)
  q <- Conj(rho) * exp(-0.5 * s_probe^2 * rowSums(W^2)) / nrow(W)
  h <- as.vector(Re(exp(1i * X %*% t(W)) %*% q))
  top <- order(h, decreasing = TRUE)[seq_len(min(3L, nrow(X)))]
  best <- NULL
  for (i in top) {
    op <- stats::optim(
      pmin(pmax(X[i, ], -6), 6),
      function(m) -sum(Re(q * exp(1i * as.vector(W %*% m)))),
      function(m) as.vector(crossprod(W, Im(q * exp(1i * as.vector(W %*% m))))),
      method = "L-BFGS-B", lower = -6, upper = 6, control = list(maxit = 60L)
    )
    if (is.null(best) || op$value < best$value) best <- op
  }
  best$par
}

## Add one component to `g` at the best mean and return the packed start.
.cf_add_component <- function(g, U, ph, p, s_probe) {
  K <- length(g$w) + 1L
  m <- rbind(g$m, .cf_place(.cf_candidates(256L, g, p), U, ph, g, s_probe))
  L <- c(g$L, list(diag(s_probe, p)))
  w <- .cf_ls_weights(
    .cf_components(U, list(w = rep(1, K), m = m, L = L)), ph
  )
  .cf_pack(w, m, L)
}

## Starting mixture: begin with one standard normal component and add one
## component at a time, refitting all of them after each addition. Returns
## the mixture and the loss after each addition.
.cf_grow <- function(K, p, U, ph, lad, bandwidths, maxit_grow) {
  S <- length(bandwidths)
  g <- list(w = 1, m = matrix(0, 1L, p), L = list(diag(p)))
  path <- mean(Mod(ph - .cf_mixture(U, g))^2)
  for (k in seq_len(K)[-1L]) {
    s_probe <- 0.5 * stats::median(unlist(lapply(g$L, diag)))
    par <- .cf_add_component(g, U, ph, p, s_probe)
    fit <- .cf_local(par, U, ph, lad, bandwidths, k, p, S, maxit_grow, 0L)
    g <- .cf_unpack(fit$par, k, p)
    path <- c(path, fit$value)
  }
  list(g = g, path = path)
}

## Relocation: remove the component whose removal costs least, add it back
## where the fit lacks mass, refit, and keep the move when the loss falls
## by at least a tenth. Repeats up to `max_moves` times.
.cf_relocate <- function(fit, U, ph, lad, bandwidths, K, p, maxit_stage,
                         maxit_lm, max_moves) {
  S <- length(bandwidths)
  moves <- 0L
  while (moves < max_moves) {
    g <- .cf_unpack(fit$par, K, p)
    C <- .cf_components(U, g)
    drop_loss <- vapply(seq_len(K), function(k) {
      w <- .cf_ls_weights(C[, -k, drop = FALSE], ph)
      mean(Mod(ph - as.vector(C[, -k, drop = FALSE] %*% w))^2)
    }, numeric(1L))
    improved <- FALSE
    for (k in order(drop_loss)[seq_len(min(2L, K))]) {
      g_rest <- list(
        w = .cf_ls_weights(C[, -k, drop = FALSE], ph),
        m = g$m[-k, , drop = FALSE],
        L = g$L[-k]
      )
      s_probe <- stats::median(unlist(lapply(g_rest$L, diag)))
      par <- .cf_add_component(g_rest, U, ph, p, s_probe)
      cand <- .cf_local(
        par, U, ph, lad, bandwidths, K, p, S, maxit_stage, maxit_lm
      )
      if (cand$value < 0.9 * fit$value) {
        fit <- cand
        improved <- TRUE
        break
      }
    }
    if (!improved) break
    moves <- moves + 1L
  }
  fit$moves <- moves
  fit
}

# ---------------------------------------------------------------------------
# internal helpers: what the characteristic function says about the target
# ---------------------------------------------------------------------------

## Wrap the user's function so that every call is checked for shape and
## finiteness and returns a complex vector.
.cf_checked <- function(cf, p) {
  function(Tm) {
    v <- cf(Tm)
    if (!(is.complex(v) || is.numeric(v)) || length(v) != nrow(Tm)) {
      cli::cli_abort(
        c(
          "{.arg cf} must return one complex (or numeric) value for each row of its argument.",
          "x" = "It returned {.obj_type_friendly {v}} of length {length(v)} for {nrow(Tm)} row{?s}."
        ),
        class = "proxymix_invalid_cf"
      )
    }
    if (any(!is.finite(v))) {
      cli::cli_abort(
        "{.arg cf} returned a value that is not finite.",
        class = "proxymix_invalid_cf"
      )
    }
    as.complex(v)
  }
}

## Frequency at which the modulus first falls to one half along each axis.
.cf_axis_scales <- function(phi, p) {
  tt <- 10^seq(-4, 4, by = 0.01)
  th <- vapply(seq_len(p), function(d) {
    Tm <- matrix(0, length(tt), p)
    Tm[, d] <- tt
    tt[which(Mod(phi(Tm)) < 0.5)[1L]]
  }, numeric(1L))
  if (anyNA(th)) {
    cli::cli_abort(
      c(
        "The modulus of {.arg cf} does not fall to one half along coordinate {which(is.na(th))[1L]}.",
        "i" = "Frequencies from 1e-4 to 1e4 were tried. A distribution on a lattice, or one with a point mass, is outside what {.fn from_cf} fits."
      ),
      class = "proxymix_invalid_cf"
    )
  }
  th
}

## Tail index along each axis: the power at which -log|cf| grows near
## zero. It is 2 when the variance is finite and alpha for a stable law.
.cf_tail_index <- function(phi, p, th) {
  vapply(seq_len(p), function(d) {
    Tm <- matrix(0, 2L, p)
    Tm[, d] <- c(0.01, 0.1) * th[d]
    v <- -log(Mod(phi(Tm)))
    log(v[2L] / v[1L]) / log(10)
  }, numeric(1L))
}

## Phase of the characteristic function at `t` along axis `d`, unwrapped by
## doubling from a small frequency: psi(2s) = 2 psi(s) + Arg(cf(2s) / cf(s)^2).
.cf_unwrapped_phase <- function(phi, p, d, t, n_double = 20L) {
  at <- function(s) {
    v <- numeric(p)
    v[d] <- s
    phi(matrix(v, 1L, p))
  }
  s <- t / 2^n_double
  psi <- Arg(at(s))
  for (k in seq_len(n_double)) {
    psi <- 2 * psi + Arg(at(2 * s) / at(s)^2)
    s <- 2 * s
  }
  psi
}

## Principal phase at `t` along axis `d`, shifted by the multiple of 2 pi
## that matches the continued phase. Unchanged when no wrap occurred.
.cf_phase_at <- function(phi, p, d, t) {
  v <- numeric(p)
  v[d] <- t
  principal <- Arg(phi(matrix(v, 1L, p)))
  turns <- round((.cf_unwrapped_phase(phi, p, d, t) - principal) / (2 * pi))
  principal + 2 * pi * turns
}

## Mean and covariance by differencing the log characteristic function at
## zero. Both must exist.
.cf_moments <- function(phi, p, th) {
  h <- 0.02 * th
  lp <- function(v) log(phi(matrix(v, 1L, p)))
  e <- function(d) {
    v <- numeric(p)
    v[d] <- h[d]
    v
  }
  mu <- vapply(seq_len(p), function(d) {
    m <- Im(lp(e(d)) - lp(-e(d))) / (2 * h[d])
    psi <- .cf_phase_at(phi, p, d, h[d])
    ## the principal logarithm wraps once mu * h passes pi
    if (abs(psi - h[d] * m) > pi) m <- psi / h[d]
    m
  }, numeric(1L))
  S <- matrix(0, p, p)
  for (d in seq_len(p)) {
    S[d, d] <- -Re(lp(e(d)) + lp(-e(d))) / h[d]^2
    for (f in seq_len(d - 1L)) {
      S[d, f] <- -Re(
        lp(e(d) + e(f)) - lp(e(d) - e(f)) - lp(e(f) - e(d)) + lp(-e(d) - e(f))
      ) / (4 * h[d] * h[f])
      S[f, d] <- S[d, f]
    }
  }
  list(mean = mu, cov = S)
}

## Checks that need the frequency scale: value 1 at zero, modulus at most
## 1, and cf(-t) equal to the conjugate of cf(t).
.check_cf_properties <- function(phi, p, th) {
  at_zero <- phi(matrix(0, 1L, p))
  if (Mod(at_zero - 1) > 1e-8) {
    cli::cli_abort(
      "{.arg cf} must equal 1 at zero; it returned {format(at_zero)}.",
      class = "proxymix_invalid_cf"
    )
  }
  Tm <- matrix(stats::rnorm(64L * p), 64L, p) %*% diag(th, p)
  v <- phi(Tm)
  if (any(Mod(v) > 1 + 1e-8)) {
    cli::cli_abort(
      "{.arg cf} has modulus above 1 (largest {signif(max(Mod(v)), 6)}), so it is not a characteristic function.",
      class = "proxymix_invalid_cf"
    )
  }
  gap <- max(Mod(phi(-Tm) - Conj(v)))
  if (gap > 1e-8) {
    cli::cli_abort(
      "{.arg cf} at {.code -t} must be the complex conjugate of {.arg cf} at {.code t}; the largest gap found is {signif(gap, 3)}.",
      class = "proxymix_invalid_cf"
    )
  }
  invisible(NULL)
}

## Coordinates for the fit: x = mu + A y. With a finite variance, `mu` and
## `A A'` are the target's mean and covariance; otherwise `mu` comes from
## the phase of the characteristic function and `A` is diagonal, set by
## the half-height frequency of each axis.
.cf_standardise <- function(phi, p, th, heavy) {
  if (!heavy) {
    mo <- .cf_moments(phi, p, th)
    ev <- eigen(mo$cov, symmetric = TRUE, only.values = TRUE)$values
    if (all(is.finite(ev)) && min(ev) > 1e-10 * max(ev)) {
      return(list(
        mu = mo$mean, A = t(chol(mo$cov)), mean = mo$mean, cov = mo$cov
      ))
    }
  }
  mu <- vapply(seq_len(p), function(d) {
    .cf_phase_at(phi, p, d, th[d]) / th[d]
  }, numeric(1L))
  list(mu = mu, A = diag(.cf_half_height / th, p), mean = NULL, cov = NULL)
}

## Frequencies in standardised coordinates: Gaussian, cycling through the
## bandwidths. Returns the matrix and the bandwidth index of each row.
.cf_frequencies <- function(M, p, bandwidths) {
  lad <- rep_len(seq_along(bandwidths), M)
  U <- matrix(stats::rnorm(M * p), M, p) *
    (.cf_half_height * bandwidths[lad] / sqrt(p))
  list(U = U, lad = lad)
}

.check_cf_arguments <- function(cf, n_dim, N, n_freq, bandwidths, n_start,
                                max_iter) {
  if (!is.function(cf)) {
    cli::cli_abort("{.arg cf} must be a function of a matrix of frequencies.")
  }
  count <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 1 && x == round(x)
  }
  if (!count(n_dim)) {
    cli::cli_abort("{.arg n_dim} must be a positive integer scalar.")
  }
  if (!count(N)) {
    cli::cli_abort("{.arg N} must be a positive integer scalar.")
  }
  if (!count(n_start)) {
    cli::cli_abort("{.arg n_start} must be a positive integer scalar.")
  }
  if (!count(max_iter)) {
    cli::cli_abort("{.arg max_iter} must be a positive integer scalar.")
  }
  if (!is.numeric(bandwidths) || length(bandwidths) < 1L ||
      anyNA(bandwidths) || any(bandwidths <= 0) || is.unsorted(bandwidths)) {
    cli::cli_abort(
      "{.arg bandwidths} must be positive numbers in increasing order."
    )
  }
  if (!count(n_freq) || n_freq < 64L * length(bandwidths)) {
    cli::cli_abort(
      "{.arg n_freq} must be an integer of at least 64 per bandwidth ({64L * length(bandwidths)})."
    )
  }
  invisible(NULL)
}

## `init` must be a mixture with `K` components in `p` dimensions.
.check_cf_init <- function(init, p, K) {
  if (is.null(init)) return(invisible(NULL))
  if (!S7::S7_inherits(init, gmm)) {
    cli::cli_abort("{.arg init} must be a {.cls gmm} object or {.code NULL}.")
  }
  if (gmm_dim(init) != p || length(init@weights) != K) {
    cli::cli_abort(
      c(
        "{.arg init} must have {.arg N} components in {.arg n_dim} dimensions.",
        "x" = "It has {length(init@weights)} component{?s} in {gmm_dim(init)} dimension{?s}; {.arg N} is {K} and {.arg n_dim} is {p}."
      )
    )
  }
  invisible(NULL)
}

## Convert `init` from x to the standardised coordinates y, where x = mu + A y.
.cf_init_standardised <- function(init, co) {
  A_inv <- solve(co$A)
  m <- t(vapply(
    init@means, function(mk) as.vector(A_inv %*% (mk - co$mu)),
    numeric(length(co$mu))
  ))
  if (length(co$mu) == 1L) m <- matrix(m, ncol = 1L)
  L <- lapply(init@covariances, function(S) {
    t(chol(A_inv %*% S %*% t(A_inv)))
  })
  list(w = init@weights, m = m, L = L)
}

# ---------------------------------------------------------------------------
# from_cf
# ---------------------------------------------------------------------------

#' Fit a Gaussian mixture to a characteristic function
#'
#' Fits a Gaussian mixture to a distribution that is known through its
#' characteristic function \eqn{\varphi(t) = E \exp(i t^\top X)}. No density
#' and no sample of the distribution are needed. The result is an ordinary
#' fitted mixture, so [rgmm()], [dgmm()], [gmm_affine()], [gmm_convolve()]
#' and the other mixture operations apply to it.
#'
#' The fit minimises the mean of
#' \eqn{|\varphi(t_j) - \varphi_g(t_j)|^2} over a fixed sample of
#' frequencies \eqn{t_j}, where \eqn{\varphi_g} is the characteristic
#' function of the mixture ([gmm_cf()]). The frequencies are Gaussian with
#' several bandwidths, scaled to the target so that `bandwidths = 1` is the
#' frequency at which \eqn{|\varphi|} has fallen to about one half. The
#' same criterion is used by Brück (2025) to train a neural sampler.
#'
#' The starting mixture is built one component at a time. Each new
#' component is placed at the mean where the fit falls furthest short of the
#' target, after both are smoothed by a Gaussian. That mean is computed from
#' the characteristic-function residual, so no density of the target is
#' needed. All components are refitted after each addition. The fit then tries
#' moving each of the two least useful components to a better mean, keeping
#' a move if it cuts the loss by at least 10%, up to `N` moves.
#'
#' @section Checking a fit:
#' The returned object reports the root mean square and the maximum of
#' \eqn{|\varphi - \varphi_g|} on `n_freq` fresh frequencies that the fit
#' did not use (`heldout_rms` and `heldout_max` in `diagnostics`).
#'
#' With `bracket = TRUE` (the default) the fit also carries a lower and an
#' upper bound for the Kolmogorov distance between the target and the fit,
#' computed by [gmm_cf_distance()] along eight directions: the first
#' `min(n_dim, 8)` coordinate axes and, when `n_dim < 8`, random unit
#' directions to make up eight (one direction when `n_dim = 1`). The bounds
#' for each direction are in `diagnostics$bracket`, the directions in
#' `diagnostics$bracket_directions`, and the largest of each over the
#' directions in the quality record as `cf_distance` and `cf_upper`
#' ([gmm_fit_quality()]). These bounds compare distribution functions; they
#' are not a Kullback-Leibler divergence. Drawing the random directions does
#' not change the state of the random number generator; with `seed` they
#' follow the seed. Set `bracket = FALSE` to skip the bounds.
#'
#' @section Heavy tails:
#' A Gaussian mixture has finite moments of every order, so it cannot
#' reproduce infinite-variance tails and fits the centre only. The tail
#' index is read from \eqn{-\log|\varphi|} near zero. It is 2 when the
#' variance is finite and equals the stability index for a stable law. When
#' it is below 1.9 a warning of class `proxymix_cf_heavy_tail` is raised.
#'
#' @section Scope:
#' The target needs a characteristic function that decays, which excludes
#' lattice distributions and distributions with a point mass.
#'
#' @param cf A function that takes an `n` by `n_dim` matrix of frequencies,
#'   one per row, and returns the characteristic function at each row as a
#'   complex vector of length `n`. A real vector is accepted for a symmetric
#'   target.
#' @param n_dim Dimension of the target.
#' @param N Number of mixture components.
#' @param n_freq Number of frequencies in the fit, at least 64 times
#'   `length(bandwidths)`. Default `8192L`.
#' @param bandwidths Increasing positive numbers, the typical length of a
#'   frequency vector in units of the target's half-height frequency (each
#'   coordinate has standard deviation `bandwidth / sqrt(n_dim)` of that).
#'   The grown start and `init` use all bandwidths together. Random starts
#'   are fitted from the lowest bandwidth upwards. Extend the range upwards
#'   for a target with fine structure and downwards for a heavy-tailed
#'   target.
#' @param n_start Number of starts. The first start is `init` if given,
#'   otherwise the grown start described above; the others are random. The
#'   fit with the lowest loss is returned. Default `1L`. With ten or more
#'   components in five or more dimensions a single start can stop at a
#'   poor fit, so four starts are a safer choice. For ten components in ten
#'   dimensions the grown start stopped at a poor fit and a random start
#'   recovered the target.
#' @param max_iter Maximum number of L-BFGS iterations per frequency stage.
#'   A Levenberg-Marquardt step with at most 200 iterations follows when the
#'   fit has at most 700 parameters. Default `300L`.
#' @param init Optional [gmm] with `N` components to start the first fit
#'   from, in place of the grown start. Default `NULL`.
#' @param seed Optional integer seed. When supplied, the random number
#'   generator state is restored on exit.
#' @param bracket Logical. Compute the Kolmogorov bounds described under
#'   *Checking a fit*. Default `TRUE`.
#'
#' @returns A [gmm_fit] with `regime` set to `"cf"`. Its `diagnostics` hold
#'   `loss` (the final loss), `loss_by_start`, `loss_path` (the loss as the
#'   grown start went from 1 to `N` components, or `NULL` when `init` is
#'   given), `heldout_rms` and
#'   `heldout_max`, `relocations` (the number of moves kept), `tail_index`
#'   (one value for each coordinate), `cf_mean` and `cf_cov` (the target's
#'   mean and covariance read from `cf`, or `NULL` when the target has no
#'   finite, positive definite covariance), `n_freq`, `bandwidths`, and,
#'   with `bracket = TRUE`, `bracket` and `bracket_directions`.
#'   `converged` is `FALSE` when the final optimisation of the returned fit
#'   stopped at its iteration limit.
#' @references
#' Brück, F. (2025). Generative neural networks for characteristic
#' functions. *Journal of Computational and Graphical Statistics*, 34(4),
#' 1456--1465. \doi{10.1080/10618600.2025.2455135}
#' @family fitting
#' @seealso [cf_target()] and [fit_proxymix()] for the same fit through the
#'   unified verb; [from_cf_elliptical()] and [from_cf_sources()] for
#'   structured targets; [gmm_cf()], [gmm_cf_distance()].
#' @export
#' @examples
#' \dontshow{if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
#'   blas_threads <- RhpcBLASctl::blas_get_num_procs()
#'   RhpcBLASctl::blas_set_num_threads(1L)
#' }}
#' ## The sum of a Gamma(3) variable and independent N(0, 0.25) noise.
#' cf <- function(t) (1 - 1i * t[, 1])^(-3) * exp(-0.125 * t[, 1]^2)
#' fit <- from_cf(cf, n_dim = 1L, N = 3L, n_freq = 1024L, seed = 1L)
#' ## the print shows the held-out error and the Kolmogorov bounds
#' fit
#' fit@diagnostics$bracket
#' \dontshow{if (exists("blas_threads")) {
#'   RhpcBLASctl::blas_set_num_threads(blas_threads)
#'   rm(blas_threads)
#' }}
from_cf <- function(
  cf,
  n_dim,
  N,
  n_freq = 8192L,
  bandwidths = c(0.3, 0.6, 1.2, 2.5, 5),
  n_start = 1L,
  max_iter = 300L,
  init = NULL,
  seed = NULL,
  bracket = TRUE
) {
  .check_cf_arguments(cf, n_dim, N, n_freq, bandwidths, n_start, max_iter)
  .check_cf_bracket_argument(bracket)
  p <- as.integer(n_dim)
  K <- as.integer(N)
  .check_cf_init(init, p, K)
  n_freq <- as.integer(n_freq)
  n_start <- as.integer(n_start)
  max_iter <- as.integer(max_iter)
  phi <- .cf_checked(cf, p)
  S <- length(bandwidths)

  run <- function() {
    # Scale, checks and coordinates -----------------------------------------
    th <- .cf_axis_scales(phi, p)
    .check_cf_properties(phi, p, th)
    tail_index <- .cf_tail_index(phi, p, th)
    heavy <- min(tail_index) < 1.9
    co <- .cf_standardise(phi, p, th, heavy)
    A_inv <- solve(co$A)
    phi_std <- function(U) {
      Tm <- U %*% A_inv
      exp(-1i * as.vector(Tm %*% co$mu)) * phi(Tm)
    }
    fr <- .cf_frequencies(n_freq, p, bandwidths)
    ph <- phi_std(fr$U)

    # Starts ------------------------------------------------------------------
    best <- NULL
    loss_by_start <- numeric(n_start)
    loss_path <- NULL
    for (i in seq_len(n_start)) {
      if (i == 1L && !is.null(init)) {
        g0 <- .cf_init_standardised(init, co)
        stages <- S
      } else if (i == 1L) {
        grown <- .cf_grow(K, p, fr$U, ph, fr$lad, bandwidths, max_iter %/% 2L)
        loss_path <- grown$path
        g0 <- grown$g
        stages <- S
      } else {
        g0 <- list(
          w = rep(1 / K, K),
          m = matrix(stats::rnorm(K * p), K, p),
          L = rep(list(diag(K^(-1 / p), p)), K)
        )
        stages <- seq_len(S)
      }
      fit <- .cf_local(
        .cf_pack(g0$w, g0$m, g0$L), fr$U, ph, fr$lad, bandwidths, K, p,
        stages, max_iter, 200L
      )
      fit$moves <- 0L
      if (fit$value > 1e-14 && K > 1L) {
        fit <- .cf_relocate(
          fit, fr$U, ph, fr$lad, bandwidths, K, p, max_iter, 10L, K
        )
        if (fit$moves > 0L) {
          moves <- fit$moves
          fit <- .cf_local(
            fit$par, fr$U, ph, fr$lad, bandwidths, K, p, S, max_iter, 200L
          )
          fit$moves <- moves
        }
      }
      loss_by_start[i] <- fit$value
      if (is.null(best) || fit$value < best$value) best <- fit
    }

    # Error on frequencies the fit did not use --------------------------------
    gs <- .cf_unpack(best$par, K, p)
    fresh <- .cf_frequencies(n_freq, p, bandwidths)
    gap <- Mod(phi_std(fresh$U) - .cf_mixture(fresh$U, gs))
    list(
      gs = gs, co = co, best = best, tail_index = tail_index, heavy = heavy,
      loss_by_start = loss_by_start, loss_path = loss_path,
      heldout_rms = sqrt(mean(gap^2)), heldout_max = max(gap)
    )
  }
  res <- if (is.null(seed)) run() else withr::with_seed(seed, run())

  if (res$heavy) {
    cli::cli_warn(
      c(
        "The target has heavy tails: the tail index read from {.arg cf} is {signif(min(res$tail_index), 3)}, below 2.",
        "i" = "The mixture has finite variance. It can match the centre of the distribution and not its far tails."
      ),
      class = "proxymix_cf_heavy_tail"
    )
  }

  # Back to the original coordinates ------------------------------------------
  A <- res$co$A
  means <- lapply(seq_len(K), function(k) {
    as.numeric(A %*% res$gs$m[k, ] + res$co$mu)
  })
  covs <- lapply(res$gs$L, function(L) tcrossprod(A %*% L))
  converged <- !isTRUE(res$best$hit_limit)
  kb <- .cf_bracket(
    gmm(weights = res$gs$w, means = means, covariances = covs), phi, p,
    bracket, seed
  )

  fit <- gmm_fit(
    weights = res$gs$w,
    means = means,
    covariances = covs,
    target = NULL,
    regime = "cf",
    diagnostics = list(
      loss = res$best$value,
      loss_by_start = res$loss_by_start,
      loss_path = res$loss_path,
      heldout_rms = res$heldout_rms,
      heldout_max = res$heldout_max,
      relocations = res$best$moves,
      tail_index = res$tail_index,
      cf_mean = res$co$mean,
      cf_cov = res$co$cov,
      n_freq = n_freq,
      bandwidths = bandwidths,
      bracket = kb$table,
      bracket_directions = kb$directions
    ),
    converged = converged,
    iterations = NA_integer_,
    call = match.call(),
    name = sprintf("from_cf[N=%d]", K),
    metadata = list(quality = list(
      regime = "cf",
      converged = converged,
      degenerate = FALSE,
      ess = NA_real_,
      ess_relative = NA_real_,
      min_component_ess = NA_real_,
      max_weight = NA_real_,
      support_fraction = NA_real_,
      kld_final = NA_real_,
      kld_approx = NA_real_,
      validation_gap = NA_real_,
      heldout_kld = NA_real_,
      cf_heldout_rms = res$heldout_rms,
      cf_distance = kb$distance,
      cf_upper = kb$upper
    ))
  )
  gmm_canonicalise(fit)
}

.check_cf_bracket_argument <- function(bracket) {
  if (!is.logical(bracket) || length(bracket) != 1L || is.na(bracket)) {
    cli::cli_abort("{.arg bracket} must be {.code TRUE} or {.code FALSE}.")
  }
  invisible(NULL)
}

## Directions for the Kolmogorov bounds of a fit.
.cf_bracket_directions <- function(p, seed = NULL) {
  axes <- diag(p)[seq_len(min(p, 8L)), , drop = FALSE]
  if (p == 1L || p >= 8L) return(axes)
  draw <- function() .cf_unit_directions(8L - p, p)
  rbind(axes, if (is.null(seed)) withr::with_preserve_seed(draw()) else
    withr::with_seed(seed, draw()))
}

## Kolmogorov bounds of `g` against the characteristic function `phi` along
## each bracket direction. NULL table and NA summaries when `bracket` is FALSE.
.cf_bracket <- function(g, phi, p, bracket, seed = NULL) {
  if (!bracket) {
    return(list(table = NULL, directions = NULL,
                distance = NA_real_, upper = NA_real_))
  }
  D <- .cf_bracket_directions(p, seed)
  rows <- lapply(seq_len(nrow(D)), function(i) {
    d <- gmm_cf_distance(g, phi, direction = D[i, ])
    data.frame(
      direction = i, distance = d$distance, upper = d$upper,
      upper_inside = d$upper_inside, truncation = d$truncation,
      cutoff_gap = d$cutoff_gap
    )
  })
  tab <- do.call(rbind, rows)
  list(table = tab, directions = D,
       distance = max(tab$distance), upper = max(tab$upper))
}

# ---------------------------------------------------------------------------
# gmm_cf
# ---------------------------------------------------------------------------

#' Characteristic function of a Gaussian mixture
#'
#' Evaluates
#' \eqn{\varphi_g(t) = \sum_k w_k \exp(i t^\top \mu_k -
#' t^\top \Sigma_k t / 2)}
#' at each row of `t`. Here `cf` stands for characteristic function; in
#' [gmm_cf_mean()], [gmm_cf_variance()] and [gmm_cf_tail_prob()] it stands
#' for counterfactual.
#'
#' @param g A [gmm] or [gmm_fit].
#' @param t Numeric matrix of frequencies with one frequency per row and as
#'   many columns as the mixture has dimensions. A vector is read as one
#'   frequency per element for a one-dimensional mixture and as a single
#'   frequency otherwise.
#'
#' @returns A complex vector with one value for each row of `t`.
#' @family ops
#' @seealso [from_cf()], [from_cf_elliptical()], [from_cf_sources()],
#'   [gmm_cf_distance()].
#' @export
#' @examples
#' g <- gmm(
#'   weights = c(0.3, 0.7),
#'   means = list(c(0, 0), c(2, 1)),
#'   covariances = list(diag(2), 0.5 * diag(2))
#' )
#' gmm_cf(g, rbind(c(0, 0), c(0.5, -1)))
gmm_cf <- function(g, t) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("{.arg g} must be a {.cls gmm} object.")
  }
  p <- gmm_dim(g)
  if (!is.numeric(t)) {
    cli::cli_abort("{.arg t} must be numeric.")
  }
  if (is.null(dim(t))) {
    t <- if (p == 1L) matrix(t, ncol = 1L) else matrix(t, nrow = 1L)
  }
  if (ncol(t) != p) {
    cli::cli_abort(
      "{.arg t} must have {p} column{?s}, one for each dimension of {.arg g}; it has {ncol(t)}."
    )
  }
  out <- complex(nrow(t))
  for (k in seq_along(g@weights)) {
    quad <- rowSums((t %*% g@covariances[[k]]) * t)
    out <- out + g@weights[k] *
      exp(1i * as.vector(t %*% g@means[[k]]) - 0.5 * quad)
  }
  out
}

# ---------------------------------------------------------------------------
# gmm_cf_distance
# ---------------------------------------------------------------------------

#' Kolmogorov distance between a characteristic-function target and a mixture
#'
#' For a direction \eqn{a}, compares the distribution of \eqn{a^\top X}
#' under the target with its distribution under the mixture and returns a
#' lower and an upper bound for the largest absolute difference between the
#' two distribution functions (the Kolmogorov distance). The differences are
#' obtained by numerical inversion of the two characteristic functions
#' (Gil-Pelaez, 1951), so neither distribution is sampled.
#'
#' If `upper` is at most \eqn{\varepsilon}, every quantile of the target's
#' projection at probability \eqn{q} lies between the mixture's quantiles at
#' probabilities \eqn{q - \varepsilon} and \eqn{q + \varepsilon}.
#'
#' `distance` is the largest difference at `n_grid` points, spaced equally
#' in mixture probability between the two values in `probs`. It is a lower
#' bound.
#'
#' `upper` is an upper bound, valid when the characteristic function of the
#' projected target is integrable and `truncation` and `cutoff_gap` are
#' small. At a point between two grid
#' points, the difference exceeds the larger end value by at most the
#' probability step, and exceeds the nearer end value by at most half the
#' spacing times a bound on the largest density difference. Outside the grid it adds the tail
#' probability `probs[1]` or `1 - probs[2]` to the difference at the nearest
#' end. Move `probs` towards 0 and 1 to shrink the outside term.
#'
#' Both bounds are computed from the numerically integrated differences, so
#' they hold up to the error of that integration, which neither includes.
#' Against closed-form distribution functions (Cauchy, normal, and a Gamma
#' plus uniform sum, by quadrature) that error was below `5e-6`, and
#' `distance` exceeded the exact Kolmogorov distance by at most `3e-6`.
#'
#' @param g A [gmm] or [gmm_fit].
#' @param cf The target's characteristic function, as in [from_cf()].
#' @param direction Numeric vector with one entry for each dimension of
#'   `g`, not all zero. Default: the first coordinate axis.
#' @param n_grid Number of grid points. Default `201L`.
#' @param probs The lowest and highest probability of the grid. Default
#'   `c(0.005, 0.995)`.
#'
#' @returns A list with `distance` (the largest difference on the grid),
#'   `at` (the grid point where it is attained), `upper` (the upper bound
#'   described above), `upper_inside` (the bound between the two quantiles
#'   in `probs`), `truncation` (a bound on the omitted part of the inversion
#'   integral from the highest frequency integrated to 20 times that
#'   frequency, added to `upper` and `upper_inside`; the part beyond is
#'   ignored), `x` (the grid), `difference` (target minus
#'   mixture distribution function at each grid point), and `cutoff_gap`
#'   (the modulus of the difference between the two characteristic
#'   functions at the highest frequency integrated). A `cutoff_gap` or
#'   `truncation` that is not small means the integral was cut off early,
#'   which happens when the target has a point mass or a very slowly
#'   decaying characteristic function. The bounds do not hold then.
#' @references
#' Gil-Pelaez, J. (1951). Note on the inversion theorem. *Biometrika*,
#' 38(3--4), 481--482. \doi{10.1093/biomet/38.3-4.481}
#' @family diagnostics
#' @seealso [from_cf()], [from_cf_elliptical()], [from_cf_sources()],
#'   [gmm_cf()].
#' @export
#' @examples
#' ## A normal target against a normal with a shifted mean.
#' cf <- function(t) exp(-0.5 * t[, 1]^2)
#' g <- gmm(weights = 1, means = list(0.1), covariances = list(matrix(1)))
#' gmm_cf_distance(g, cf)$distance
#' 2 * pnorm(0.05) - 1
gmm_cf_distance <- function(
  g,
  cf,
  direction = NULL,
  n_grid = 201L,
  probs = c(0.005, 0.995)
) {
  if (!S7::S7_inherits(g, gmm)) {
    cli::cli_abort("{.arg g} must be a {.cls gmm} object.")
  }
  if (!is.function(cf)) {
    cli::cli_abort("{.arg cf} must be a function of a matrix of frequencies.")
  }
  p <- gmm_dim(g)
  if (is.null(direction)) direction <- c(1, rep(0, p - 1L))
  if (!is.numeric(direction) || length(direction) != p ||
      any(!is.finite(direction)) || all(direction == 0)) {
    cli::cli_abort(
      "{.arg direction} must be {p} finite number{?s}, not all zero."
    )
  }
  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs) ||
      probs[1L] <= 0 || probs[2L] >= 1 || probs[1L] >= probs[2L]) {
    cli::cli_abort(
      "{.arg probs} must be two increasing probabilities strictly inside (0, 1)."
    )
  }
  n_grid <- as.integer(n_grid)
  if (length(n_grid) != 1L || is.na(n_grid) || n_grid < 2L) {
    cli::cli_abort("{.arg n_grid} must be an integer of at least 2.")
  }
  phi <- .cf_checked(cf, p)
  a <- as.numeric(direction)

  # The mixture's projection ---------------------------------------------------
  w <- g@weights
  mk <- vapply(g@means, function(m) sum(m * a), numeric(1L))
  sk <- sqrt(vapply(
    g@covariances, function(S) sum(a * (S %*% a)), numeric(1L)
  ))
  pmix <- function(x) {
    vapply(x, function(z) sum(w * stats::pnorm(z, mk, sk)), numeric(1L))
  }
  qmix <- function(q) {
    stats::uniroot(
      function(x) pmix(x) - q, range(mk) + c(-1, 1) * 40 * max(sk),
      tol = 1e-10
    )$root
  }
  qq <- vapply(c(0.25, 0.5, 0.75), qmix, numeric(1L))
  centre <- qq[2L]
  sc <- (qq[3L] - qq[1L]) / 1.349
  ## grid at equal steps of the mixture's distribution function
  pr <- seq(probs[1L], probs[2L], length.out = n_grid)
  x <- vapply(pr, qmix, numeric(1L)) - centre

  # Difference of the two characteristic functions, centred -------------------
  del <- function(s) {
    phi(outer(s, a)) * exp(-1i * s * centre) -
      as.vector(exp(1i * outer(s, mk - centre) - 0.5 * outer(s^2, sk^2)) %*% w)
  }
  ## integrate up to where the difference has died out
  s_max <- 4 / sc
  while (s_max < 400 / sc &&
         max(Mod(del(s_max * c(1, 1.1, 1.2, 1.35, 1.5)))) > 1e-7) {
    s_max <- 1.5 * s_max
  }

  # Gil-Pelaez integral: log-spaced near zero, then the midpoint rule ---------
  s0 <- 0.1 / sc
  v <- seq(log(1e-8 * s0), log(s0), length.out = 400L)
  dv <- v[2L] - v[1L]
  s_lo <- exp(v)
  ds <- min(pi / (8 * max(abs(x))), s0 / 8)
  s_hi <- seq(s0 + ds / 2, s_max, by = ds)
  d_lo <- del(s_lo)
  d_hi <- del(s_hi)
  difference <- vapply(x, function(z) {
    lo <- Im(exp(-1i * s_lo * z) * d_lo)
    hi <- Im(exp(-1i * s_hi * z) * d_hi) / s_hi
    -(sum(lo) * dv - 0.5 * (lo[1L] + lo[400L]) * dv + sum(hi) * ds) / pi
  }, numeric(1L))

  # Upper bound -----------------------------------------------------------------
  ## frequencies from the cutoff to 20 times the cutoff, log-spaced
  v_t <- seq(log(s_max), log(20 * s_max), length.out = 200L)
  dv_t <- v_t[2L] - v_t[1L]
  d_t <- Mod(del(exp(v_t)))
  truncation <- sum(d_t) * dv_t / pi
  ## bound on the largest difference of the two densities
  density_gap <- (
    sum(Mod(d_lo) * s_lo) * dv + sum(Mod(d_hi)) * ds +
      sum(d_t * exp(v_t)) * dv_t
  ) / pi
  i <- which.max(abs(difference))
  ## bounds between grid points; see ?gmm_cf_distance
  between <- pmin(0.5 * diff(x) * density_gap, diff(pr))
  at_ends <- pmax(abs(difference[-n_grid]), abs(difference[-1L]))
  upper_inside <- max(at_ends + between) + truncation
  upper_outside <- max(
    probs[1L] + abs(difference[1L]),
    1 - probs[2L] + abs(difference[n_grid])
  ) + truncation

  list(
    distance = abs(difference[i]),
    at = x[i] + centre,
    upper = max(upper_inside, upper_outside),
    upper_inside = upper_inside,
    truncation = truncation,
    x = x + centre,
    difference = difference,
    cutoff_gap = Mod(del(s_max))
  )
}
