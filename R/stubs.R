## Tier-2 stubs.
##
## Each stub ships with a stable signature, full roxygen documentation, and
## a body that signals a `proxymix_not_yet_implemented` condition pointing
## the user at `vignettes/roadmap.Rmd`. The signature-stability test in
## `tests/testthat/test-stubs.R` exists to prevent silent API drift before
## the bodies graduate to Tier 1.

stub_not_yet_implemented <- function(fn_name) {
  cli::cli_abort(
    c("`{fn_name}()` is not yet implemented.",
      "i" = "See {.file vignettes/roadmap.Rmd} for the planned scope and timing."),
    class = "proxymix_not_yet_implemented",
    .frame = parent.frame()
  )
}

#' From an aggregate likelihood to a Gaussian-mixture proxy (stub)
#'
#' **Tier-2 stub** — signature stable; body not yet implemented.
#'
#' Plan: fit a Gaussian-mixture proxy `f_theta(x)` so that, when used as
#' the latent-level density inside an aggregate-likelihood downscaling
#' framework `g(y) = E_{x | y}[f(x)]`, the resulting downscaling likelihood
#' has tractable closed-form marginalisation. Targets the
#' kernel-downsizing application described by Sejdinovic et al.
#'
#' @param y A numeric matrix or vector of aggregate-level observations.
#' @param latent_aggregator A function mapping latent `x` to aggregate
#'   space.
#' @param N Number of Gaussian components in the proxy.
#' @param ... Future configuration arguments.
#'
#' @returns A [gmm_fit] (when implemented).
#' @family stubs
#' @export
#' @examples
#' try(from_aggregate_likelihood(matrix(0, 1, 1),
#'                               latent_aggregator = identity,
#'                               N = 2L))
from_aggregate_likelihood <- function(y, latent_aggregator, N = 3L, ...) {
  stub_not_yet_implemented("from_aggregate_likelihood")
}

#' KLD-EM under collider-induced conditional-independence constraints (stub)
#'
#' **Tier-2 stub** — signature stable; body not yet implemented.
#'
#' Plan: a regime-(iii) KLD-EM that projects each iteration onto the
#' manifold of joint densities respecting a DAG-implied set of
#' conditional-independence constraints. Targets the
#' collider-regularised regression direction described by Sejdinovic et al.
#'
#' @param target A [gmm_target] with a non-NULL `log_density`.
#' @param dag A description of the DAG structure (planned: an adjacency
#'   matrix or a `dagitty` object).
#' @param N Number of mixture components.
#' @param ... Future configuration arguments.
#'
#' @returns A [gmm_fit] (when implemented).
#' @family stubs
#' @export
#' @examples
#' try(fit_kld_em_collider(banana_target(), dag = matrix(0, 2, 2), N = 2L))
fit_kld_em_collider <- function(target, dag, N = 3L, ...) {
  stub_not_yet_implemented("fit_kld_em_collider")
}

#' Gaussian-mixture samples to APSIM scenario tables (stub)
#'
#' **Tier-2 stub** — signature stable; body not yet implemented.
#'
#' Plan: convert samples from a fitted [gmm_fit] into the tabular format
#' consumed by an APSIM scenario runner, optionally honouring a
#' user-supplied schema mapping mixture coordinates to APSIM variable names
#' and units.
#'
#' @param fit A [gmm_fit].
#' @param n Number of scenarios to generate.
#' @param schema A list mapping mixture coordinates to APSIM variable names
#'   and units.
#' @param ... Future configuration arguments.
#'
#' @returns A data frame in APSIM scenario format (when implemented).
#' @family stubs
#' @export
#' @examples
#' x <- matrix(stats::rnorm(200), ncol = 2)
#' fit <- fit_proxymix(gmm_target_from_samples(x), N = 2L, max_iter = 10L)
#' try(to_apsim_scenarios(fit, n = 100L, schema = list()))
to_apsim_scenarios <- function(fit, n = 100L, schema = list(), ...) {
  stub_not_yet_implemented("to_apsim_scenarios")
}

#' Wrap an expensive simulator as a `gmm_target` (stub)
#'
#' **Tier-2 stub** — signature stable; body not yet implemented.
#'
#' Plan: probe an expensive simulator over a Latin-hypercube design,
#' build a kernel-density estimate (or empirical-likelihood surface) on
#' the simulator outputs, and expose the result as a [gmm_target] with
#' an evaluable `log_density`.
#'
#' @param simulator A function `function(x)` mapping inputs to outputs.
#' @param design An `n` by `p` matrix of input design points.
#' @param ... Future configuration arguments.
#'
#' @returns A [gmm_target] (when implemented).
#' @family stubs
#' @export
#' @examples
#' try(from_simulator(simulator = identity,
#'                    design = matrix(stats::rnorm(20), ncol = 2)))
from_simulator <- function(simulator, design, ...) {
  stub_not_yet_implemented("from_simulator")
}
