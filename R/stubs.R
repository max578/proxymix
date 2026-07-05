## Planned features (not yet implemented). INTERNAL as of v0.13.0: the
## placeholders keep stable signatures and classed not-yet-implemented
## errors (guarded by tests/testthat/test-stubs.R via `:::`), but are no
## longer exported -- an exported function whose only behaviour is to
## error reads as vapourware on the public surface.

stub_not_yet_implemented <- function(fn_name) {
  cli::cli_abort(
    c("`{fn_name}()` is not yet implemented.",
      "i" = "It is a planned feature; see the package NEWS for status."),
    class = "proxymix_not_yet_implemented",
    .frame = parent.frame()
  )
}

## from_aggregate_likelihood(): planned -- fit a proxy f_theta(x) whose
## aggregate-likelihood pushforward g(y) = E[x|y][f(x)] stays closed form.
from_aggregate_likelihood <- function(y, latent_aggregator, N = 3L, ...) {
  stub_not_yet_implemented("from_aggregate_likelihood")
}

## fit_kld_em_collider(): planned -- KLD-EM projected onto densities
## satisfying a DAG-implied conditional-independence set.
fit_kld_em_collider <- function(target, dag, N = 3L, ...) {
  stub_not_yet_implemented("fit_kld_em_collider")
}

## to_apsim_scenarios(): planned -- samples from a fit to a scenario table
## for an external simulator, under a user schema.
to_apsim_scenarios <- function(fit, n = 100L, schema = list(), ...) {
  stub_not_yet_implemented("to_apsim_scenarios")
}

## from_simulator(): planned -- probe an expensive simulator over a design,
## expose a KDE / empirical-likelihood surface as an evaluable gmm_target.
from_simulator <- function(simulator, design, ...) {
  stub_not_yet_implemented("from_simulator")
}
