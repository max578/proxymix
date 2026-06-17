# proxymix decision-layer benchmarks -- results

Reproduce with `inst/benchmarks/uplift_benchmarks.R` (set `PROXYMIX_BENCH_DATA`
to the folder holding `hillstrom.csv` and `criteo-uplift-v2.1.csv.gz`).
Run 2026-06-17. Comparators: `grf` 2.6.1 causal forest, S/T/X meta-learners
(grf regression forests as base learners), `DoubleML` 1.0.2 IRM for the ATE.
These are independent implementations (Independent Oracle Principle).

## Synthetic PEHE -- proxymix vs grf, known ground truth (lower is better)

| DGP | proxymix | grf | ratio |
|---|---|---|---|
| linear | 0.0355 | 0.0834 | 0.43 |
| crossing | 0.0324 | 0.1518 | 0.21 |
| nonlinear | 0.0473 | 0.0965 | 0.49 |
| interaction | 0.0359 | 0.1417 | 0.25 |

proxymix has lower PEHE than grf on every process (2x-5x better). On smooth,
low-dimensional structure the parametric mixture is more statistically
efficient than a forest. **S1 met and exceeded** (target: within 10% of grf,
not worse on any).

## Scoring throughput (100,000 units, single core)

| method | time | units/s |
|---|---|---|
| proxymix `proxy_cate` | 0.216 s | 462,963 |
| grf `predict` | 1.279 s | 78,186 |

Speed-up 5.9x. **S3 met** (kill criterion was < 3x; perf target was
>= 1e5 units/s -- 4.6e5 achieved). Below the aspirational 10x headline; the
closed-form O(K) scoring is realised by the vectorised serving path.

## Hillstrom (MineThatData) -- Womens E-Mail vs No E-Mail, outcome = visit

n_train = 29,885, n_test = 12,808, 9 covariates. Qini coefficient (higher is
better):

| method | Qini | AUUC |
|---|---|---|
| **proxymix** | **53.18** | **212.42** |
| T-learner | 50.47 | 209.71 |
| S-learner | 48.57 | 207.81 |
| X-learner | 46.36 | 205.60 |
| grf | 39.97 | 199.21 |

proxymix has the **best** Qini and AUUC. ATE: proxymix 0.0438, grf 0.0427,
DoubleML(IRM) 0.0444 -- proxymix's average effect matches the debiased
DoubleML estimate to within 1.4%. **S2 met and exceeded on Hillstrom.**

## Criteo Uplift v2.1 -- outcome = visit (an honest loss)

Subsampled 250,000 of 25,309,483 rows (logged, not silent); n_train = 150,000,
n_test = 100,000, 12 covariates. Qini coefficient:

| method | Qini | AUUC |
|---|---|---|
| X-learner | 249.19 | 657.64 |
| grf | 243.98 | 652.44 |
| T-learner | 238.56 | 647.02 |
| proxymix | 228.38 | 636.83 |
| S-learner | 216.80 | 625.25 |

On Criteo proxymix is **4th of 5** -- 8.4% below the best (X-learner), which
misses the within-5%-of-best target. Its ATE (0.0038) underestimates the
DoubleML(IRM) ATE (0.0085). On a 12-dimensional, very-low-event-rate
(~4.7% visit) target with weak, diffuse heterogeneity, the Gaussian-mixture
reading is competitive but not best -- tree ensembles capture the sparse signal
better. **We report this loss; we do not claim universal dominance** (spec
section 7).

## Summary

proxymix's decision layer is strongest where the mixture structure matches the
data: smooth, lower-dimensional problems (synthetic, Hillstrom), where it
matches or beats causal forests and meta-learners while scoring ~6x faster and
supplying a closed-form identification audit. On high-dimensional, sparse-signal
problems (Criteo) it is competitive but not best. The edge is interpretable,
auditable, fast closed-form uplift -- not universal accuracy dominance.

## SA gate (spec section 16)

| Kill criterion | Verdict |
|---|---|
| (a) PEHE materially worse than grf on a majority of DGPs | NOT triggered -- proxymix beats grf on all 4 synthetic DGPs |
| (b) closed-form speed advantage < 3x | NOT triggered -- 5.9x |
| (c) identification cannot be made honest | NOT triggered -- refusals + confounding gap + reported loss |
| (d) `R CMD check` cannot reach 0/0/0 | NOT triggered -- 0/0/2 environmental only |

Verdict: **proceed.** The Criteo result is a bounded, reported loss, not a kill
condition. The two operators stay Tier 1.
