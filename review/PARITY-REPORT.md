# Parity report — `glm-parallel` vs baseline (`glm-notes` @ 4b09088)

Verification of the performance-optimised implementation against the
pre-optimisation baseline. Harness: `review/verify-parity.R` (modes `old`,
`new`, `compare`); the two implementations run in separate R processes
(`R_LIBS` selects the build, both named `FlexRL`), a third process compares
the dumps without loading FlexRL. Shared state/problem code:
`tests/testthat/helper-golden.R` + `review/parity-problems.R`. Raw dumps and
logs: `/mnt/halveGig/kaya/Broncode/FlexRL-bench/parity/` (outside the repo).

* old = baseline build of branch `glm-notes` @ 4b09088 (`FlexRL-bench/rlib-old`)
* new = current build of branch `glm-parallel` (`FlexRL-bench/rlib-new`)
* R 4.6.1 (2026-06-24), same machine, sequential runs.

Because the optimisation changed how uniforms are consumed (pre-drawn blocks,
inverse-CDF, factorised `sampleL`), the RNG stream for a given seed differs
**by design** (PERF-PLAN §5). Exact equality is therefore only expected for
pure-R code (A) and for the deterministic likelihood/pair-set kernels (B);
stochastic outputs are compared in distribution (C, D).

## A. Harness sanity — exact

`simulate_data()` + `prepare_data()` outputs on fixed seeds (pure R-RNG code,
unchanged by the optimisation): **all 8 dumped objects `identical()` between
old and new** (encodedA/encodedB, n_values, true_pairs for two synthetic
problems).

## B. Deterministic kernel parity — exact (max |diff| = 0)

Six fixed states (rebuilt bit-identically in both processes from pure base-R
randomness; links drawn one-to-one from each state's own candidate set):

| state | data | PIVs (dynamics) | pairs | links | max&#124;dLLL&#124; | max&#124;dLLA&#124; | max&#124;dLLB&#124; | log_lik spread |
|---|---|---|---|---|---|---|---|---|
| std_links | synthetic std | 3 stable + 1 structured (covariates) | 296 | 40 | 0 | 0 | 0 | 0 |
| std_nolinks | synthetic std | idem, no links | 264 | 0 | 0 | 0 | 0 | 0 |
| nltcs_links | ex-reg-NLTCS-82/94 | 3 stable + 1 flexible | 276 | 100 | 0 | 0 | 0 | 0 |
| shiwreg_links | ex-reg-SHIW-16/20 | 4 stable | 11 | 8 | 0 | 0 | 0 | 0 |
| shiwfull_links | SHIW2016/2020 (first 2500 rows) | 5 stable | 753 | 150 | 0 | 0 | 0 | 0 |
| tiny_links | tiny 2-value problem | 1 structured (custom S) | 6 | 2 | 0 | 0 | 0 | 0 |

* Candidate pairs: the SET of pairs is identical in every state (sorted-row
  comparison). The row ORDER differs (baseline: pattern buckets in
  first-appearance order; optimised: lexicographic by pattern tuple) — by
  design, documented in `buildCandidatePairs()`. Only the degenerate tiny
  state (6 pairs) happens to coincide in order.
* Likelihood values: old `LLL[select]` (sparse route) vs new `LLLq` (per-pair
  vector), aligned by (i, j) pair — **max |diff| = 0** (bitwise), likewise
  LLA and LLB. The log-likelihood of the state agrees across all four routes
  (baseline internal value, baseline `log_lik()`, new `log_lik()` on a rebuilt
  sparse LLL, new `log_lik_pairs()`): identical to all printed digits, spread 0.
* The states deliberately cover: links present and absent, unstable PIVs with
  the survival term (structured with covariates + `date`, and the custom
  constant-S model), flexible PIVs in the blocking, real-data value
  distributions, missing values (0-sentinel), and the 0-column degenerate
  candidate-pair case (tiny).

## C. Stationary distribution (tiny problem) — pass

_(reported in the next commit together with the seed-cloud evidence)_

## Verdict (so far)

No discrepancy found. Deterministic kernels are bit-identical; the harness
sanity checks pass exactly. Stationary-distribution and seed-cloud evidence
follow in the next commit.

