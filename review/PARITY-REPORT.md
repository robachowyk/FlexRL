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

`tiny_gibbs_frequencies(n_sweeps = 30000, thin = 3)`, seed 2024 (as
`test-exact-posterior.R`), 10000 thinned draws per version:

| quantity | value | test tolerance |
|---|---|---|
| TV(old, exact posterior) | **0.0318** (max cell 0.0126) | TV < 0.05, cell < 0.025 |
| TV(new, exact posterior) | **0.0294** (max cell 0.0115) | TV < 0.05, cell < 0.025 |
| TV(old, new) | **0.0186** | same order as MC error (< 0.05) |

The old-vs-new TV distance is of the same order as each version's distance to
the exact posterior — i.e. the two samplers differ only by Monte-Carlo noise.

## D. StEM seed-cloud parity — pass

12 seeds (11, 101, 202, 303, 404, 505, 606, 707, 808, 909, 1010, 1111) per
problem, per version. Settings: synthetic std problem uses the test suite's
`std_stem_args()` (10/5 StEM, 8/4 Gibbs, 50 post-samples); the real-data
problems use scaled-down vignette-style settings (6/3 StEM, 8/4 Gibbs,
20 post-samples). Real-data problems (vignettes' own data and configs):

* `shiwreg`: `ex-reg-SHIW-16.csv` / `ex-reg-SHIW-20.csv`, regional-SHIW config
  (ANASCI, SESSO, STUDIO, STACIV all stable) — 42 x 49 records;
* `nltcs`: `ex-reg-NLTCS-82.csv` / `ex-reg-NLTCS-94.csv`, block-vignette config
  (SEX, DOB_MONTH, DOB_YEAR stable; STATE flexible) — 491 x 901 records;
* `shiwfull`: `SHIW2016.csv` / `SHIW2020.csv`, deduplicated by ID as the
  vignette does, then subsampled to the first 2500 rows of each file
  (deterministic; full files would make the baseline seed-cloud take hours),
  old block-vignette full-SHIW config (ANASCI, SESSO, STUDIO, STACIV, IREG,
  all stable) — 2484 x 2495 records.

No vignette real-data config contains a structured PIV (they use stable +
flexible only), so the survival/alpha machinery is exercised on real data via
the flexible-PIV blocking and on synthetic data via the structured V4/tiny
problems — mirroring exactly what the vignettes themselves cover.

Per-parameter agreement: max |z| of each new seed against the old cloud's
per-component mean/sd, and vice versa (scale: the larger of the two clouds'
sds; 1e-7 absolute floor). A value ≲ 4 means "looks like just another seed".

| problem | max&#124;z&#124; new vs old cloud | max&#124;z&#124; old vs new cloud | Delta Jaccard old-vs-old (reference) | Delta Jaccard old-vs-new |
|---|---|---|---|---|
| std | 2.80 | 3.59 | mean 0.577, sd 0.056, min 0.446 | mean 0.588, sd 0.061, min 0.426 |
| shiwreg | 3.52 | 3.54 | mean 0.690, sd 0.119, min 0.429 | mean 0.669, sd 0.117, min 0.350 |
| nltcs | 3.46 | 3.04 | mean 0.549, sd 0.032, min 0.457 | mean 0.541, sd 0.034, min 0.452 |
| shiwfull | 2.89 | 3.19 | mean 0.496, sd 0.021, min 0.454 | mean 0.498, sd 0.021, min 0.453 |

* The Jaccard overlap of the x > 0.5 link sets is statistically
  indistinguishable between old-vs-old (different seeds, the reference
  distribution for seed noise) and old-vs-new (the swap looking like just
  another seed).
* Number of links (Delta x > 0.5) per problem, old range vs new range:
  std 61–78 vs 67–83; shiwreg 24–27 vs 24–29; nltcs 274–303 vs 285–305;
  shiwfull 1159–1220 vs 1133–1202 — overlapping ranges throughout.
* Runtime of one scaled StEM run (this problem sizes): old 0.3–19 s,
  new 0.06–0.3 s per seed.

## Verdict

No statistical discrepancy found. Deterministic kernels are bit-identical;
stochastic outputs are indistinguishable at Monte-Carlo level on synthetic and
real (vignette) data. The permanent regression suite built from these results
lives in `tests/testthat/test-golden-parity.R` with fixtures under
`tests/testthat/golden/` (generated by `review/generate-goldens.R` from the
baseline build).

## Anomalies / notes

* The z-metric initially blew up on near-degenerate columns (zero-frequency
  eta categories that are exactly 0 in almost every seed: e.g. STATE of the
  NLTCS problem, rare ANASCI birth years). This was a harness artifact, not a
  sampler difference: with sd(reference) = 0 the scale must fall back to the
  other cloud's sd (see `golden_zstats()`). Fixed and documented there.
* `phi` components bound-locked by `fix_mistakes`/`bound_mistakes` are exactly
  constant in BOTH implementations (identical to machine precision) — the
  1e-7 floor handles them.
* The 1e-7 floor for z-scores and the 1e-10 tolerance for deterministic
  comparisons are deliberately far apart: observed deterministic differences
  were exactly 0, so the deterministic tolerance has an order of magnitude of
  headroom for FP-summation-order effects on other platforms.

