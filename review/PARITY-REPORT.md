# Parity report — `glm-parallel` vs baseline (`glm-notes` @ 4b09088)

> Post-Phase-3a note: the likelihood accumulation now lives in the C++ kernel
> `pairLogLik()` (commit 0035346). Its output order matches the old R
> arithmetic, and the golden-parity suite (`test-golden-parity.R`, tolerance
> 1e-10) re-certified sections B's checks against the same baseline fixtures —
> the numbers below refer to the Phase-1 implementation; the golden suite is
> the living guarantee., and vs upstream `master` (@ fe9551c)

Verification of the performance-optimised implementation against the
pre-optimisation baseline, and against the repository's upstream `master`
branch. Harness: `review/verify-parity.R` (modes `old`, `new`, `compare`; an
optional label argument lets any build stand in for either side, e.g.
`Rscript review/verify-parity.R old <outdir> master`). The implementations run
in separate R processes (`R_LIBS` selects the build, all named `FlexRL`), a
third process compares the dumps without loading FlexRL. Shared state/problem
code: `tests/testthat/helper-golden.R` + `review/parity-problems.R`. Raw dumps
and logs: `/mnt/halveGig/kaya/Broncode/FlexRL-bench/parity/` (outside the
repo).

* old = baseline build of branch `glm-notes` @ 4b09088 (`FlexRL-bench/rlib-old`)
* new = current build of branch `glm-parallel` (`FlexRL-bench/rlib-new`)
* master = upstream build of `origin/master` @ fe9551c
  (`FlexRL-bench/rlib-master`, worktree `/mnt/halveGig/kaya/Broncode/FlexRL-master`)
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

---

# vs upstream master (fe9551c)

## Topology — what master-vs-parallel actually compares

`origin/master` @ fe9551c is `77daff8` + 7 commits that touch only
`vignettes/`, `.gitignore`, `.Rbuildignore`, pkgdown files and DESCRIPTION
(the DESCRIPTION diff is packaging metadata only: `VignetteBuilder: knitr,
quarto` + Suggests reorder; R/, src/, NAMESPACE, tests are untouched). So
**master's code is identical to `77daff8`**, i.e. `glm-notes` @ 4b09088 minus
the local review commit. That commit added two behavioural fixes (plus
defensive guards):

* **#A — eta-indexing bias in the linked-pair H-sampler** (`simulateH`,
  src/FlexRLfunctions.cpp): master expands `eta_choice =
  eta_k[choice_set(_,1) - 1]` — eta indexed by the **B** grid value — while
  the model (paper Eq. 4, and `simulateD()`'s own likelihood, which uses
  `eta[h^A]`) indexes eta by **A**. The fix indexes by A.
* **#I-a — interned-constant corruption in `sampleD`** (src/): master writes
  through the `LogicalVector& sumRowD/sumColD` arguments; when a caller passes
  a length-1 logical (R interns it as the global TRUE/FALSE constant), the
  write corrupts R's global constant session-wide. The fix clones the
  arguments. Trigger in master's R code: `sampleD(..., sumRowD > 0, sumColD >
  0)` with `nA == 1` or `nB == 1`.
* **0-column guard in `pasteIntoPattern`** (src/): with **no stable PIVs** the
  pattern matrix has 0 columns and master reads `A(i,0)` out of bounds; the
  fix returns all-equal patterns (blocking degenerates to the full Cartesian
  product). Master's R-side `simulate_data()` also has an `x`-vs-`k`
  column-index bug in the structured-PIV instability loop (stale loop variable
  from earlier loops) — latent whenever the only structured PIV is the last
  PIV, as in every problem used here.

Consequently master-vs-parallel = **review bug fixes + performance changes
combined**. The clean attribution control is master-vs-`glm-notes` (isolates
the fixes); the notes-vs-parallel comparison (sections A–D above) already
established that the optimisation alone is seed-noise-only.

Master was built in its own worktree/library (`FlexRL-bench/rlib-master`),
confirmed to load with the pre-`unif` kernel signatures (`"unif" %in%
names(formals(sampleNL))` = FALSE, `sampleD` has the `loglik` formal, no
`buildCandidatePairs`). The parity dump ran as `old` mode labelled `master`
(one process, sections A–D; no crash, no hang, no erratic behaviour — see
#I-a note below). Comparisons: `compare` mode with labels `master`/`new` and
`master`/`old`; logs `parity-compare-master-vs-{new,old}.log` in the parity
directory.

## A. Harness sanity — exact

All 8 dumped `simulate_data()`/`prepare_data()` objects `identical()` between
master and new (and between master and old — same variant). The latent
`simulate_data` x-vs-k bug above does not manifest: for a single structured
PIV that is also the last PIV, the stale `x` equals `k`, so the executed
statements are identical.

## B. Deterministic kernel parity — exact except the zero-stable-PIV degenerate state

Five of the six fixed states are **bitwise identical master-vs-new and
master-vs-notes** (max |dLLL| = max |dLLA| = max |dLLB| = 0, log_lik spread 0,
candidate-pair SETs equal; master-vs-notes also shares the row order, since
both use the string pipeline).

`tiny_links` **FAILS against both** — this is a genuine master defect, not
noise:

* The tiny state has **no stable PIVs** (its single PIV is structured), so
  blocking runs `pasteIntoPattern` on a 0-column matrix. Master reads out of
  bounds and returns garbage patterns: its candidate set is the 2 diagonal
  pairs `{(1,1), (2,2)}` instead of the correct full Cartesian product of
  2 × 3 = 6 pairs (what glm-notes' guard and the optimised
  `buildCandidatePairs()` both produce).
* The state's links are (1,1) and (2,3). Pair (2,3) is **not in master's
  corrupted candidate set**, so master's LLL matrix has no entry there and the
  state's log-likelihood loses exactly that contribution: master −8.233064 vs
  glm-notes/parallel −13.935662; difference 5.702598 = the missing
  LLL[2,3]. (LLL/LLA/LLB values at the pairs master does contain are bitwise
  identical to the other builds'.)

## C. Tiny-problem stationary distribution — master far off, dominated by the candidate-set defect

`tiny_gibbs_frequencies(n_sweeps = 30000, thin = 3)`, seed 2024, 10000 thinned
draws per version:

| quantity | value |
|---|---|
| TV(master, exact posterior) | **0.6371** (max cell 0.2677) |
| TV(glm-notes, exact posterior) | 0.0318 (previous run) |
| TV(new, exact posterior) | 0.0294 |
| TV(master, glm-notes) | 0.6372 |
| TV(master, new) | 0.6388 |

Mechanism: the tiny problem has no stable PIVs, so **every** Gibbs sweep in
master recomputes its candidate set through the unguarded `pasteIntoPattern`.
Master's chain visited only **4 distinct configurations out of 13** (the
diagonal-only ones reachable from the corrupted {(1,1),(2,2)} candidate set),
while glm-notes visits all 13 — 9 configurations are unreachable, which alone
accounts for a TV of this magnitude. This number is therefore **not** a clean
measurement of the #A eta bias (which is also active here); it documents the
candidate-set defect.

## #A in isolation — kernel-level measurement (`review/verify-eta-bias.R`)

To separate the eta bias from the candidate-set defect, `sampleL` (the
linked-pair H-sampler kernel) was called directly — no candidate pipeline
involved — for the two canonical tiny-problem cases (S = 0.6, eta = c(0.7,
0.3), phi as `tiny_pars()`), 200 000 draws per arm, against the exact
conditional P(a, b | G, Δ=1) ∝ eta[a] · [S·1{a=b} + (1−S)/(n−1)·1{a≠b}] ·
P(G^A|a) · P(G^B|b):

| arm | TV(empirical, oracle), case GA=0/GB=1 | case GA=1/GB=0 |
|---|---|---|
| baseline kernel, `eta[h^B]` (master's shipped `simulateH`) | **0.1773** | **0.1792** |
| baseline kernel, `eta[h^A]` (corrected; glm-notes' shipped `simulateH`) | 0.0012 | 0.0019 |
| optimised kernel (factorised, eta[h^A] internal) | 0.0026 | 0.0027 |

* Master's shipped indexing deviates from the exact conditional by TV ≈ 0.18
  (e.g. P(a=2,b=1) = 0.371 vs oracle 0.193 — mass tilted toward the
  eta-heavy B value). The one-line indexing change removes it.
* The corrected arm produces **identical numbers under the master and
  glm-notes builds** (same kernel, same seed) — confirming the only
  behavioural difference between the two builds' H-sampler is the indexing
  line.
* The optimised kernel matches the oracle at MC-noise level (its RNG stream
  differs by design).

Why the bias is invisible on most problems: `survivalpSameH` is 1 for stable
and flexible PIVs (only structured PIVs get a survival model), and with S = 1
only the diagonal (h, h) grid cells carry weight — where `eta[h^B] ==
eta[h^A]`. The bias is a **no-op unless the problem has a structured PIV**
with S < 1. That is why the master-vs-glm-notes seed clouds below are bitwise
identical on the three real-data problems.

## D. StEM seed-cloud parity — master vs new, with the master-vs-glm-notes control

12 seeds (11 … 1111) per problem, same settings as §D above.

| problem | max&#124;z&#124; new vs master cloud | max&#124;z&#124; master vs new cloud | max&#124;z&#124; notes vs master cloud | max&#124;z&#124; master vs notes cloud | clouds master vs notes bitwise identical? |
|---|---|---|---|---|---|
| std | 2.75 | 3.09 | 2.99 | 3.52 | no (structured V4 → #A active) |
| shiwreg | 3.52 | 3.54 | 3.14 | 3.14 | **yes** |
| nltcs | 3.46 | 3.04 | 2.74 | 2.74 | **yes** |
| shiwfull | 2.89 | 3.19 | 2.62 | 2.62 | **yes** |

Delta (x > 0.5) Jaccard, mean (sd, min):

| problem | within-master (reference) | master-vs-new | master-vs-notes |
|---|---|---|---|
| std | 0.623 (0.050, 0.526) | 0.597 (0.058, 0.450) | 0.596 (0.058, 0.447) |
| shiwreg | 0.690 (0.119, 0.429) | 0.669 (0.117, 0.350) | 0.715 (0.142, 0.429) |
| nltcs | 0.549 (0.032, 0.457) | 0.541 (0.034, 0.452) | 0.586 (0.129, 0.457) |
| shiwfull | 0.496 (0.021, 0.454) | 0.498 (0.021, 0.453) | 0.538 (0.141, 0.454) |

(For the three bitwise-identical problems the master-vs-notes "cross" numbers
are self-comparison artifacts — identical clouds; the cross Jaccard is the
within-cloud distribution plus the 12 same-seed diagonal pairs at Jaccard 1,
e.g. shiwreg (66·0.690 + 12·1)/144 = 0.716 ✓.)

The full 12-seed × 4-problem master cloud completed without reduction (no
crash, no hang). Number of links (Delta x > 0.5), master range vs new range:
std 64–78 vs 67–83; shiwreg 24–27 vs 24–29; nltcs 274–303 vs 285–305;
shiwfull 1159–1220 vs 1133–1202 — overlapping ranges throughout.

Attribution:

* **Real-data problems (no structured PIV): master ≡ glm-notes bitwise.** The
  master-vs-new differences there are *purely* the optimisation and are
  seed-noise-level, matching the notes-vs-new reference (the numbers are
  literally the same comparison, since the master and notes clouds coincide).
* **std (structured V4): #A is active, but its effect on the final StEM
  parameters at these settings stays within the seed-noise band**
  (max |z| ≤ 3.52 < 4 in both master-vs-new and master-vs-notes; cross
  Jaccards statistically indistinguishable from the within-cloud reference).
  The kernel-level bias is real and large (TV ≈ 0.18 per linked-pair draw,
  above), but it averages out over 10 StEM iterations × 50 posterior samples
  at this problem size.
* **tiny (structured V1, S = 0.6, and additionally zero stable PIVs): both
  defects active; the candidate-set defect dominates** (§C).

## Benchmark

`R_LIBS=rlib-master Rscript review/bench-stem.R` (scale 1: nA=10000,
nB=15000, n_links=7500, 22 sweeps), main repo dir:

| build | StEM end-to-end | per sweep | profiled run top self-time |
|---|---|---|---|
| master @ fe9551c | **241.24 s** | 10.965 s | `make.unique` 63.8% (`simulateD` 94.4% total) |
| glm-notes @ 4b09088 (baseline) | 238.7 s | 10.85 s | `make.unique` 63.9% (`simulateD` 94.3% total) |
| glm-parallel (current) | 8.29–8.39 s | 0.38 s | `simulateD` 72% |

Master ≈ glm-notes as expected (the review commit was correctness-only; ~1%
is run noise). The `make.unique` string-pipeline bottleneck that Phase 1.4
removed is intact in master.

## Verdict (vs master)

* Deterministic kernels: **bitwise exact everywhere except the
  zero-stable-PIV degenerate state**, where master's unguarded
  `pasteIntoPattern` corrupts the candidate set (fixed by the review commit's
  0-column guard). All observed likelihood values agree bitwise on every pair
  master does compute.
* Stochastic differences between master and glm-parallel are **fully
  attributable to the review-phase bug fixes**: on problems without a
  structured PIV master is bitwise identical to glm-notes, so the entire
  master-vs-parallel difference there is the (seed-noise-level) optimisation
  effect; on the structured-PIV problems the additional #A effect is
  kernel-proven (TV ≈ 0.18 on the linked-pair conditional) but stays below
  seed noise at the parameter level for std, and is swamped by the
  candidate-set defect on tiny. glm-notes/glm-parallel remain **closer to the
  exact posterior** everywhere both could be measured (TV 0.032/0.029 vs
  0.637 on tiny).
* #I-a did **not** trigger in any of these runs: it needs a length-1
  `sumRowD`/`sumColD` (i.e. a single-record file), and every problem here has
  ≥ 2 records per file (`sumRowD > 0` is then length ≥ 2, not an interned
  constant). No master process crashed, hung, or behaved erratically; all
  runs completed in single passes. The bug remains real (demonstrated by the
  review-phase session-corruption repro) but is outside the envelope of these
  workloads.
* The goldens must **not** be regenerated from master: master carries the #A
  bias and the unguarded `pasteIntoPattern`; the fixtures correctly encode the
  fixed baseline.

