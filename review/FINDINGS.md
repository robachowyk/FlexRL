# FlexRL review — FINDINGS

Review of FlexRL 1.0.1 (StEM record linkage + FDP estimation), Oct 2026.
Scope: statistical correctness of the Gibbs/StEM machinery and the C++ kernels,
verified against exact posterior computations and synthetic-data simulations.
All fixes are **local** to this working copy; nothing was pushed upstream.

Evidence artifacts live in `~/tmp/flexrl-check/phase0/` (logs, probe scripts).

---

## 0. Summary

| # | Finding | Severity | Status |
|---|---------|----------|--------|
| A | η indexed by the wrong margin in `sampleH`/`sampleL` | high (statistical) | **FIXED** + tests |
| B | `simulate_data` applied "changes" to the wrong PIV column | medium | **FIXED** + test |
| C | `compute_augmRL_FDP_synth` scaling error | high (statistical) | **FIXED** + verified |
| D/H | `pasteIntoPattern` OOB read on 0-column matrices | crash | **FIXED** + test |
| E | diyar wrapper: entity_id leakage, crash on `change`, weak-links behaviour | medium | documented, tests pin it |
| F | multilink wrapper errors when `duplicates = c(1,0)` | low | documented, test pins it |
| G | config validation accepted character numerics → cryptic crash in StEM | low | **FIXED** + test |
| I-a | `sampleD` wrote through interned R constants (session-wide corruption) | critical | **FIXED** (C++ clone) + test |
| I-b | square `Delta`: `Matrix(0,n,n)` → diagonalMatrix → `summary()` → crash | crash (nA = nB) | **FIXED** + test |
| J | `link_with_fedmatch` returned character row *names* as indices | high (silent wrong results) | **FIXED** + suite |
| K | StEM accepted hand-edited encoded data (NA) → Rcpp INT_MAX crash | medium | **FIXED** (validation) + test |
| L | StEM post-burn-in average includes the burn-in draw (`x[burnin:iter]`) | low (off-by-one) | documented |
| M | Robustness probes: 6 degenerate inputs now fail loudly & clearly | — | battery run |

Independent scientific validation (all passing, in `tests/testthat/`):
- Full Gibbs chain matches the **exact joint posterior** on the tiny problem
  (max cell diff 0.0126, TV 0.0318; `test-exact-posterior.R`).
- Exact conditional check of `sampleL` against the model conditional
  (`test-gibbs-conditionals.R`), plus pure-R replication of the kernels,
  exact pswitch checks, power-iteration checks (`test-gibbs-invariants.R`, 257 expectations).
- Parameter recovery on structured-dynamics data (Sim 1): γ 0.733 vs true 0.75;
  η recovered for V1–V3; **φ biased upward at low typo rates** (see §3.3).

---

## 1. Fixed findings (details)

### #A — η indexed by the wrong margin in the H-sampler (high, statistical)
`sampleH` built the joint (H^A, H^B) choice grid and passed `eta` expanded along
the **B** column, while the model (paper Eq. 4; `simulateD`'s likelihood) defines
η as the marginal distribution of **H^A**. With non-uniform η this reweights the
posterior of every linked pair's true values. Fixed in `src/FlexRLfunctions.cpp`
(`sampleH`, `eta_choice = eta_k[choice_set(_,0) - 1]`).
Verified by `test-exact-posterior.R` (exact joint posterior) and
`test-gibbs-conditionals.R` (conditional-level oracle). Note: the
conditional-level test initially replicated the *pre-fix* calling convention and
kept "passing" only because my early suite harness conflated exit code 0 with
"no failures" (see §5); it now pins the post-fix convention.

### #B — `simulate_data` changed the wrong column (medium)
The loop variable for applying "changes" was stale, so the dynamics of the
*last* configured PIV were applied to the structured PIV. Regression test:
`test-simulate-data.R` ("REGRESSION (FINDINGS #B)").

### #C — `compute_augmRL_FDP_synth` scaling (high, statistical)
The synthetic-data FDP estimator used a mis-scaled ratio; fixed and verified
end-to-end against known-FDP synthetic setups (Sim 2b/2c probes). Post-fix
behaviour on stable-dynamics data is informative and tracks the score-based
estimator; on structured-dynamics data it remains degenerate (see §3.2).

### #D/#H — `pasteIntoPattern` out-of-bounds (crash)
The C++ pattern builder read past the end for 0-column matrices. Guarded;
regression test `test-paste-pattern.R`.

### #I-a — `sampleD` corrupted R's interned logical constants (critical)
`sampleD(..., LogicalVector & sumRowD, LogicalVector & sumColD)` **wrote** through
its reference arguments (`sumRowD(i) = true`). When a caller passed a *scalar*
logical — which happens in `simulateD`'s wrapper whenever `nA == 1` or `nB == 1`,
because `sumRowD > 0` is then a length-1 comparison — R hands over its interned
`TRUE`/`FALSE` constant, and the write flipped **R's global FALSE for the whole
session**: `1 > 2` evaluated to `TRUE` afterwards, and under
`gctorture2(step = 50)` R signalled *"internal FALSE value has been modified"*
in `sampleD`. The cryptic downstream cascade (condition-system errors,
`finally` missing, `loadNamespace`/`dynGet` failures, non-determinism) was all
consequences of this single corruption; cli/OpenMP/Matrix were exonerated.
Fixed by `Rcpp::clone()`-ing both arguments inside `sampleD`
(`src/FlexRLfunctions.cpp`, ~lines 206–290). Regression test:
`test-stem-tiny.R` (direct canary + full StEM runs on (1,5)).

### #I-b — square Delta crashed StEM (`'dimnames' applied to non-array`)
`Matrix::Matrix(0, n, n, sparse = TRUE)` returns a **diagonalMatrix** for square
inputs. `Matrix::summary()` of a diagonalMatrix is a `diagSummary` (a summary of
the *diagonal values*, not the i/j/x triplet list StEM needs), and
`as.data.frame()` on it fails. This crashed StEM whenever
`nrow(encodedA) == nrow(encodedB)` **and** Delta stayed empty or diagonal-only —
exactly the natural (2,2) all-linked case; non-square problems never hit it.
Fixed by constructing Delta (and LLL) with
`Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = ...)`
which is always a general `dgCMatrix`. (Note: `doDiag = FALSE` alone is *not*
enough — it yields a symmetric `dsCMatrix`, which would mirror off-diagonal
assignments.) Regression test: `test-stem-tiny.R` (square shapes (2,2), (3,3)).

### #J — `link_with_fedmatch` returned character row *names* as indices (high)
The wrapper keyed records by `rownames(dataA)` (character!) and returned those
strings as `idxA`/`idxB`. Consequences: (i) range checks and downstream numeric
use break (lexicographic comparison: `"9" <= "150"` is FALSE); (ii) whenever the
encoded data has been filtered — which `prepare_data`'s support restriction
routinely does — rownames are the *original* row numbers, so the returned
"indices" point at the **wrong rows silently**. Fixed: key by `seq_len(nrow())`
and coerce with `as.integer()`. Caught by `test-wrappers.R`.

### #K — StEM crashed cryptically on hand-edited encoded data (medium)
The encoded-PVI contract is "integers in 0..n_values, 0 = missing"
(`prepare_data` maps NA → 0). Feeding `NA` directly into an encoded column (a
plausible power-user edit) reached the C++ Gibbs kernels as `NA_INTEGER`;
`G(i) - 1` then overflowed to INT_MAX and Rcpp aborted with
`Index out of bounds: [index=2147483647; extent=4]`. StEM now validates the
encoded columns up front and fails with a clear message (positional `n_values`
lookup, matching the convention used by the rest of StEM). Regression test:
`test-stem-tiny.R`.

### #G — character numerics accepted by config validation (low, fixed)
`is_numeric_or_na()` accepted strings like `"0.02"` (greps `^[0-9.]+$`), which
later crashed StEM with `non-numeric argument to binary operator` inside
`1 - PIVs_fixMistakes[, k]`. Tightened to require actual numerics; the config is
now rejected with the intended message. Test: `test-prepare-data.R`.

---

## 2. Documented, not fixed (author's call)

### #L — post-burn-in averaging off-by-one (low)
All post-burn-in averages use `x[StEM_burnin:StEM_iter]`, which **includes** the
burn-in-th draw instead of starting one later (`(StEM_burnin + 1):StEM_iter`).
Effect is small for long chains; more consequential for the short chains used in
the tiny-problem tests. Documented; not changed because the "correct" convention
(0- or 1-based semantics of `StEM_burnin`) is a definitional choice for the author.

### #E — diyar wrapper (medium)
Three traps, pinned by `test-wrappers.R`: (i) if the input contains `entity_id`,
the wrapper silently uses it, inflating precision; (ii) it errors on
`prepare_data` output containing the `change` column; (iii) with ≤3 attributes
it returns essentially no links.

### #F — multilink wrapper (low)
Errors when `duplicates = c(1,0)` is requested; pinned by test.

### Stale `_DeltaMap` across chains
The C++ Delta lookup map persists between StEM chains; `initDeltaMap()` is
called per chain in the package code, but any *external* driver (e.g. custom
scripts calling the kernels directly) must remember to reset it. A cross-chain
run without the reset silently reuses stale pattern→index mappings.

### `gibbs_burnin = 0` auto-burnin trap
`gibbs_burnin = 0` triggers an automatic burn-in rather than literally zero;
surprising but intentional; pinned by `test-gibbs-invariants.R`.

### nlminb convergence unchecked
The M-step's `stats::nlminb()` return code (convergence flag) is not inspected;
a non-converged M-step passes silently. Low impact in practice (StEM re-estimates
every iteration) but worth a warning.

### Autocorrelated `n_post_sample` draws
The `n_post_sample` Delta draws are a single continuous chain tail, not thinned
independent draws; fine for posterior-mean estimation, but users treating them
as independent samples (e.g. for interval estimates) will understate uncertainty.

### Vignette hygiene (R CMD check NOTE)
`vignettes/OLDFlexRL-*.Rmd` (titles equal to their file names) collide with the
current vignettes' `%\VignetteIndexEntry` values, causing the
"Duplicated vignette titles" NOTE; `.orig` and `.html` artifacts also ship in
`vignettes/`. Suggestion: move the OLD vignettes out of `vignettes/` (or give
them distinct index entries) and drop the `.orig` files.

---

## 3. Simulation findings (behaviour, not bugs)

### 3.1 Parameter recovery (Sim 1, structured-dynamics data)
γ recovered at 0.733 vs true 0.75; η recovered for V1–V3; the structured PIV's
survival parameters are identified. Full log: `phase0/sim1_recovery.R`.

### 3.2 φ collapse at low typo rates (identifiable limitation)
At typo rates ≲1%, φ is biased **upward** (the EM-style imputation absorbs
mistakes into the η/Δ estimates); it recovers once typo rates reach ~10%.
Related: `FDP_score` is optimistic while φ is collapsed, and `FDP_synth`
degenerates on structured-dynamics data (the `date = 0` handicap in the
synthetic-data path silently yields FDP ≈ 0.000 — worth a warning for users).
On stable-only data `FDP_synth` is informative.

### 3.3 Deranged-input battery (Sim 3, 15 probes)
Final state (after fixes #I-a/#I-b/#J/#K/#G): **9 pass, 6 fail loudly with
clear messages** — no crashes, no silent garbage:
all-NA hand-edited column → clear contract error; character-numerics config →
validation message; n_values = 1 (generator) → clear message; structured-only
config with wrong `cond_hazard_params` length → clear message; `PIVs_config` as
vector → clear message; plus passing degenerate cases: nA=1, nB=1, constant PIV
(post-restriction), disjoint supports, duplicate rows, all-linked, zero links,
50% missingness, nB = 20·nA, entity_id collisions.
Script/log: `phase0/sim3_deranged.R`, `phase0/sim3_results.log`.

---

## 4. R CMD check (--as-cran)

`FlexRL_1.0.1.tar.gz`, R 4.6.1: **1 ERROR, 1 WARNING, 4 NOTEs**.
- ERROR/WARNING (PDF manual): **environmental** — no `pdflatex` on this machine.
- NOTE (compilation flags): distro hardening flags — environmental.
- NOTE (non-standard `FlexRL-manual.tex` in check dir): artifact of the failed PDF step.
- NOTE (**duplicated vignette titles**): real, see §2 (vignette hygiene).

---

## 5. Test suite (final state)

14 files in `tests/testthat/`, run per-file with correct failure accounting:
**443 → 446 expectations passing, 0 failures, 0 errors** (~66 s; the
exact-posterior file dominates at ~45 s).

New/changed files this review:
- `test-stem-tiny.R` — #I-a canary, #I-b square shapes, #K validation.
- `test-gibbs-conditionals.R` — η convention updated to post-#A (was silently
  stale; see §1 #A).
- `test-prepare-data.R` — #G flipped from "pins the hole" to "asserts the fix".
- `test-survival-model.R`, `test-simulate-data.R`, `test-synthesise.R`,
  `test-wrappers.R` — test-code fixes (see below).

**Harness lesson (important for future runs):** the earlier "full suite green"
claims were produced by a harness that equated a zero exit status with
"passing". `test_file()`/`test_dir()` do not exit non-zero per failing
expectation in all configurations, so **six files carried failing expectations
for a while, masked as "ok"** (including a non-existent `expect_finite`
expectation and a `cbind(1, 2)` that builds a 1×2 matrix, not two records).
Always count `failed`/`error` from the result object (see the audit snippets in
`phase0/suite_audit*.log`).

Test-code bugs fixed during triage (for the record): `expect_finite` →
`expect_true(is.finite(...))`; survival-model test built `cbind(1, 2)` (one
record, two covariates) instead of `cbind(1, c(1, 2))` (two records) and had
wrong expected values; `simulate_data`'s `proba_same_H` is `NA` for
stable-only configs (documented contract) — test asserted a matrix;
`prepare_data`'s `cov1`/`cov2` refer to the **first/second argument** (swapped
to covA/covB on file swap) — test had them inverted; `same_mistakes = TRUE`
(default) requires identical per-file mistake rates — the swap-reversal test
needs `same_mistakes = FALSE`; `synthesise` may drop real rows via the support
restriction — test now asserts subset semantics; `synthpop::syn` needs ≥2
columns — single-PIV test used two; BRL fixture with ~20 distinct patterns for
120 records is too ambiguous for a 300-iteration chain (BRL links 60/60 with
discriminating values — wrapper indexing verified correct).
