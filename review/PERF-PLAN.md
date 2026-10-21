# FlexRL — performance & parallelization plan

Branch: `glm-parallel`. Living document: update the worklog (§8) as work lands.

## 0. Goal & constraints

Goal: cut the wall-clock of `StEM()`'s Gibbs machinery and open the door to using
all cores (dev machine: 32). Constraint hierarchy, in order:

1. **Statistical exactness first.** The Gibbs stationary distribution must not
   change (beyond MC error). `tests/testthat/test-exact-posterior.R` is the
   arbiter — it pins the sampler's frequencies against the exact posterior of a
   fully enumerable tiny problem.
2. **Reproducibility.** `set.seed()` must keep controlling everything, and two
   runs with the same seed must be bit-identical (`test-gibbs-invariants.R`).
3. **CRAN-friendliness.** No heavy new dependencies; OpenMP (base R tooling) is
   the only anticipated addition. Examples and tests default to 1 thread
   (CRAN check allows ≤ 2).

## 1. Inventory: code that iterates pairs / performs updates

| Where | Lines | Iterates over | Mutates | Order-dependent? |
|---|---|---|---|---|
| `src/FlexRLfunctions.cpp` `sampleD()` | 206–289 | candidate pairs `S`, row by row | `_DeltaMap` (global static map), `rowD`, `colD`, `loglik`, `nlinkrec` | **yes** — via one-to-one feasibility + link count only (§2) |
| `src` `sampleNL()` | 302–336 | non-linked records of one PIV | latent true values `H(i)` | **no** — conditionally independent |
| `src` `sampleL()` | 355–419 | linked pairs of one PIV | joint latent `(H^A, H^B)(i)` | **no** — conditionally independent |
| `src` `sampleH()` | 479–536 | PIVs `k` | — | **no** — PIVs independent given links/params |
| `src` `indexPatterns()` / `pairPatterns()` | 42–100 | pattern buckets | builds `select` | no |
| `R/FlexRecordLinkage.r` `simulateD()` | 523–629 | candidate pairs (vectorized R) | `LLA`, `LLB`, `LLL` (nA×nB sparse) | no — pure function of the state |
| `R` `StEM()` loops | 1086–1280 | StEM iters / Gibbs sweeps / post-samples | the chain itself | **yes — inherently sequential**; parallelism can only live *inside* one sweep |
| `R` `naive_linkage()` | 1340+ | pattern buckets | — | out of scope (deterministic, not the sampler) |

Known per-kernel costs — **measured** on the bench problem (`review/bench-stem.R`,
scale 1: nA=10000, nB=15000, n_links=7500, 22 sweeps):

- **Baseline: 238.7 s end-to-end (10.85 s/sweep).** Profiled run: `simulateD`
  accounts for **94.3%** of total time. Inside it:
  - `[.data.frame` 68.1% total, of which **`make.unique` 63.9% self** —
    data-frame row subsetting with ~million-row index vectors
    (`data$encodedA[select[, 1], k, drop=FALSE]` etc., done ~2× per PIV per
    sweep); R's `[.data.frame` builds and de-duplicates the row names of the
    extracted result (micro-repro confirmed: 94% of `[.data.frame` time in
    `make.unique`). Fix: grab plain column vectors once per sweep and index
    them directly — no data-frame machinery in the hot path.
  - `[<-` 11.0% (`replCmat4`/`.TM.repl.i.mat`) — the `LLL[select] <- …`
    sparse sub-assignment.
  - `cbind` 10.5% — hazard-covariate assembly per PIV.
- **"Drawing Delta" ≈ 11 s per post-sample** (the per-link
  `Delta[i, j] <- Delta[i, j] + 1` single-cell dgCMatrix sub-assign, ~1.4 ms
  each) ≈ **46% of end-to-end** at n_post_sample=10. Fix: accumulate triplets
  and build/add one `sparseMatrix` per post-sample (Phase 1.5).
- The earlier hypothesis that `Rcpp::sample()` per record/pair dominates is
  **wrong at this scale**: `sampleH` total = 0.9% (the calls live inside its
  `.Call`). Still worth fixing (it is the prerequisite for threading and grows
  with file size), but it is not the bottleneck.
- `sampleD`'s scan does not even reach the top-15 self-time list at this scale.

Note: `Rprof.out` at the repo root profiles `link_with_BRL` (a competitor-package
comparison run), **not** FlexRL's own sampler — it is not evidence about our hot
path.

## 2. Does the update order matter? (hypothesis, checked)

**H-sampling: order is irrelevant — exactly.** In `sampleNL` each record's draw
conditions only on `(G(i), eta, phi)`; in `sampleL` each linked pair's draw
conditions only on `(GA(i), GB(i), survivalpSameH(i), eta, phi)`. No draw reads
another's output → embarrassingly parallel, and parallel execution is
statistically *identical* (given per-draw independent uniforms).

**D-sampling: order matters, but only through two narrow channels.** The
acceptance odds in `sampleD` simplify. With
`sumlogdensity = log(exp(loglik) + exp(loglikNew))`:

    pswitch = exp(loglikNew - sumlogdensity)
            = exp(Δ)/(1 + exp(Δ)),   Δ = loglikNew − loglik

so **`loglik` cancels entirely** in the acceptance probability. The per-pair
delta is

    Δ_add(q) = −LLB(j) − LLA(i) + LLL(q) + log γ − log(1−γ) − log(nB − nlinkrec)
    Δ_rem(q) = +LLB(j) + LLA(i) − LLL(q) − log γ + log(1−γ) + log(nB − nlinkrec + 1)

which contains only quantities **fixed for the whole sweep** (`LLL`, `LLA`,
`LLB`, `γ` — `gamma` arrives as `pLink = rep(gamma, nA)`) plus the running link
count. Therefore pair *q*'s decision depends on earlier pairs only via:

1. **feasibility** — `rowD(i)` / `colD(j)`, i.e. the one-to-one constraint, and
2. **`nlinkrec`** — the bipartite-matching normalizer term.

Everything else (the expensive-looking log-likelihood bookkeeping) is
pair-local and can be hoisted/precomputed. This is what makes both the cheap
optimizations (§3 option C) and the parallel schemes (§3 options A/B) possible.

Corollary used later: `sampleD` does not actually need `loglik` as an *input* —
it only accumulates it for output. The final value is reproducible from scratch
as (matching + unlinked-record terms) + Σ LLL over finally-linked pairs + Σ
deltas of accepted flips.

## 3. Options for `sampleD`

| Option | Per-pair scan cost | Sequential span | Exact? | Effort |
|---|---|---|---|---|
| today | ~6 transcendentals + 1 `runif` | O(n_pairs) | yes | — |
| **C: precompute deltas, keep the serial scan** | 1 exp + 1 `runif` + flops (count term maintained incrementally) | O(n_pairs) — unchanged | yes — same kernel | hours |
| **A: exact conflict-free batching** | proposal uniform + parallel reduction per batch; 1 acceptance uniform / batch | O(#batches) + O(n_pairs / cores) | yes — valid π-invariant block Metropolis | days |
| B: stale-count parallel draws | as A | fully parallel | **no** — bias through the matching-count term | days |

Notes:

- **Option A details.** Color the pair-conflict graph (two candidate pairs
  conflict when they share an A- or a B-record). Per pattern bucket the
  candidate pairs form the edge set of K_{p,q}; edge-coloring gives batches by
  `(rank_A + rank_B) mod max(|A|,|B|)`. Within a batch pairs are disjoint, so
  flips cannot collide through feasibility; the only remaining coupling is
  `nlinkrec`, whose joint contribution to a batch flip depends only on
  `(k_add, k_rem, nlinkrec)`:
      ΔLL(batch) = Σ_{q added} Δ_add(q) + Σ_{q removed} Δ_rem(q)
                   − [logCfg(nlinkrec + k_add − k_rem) − logCfg(nlinkrec)]
  with `logCfg(n) = Σ_{r<n} log(nB − r)`. Proposal: flip each actionable pair
  w.p. ½ (symmetric — the actionable set is invariant within a conflict-free
  batch); accept the whole batch with `min(1, exp(ΔLL))`. Per-pair deltas are
  precomputed once per sweep in parallel.
- **Option A risk:** #batches ≈ largest pattern-bucket size; heavily skewed
  buckets (many records sharing one true stable pattern) erode the gain.
- **Option B is rejected** for now: freezing `nlinkrec` per batch biases the
  kernel, visibly on small `nB` (the tiny validation problem has nB = 3, where
  the count term is huge). The package's selling point is exactness and
  `test-exact-posterior.R` would flag it.
- **Decision rule** (applied after Phases 1–2): implement A only if the
  `sampleD` scan is still > 15–20% of total sweep time once the quick wins have
  landed; otherwise the complexity is not worth single-digit-percent gains.

## 4. Phased roadmap

### Phase 1 — quick wins (single-threaded; distributions identical; RNG streams change)

1. **Inverse-CDF categorical sampler fed by pre-drawn uniforms**, replacing
   `Rcpp::sample()` per record/pair in `sampleNL` / `sampleL`. Uniforms are
   drawn per sweep from R's RNG (sequential, `set.seed`-controlled) and passed
   down into the kernels — also the prerequisite for threading (R's RNG is not
   thread-safe).
2. **Factorize `sampleL`.** Sample `H^A` from its O(nval) marginal
      P(a) ∝ eta[a] · pGA(a|GA) · [S·pGB(a|GB) + (1−S)/(nval−1)·(1−pGB(a|GB))]
   then `H^B | H^A` ∝ pjoint(a,b) · pGB(b|GB) — an exact factorization of the
   same joint conditional (test-gibbs-conditionals' oracle still applies).
   Per pair O(nval²) → O(nval); per-pair `NumericVector` allocations and the
   `ExpandGrid` support grid disappear; `sampleL` returns the joint (H^A, H^B)
   matrix directly.
3. **Kill the nA×nB sparse `LLL`.** `simulateD` accumulates a plain per-pair
   vector `LLLq` (one pass, no `Matrix` sub-assignment, single `anyNA`).
   `sampleD` drops its `loglik` input and reconstructs the value from scratch
   (§2 corollary), returning `linksQ` (q-indices of the final links) so the
   R-side sanity check can recompute the likelihood without touching a sparse
   matrix. The exported `log_lik()` keeps its signature (public API);
   `simulateD` switches to an internal `log_lik_pairs()`.
4. **`buildCandidatePairs()`.** One C++ call replacing the
   `pasteIntoPattern → unique → factor → indexPatterns → pairPatterns` string
   pipeline in `simulateD` (lexicographic tuple sort → deterministic bucket
   order; handles the 0-column degenerate case). The old pipeline stays in
   place for `naive_linkage()`.
5. **Delta accumulation via triplets** *(discovered by the baseline profile —
   ~46% of end-to-end)*: in the `n_post_sample` loop, replace the per-link
   `Delta[i, j] <- Delta[i, j] + 1` single-cell sparse sub-assignment with one
   `Matrix::sparseMatrix(i, j, x = 1, dims = …)` build + add per post-sample.

### Phase 2 — measure

- `review/bench-stem.R`: fixed synthetic problem, end-to-end `StEM()` time plus
  a per-kernel Rprof breakdown. Run before/after each phase; numbers recorded
  in §8.
- Apply the §3 decision rule for Option A.

### Phase 3 — OpenMP (conditional on Phase 2)

- `src/Makevars` with `$(SHLIB_OPENMP_CXXFLAGS)`; `threads` argument to
  `StEM()` (default 1).
- Parallel over records/pairs — exactly order-independent: `sampleNL`,
  `sampleL`, `LLLq` construction, `buildCandidatePairs` bucket fill.
- `sampleD` exact batching (Option A) if the decision rule triggered.
- Uniforms still pre-drawn sequentially per sweep; thread↔slot assignment is
  static → bit-reproducible run-to-run.

## 5. RNG strategy (decided)

- All randomness keeps flowing from R's RNG: `Rcpp::runif()` blocks are drawn
  sequentially per sweep inside the C++ kernels, then consumed deterministically
  (record/pair `q` uses its own slot; factorized `sampleL` consumes 2 slots per
  pair). No thread touches the RNG in Phase 3 either — the block is drawn before
  each parallel section.
- `set.seed()` semantics preserved; two runs at the same seed remain
  bit-identical (static thread↔slot assignment keeps parallel runs
  deterministic).
- **Streams are not comparable to FlexRL ≤ 1.0.1**: inverse-CDF consumes
  uniforms differently than `Rcpp::sample()`, and factorized `sampleL` uses two
  uniforms per pair instead of one. Distributions are identical; tolerance-based
  tests are unaffected; reproducibility tests compare two runs of the *same*
  code. → record in `NEWS.md`.

## 6. Validation strategy

- Existing suite, unchanged tolerances:
  - `test-exact-posterior.R` — stationary distribution vs the exact posterior
    on the tiny problem. **The arbiter for any kernel change.**
  - `test-gibbs-invariants.R` — one-to-one constraint at all times;
    bit-reproducibility of two same-seed runs.
  - `test-gibbs-conditionals.R` — kernel-level oracles for `sampleL` /
    `sampleNL` (updated to the new signatures in Phase 1).
  - `test-stem-tiny.R` — interned-constant regression canary (`1 > 2`).
- **Golden-parity regression suite** (added after the parity study, see
  `review/PARITY-REPORT.md`): `tests/testthat/test-golden-parity.R` validates
  the current implementation against fixtures generated by the
  pre-optimisation baseline (`review/generate-goldens.R`, baseline
  `glm-notes` @ 4b09088, stored under `tests/testthat/golden/`, ~264 KB):
  - fixed-state likelihood artifacts and candidate-pair SETs: exact
    (tolerance 1e-10; observed diff 0) — LLL aligned by (i, j) pair because
    the two implementations order candidate pairs differently by design;
  - 12-seed StEM clouds on the synthetic std problem and three real-data
    problems from `vignettes/` (ex-reg-SHIW, ex-reg-NLTCS, SHIW2016/2020
    first 2500 rows): per-parameter |z| <= 5 against the baseline cloud
    (observed 2.8–3.6) and Delta (x > 0.5) Jaccard mean within 0.05 of the
    old-vs-old seed-noise reference (observed |d| <= 0.021);
  - tiny-problem frequencies (30k sweeps, thin 3, seed 2024): TV to the
    baseline frequencies <= 0.05 (observed 0.0186) plus the
    exact-posterior checks.
  The goldens make the RNG-stream change (§5) permanently measurable: any
  future kernel change that shifts the stationary distribution fails the
  suite even though streams keep differing.
- New checks:
  - `linksQ` ↔ `links` consistency: every returned link's q must map back to
    the same (i, j) pair; `length(linksQ) == nrow(links)`.
  - Benchmark numbers recorded in §8 (no hard thresholds; watch for regressions
    between phases).
- Manual: `Rprof()` before/after each phase on the bench problem.

## 7. Risks / notes

- `sampleL` factorization changes RNG consumption (2 per pair instead of 1) —
  distribution-preserving, but the conditional-level test must be re-checked
  against the same oracle.
- `sampleD`'s loglik reconstruction changes the FP summation order (~1e-12
  drift); the R-side sanity check rounds to 3 decimals, as before.
- Zero-weight edge cases (e.g. all-zero `eta` after a degenerate M-step):
  inverse-CDF falls back to the last index instead of R's `sample()` error.
  The package already warns about zero-frequency categories; documented, not
  guarded further.
- `_DeltaMap` remains a global static — thread-hostile; restructure in Phase 3
  if Option A triggers (batch flips touch disjoint keys, so a per-batch
  sequential apply or per-row locking both work).
- `pasteIntoPattern` / `indexPatterns` / `pairPatterns` / `ExpandGrid` /
  `generateSequence` stay exported: `naive_linkage()` and the tests use them.
- `Rprof.out` at the repo root is untracked and unrelated — leave alone.

## 8. Worklog

Baseline (branch `glm-parallel`, pre-change, R 4.6.1, bench scale 1):
- StEM end-to-end **238.7 s** for 22 sweeps (**10.85 s/sweep**), nA=10000,
  nB=15000, n_links=7500, 4 PIVs.
- Profile: `simulateD` 94.3% of total; `[.data.frame` 68.1% (`make.unique`
  63.9% self); Matrix `[<-` 11.0%; `cbind` 10.5%; "Drawing Delta" ≈ 11 s per
  post-sample (~46% of end-to-end).
- Consequence: Phase 1 priority order adjusted — the LLLq rewrite (1.3) and
  the Delta-triplet fix (1.5) target ~90% of measured time; `Rcpp::sample`
  replacement (1.1/1.2) stays for exactness-neutral speedup and threading
  readiness.

- [x] Phase 1.1 — inverse-CDF sampler + `sampleNL` on pre-drawn uniforms
- [x] Phase 1.2 — factorized `sampleL` (joint (H^A, H^B) return)
- [x] Phase 1.3 — per-pair `LLLq`, named-column access, `sampleD` loglik
      reconstruction + `linksQ` (kills the `make.unique` bottleneck)
      → **mid-phase bench: 238.7 s → 20.2 s end-to-end (11.8×); profiled run
      132.4 s → 10.2 s. Remaining: Delta drawing `[<-` 37% (Phase 1.5),
      `simulateD` 45% total — of which the string candidate pipeline and
      `cbind` covariate assembly (Phase 1.4), `sampleD` scan 3%.**
- [x] Phase 1.4 — `buildCandidatePairs()` replaces the string pipeline
- [x] Phase 1.5 — Delta accumulation via triplets (kills "Drawing Delta")
- [x] Phase 2 — bench numbers:

| Measurement | Baseline | After Phase 1 |
|---|---|---|
| StEM end-to-end (22 sweeps) | 238.7 s | **8.29 s (28.8×)** |
| per sweep | 10.85 s | 0.377 s |
| profiled run (2 iters, 5 draws) | 132.4 s | 4.8 s |
| dominant self-time | `make.unique` 64% (data-frame row subsetting) | `simulateD` 72% (R-level per-PIV passes) |
| `sampleD` scan share | < 3% | 7.1% |

      **Decision rule outcome (§3): `sampleD` at 7.1% < 15–20% ⇒ Option A
      (exact batching) NOT triggered.** The remaining time is `simulateD`'s
      vectorized per-PIV passes over the candidate pairs (LLLq accumulation,
      survival-model `cbind`/`ph`), i.e. memory-bandwidth-bound R arithmetic —
      the natural Phase 3 lever if wanted is moving that loop to C++ (and then
      OpenMP over pairs), not the batch sampler.
- [x] Parity verification (old vs new), harness `review/verify-parity.R`,
      results in `review/PARITY-REPORT.md`:
      - A: `simulate_data()`/`prepare_data()` outputs `identical()` old vs new.
      - B: 6 fixed states (synthetic std ± links, ex-reg-NLTCS, ex-reg-SHIW,
        SHIW2016/2020 first 2500 rows, tiny): candidate-pair SETs equal, LLL /
        LLA / LLB bitwise identical, `log_lik()` = `log_lik_pairs()` = baseline
        internal value (max |diff| = 0).
      - C: tiny problem, 30k sweeps/thin 3/seed 2024: TV(old, exact) = 0.0318,
        TV(new, exact) = 0.0294, TV(old, new) = 0.0186 (MC noise).
      - D: 12-seed StEM clouds, synthetic std + 3 real-data problems:
        max |z| of each cloud against the other 2.80–3.59; Delta (x > 0.5)
        Jaccard old-vs-old vs old-vs-new indistinguishable (means 0.50–0.69,
        |Δmean| <= 0.021).
      - Harness artifact fixed en route: cloud z-scores must scale
        near-degenerate columns (zero-frequency eta categories) by the larger
        of the two clouds' sds, and per-column standardisation needs `sweep()`
        (plain matrix `/` recycles column-major).
- [x] Golden-parity regression suite: `review/generate-goldens.R` (run against
      the baseline build) -> `tests/testthat/golden/` (264 KB);
      `tests/testthat/test-golden-parity.R` validates the current package
      against the stored baseline outputs (exact kernels, cloud z/Jaccard,
      tiny TV). Suite after adding it: **550 PASS / 0 FAIL** (~31 s).
- [x] Comparison vs upstream `master` (@ fe9551c; harness label mode +
      `review/verify-eta-bias.R`), results in `review/PARITY-REPORT.md`
      "vs upstream master" section. Topology: master = `77daff8` (code) =
      glm-notes minus the review fixes, so master-vs-parallel = fixes + perf.
      - A: exact. B: 5/6 states bitwise identical vs new AND vs notes;
        `tiny_links` fails vs both — master's unguarded `pasteIntoPattern`
        reads out of bounds with zero stable PIVs and returns a 2-pair
        candidate set ({(1,1),(2,2)} instead of the 2x3 product), dropping the
        state's (2,3) link (log-lik −8.233 vs −13.936; Δ = the missing
        LLL[2,3] = 5.7026).
      - C tiny (30k sweeps, seed 2024): TV(master, exact) = 0.6371 vs
        TV(new, exact) = 0.0294 (notes 0.0318) — master visited only 4/13
        configurations (candidate-set defect), not a clean #A measurement.
      - #A isolated (`verify-eta-bias.R`, sampleL vs exact conditional,
        200k draws): master's shipped eta[h^B] indexing TV = 0.177/0.179;
        corrected eta[h^A] arm (identical under master and notes builds) and
        the optimised kernel both at MC noise (TV 0.001–0.003). Bias is a
        no-op for stable/flexible PIVs (S = 1 ⇒ diagonal-only support where
        eta[h^B] == eta[h^A]) — hence master ≡ notes BITWISE on the three
        real-data problems, and on std the effect stays below seed noise
        (max |z| ≤ 3.52 vs new and vs notes).
      - Bench master: **241.24 s** end-to-end (baseline 238.7 s, current
        8.29 s), profile `make.unique` 63.8% — the Phase-1.4 target intact
        upstream.
      - #I-a did not trigger anywhere (needs nA = 1 or nB = 1; all problems
        have ≥ 2 records per file); no master run crashed or hung.
- [x] Phase 2b — big-scale profiling (`bench-stem.R 5` and `10`):

| | scale 5 (50k/75k, 37.5k links, ~14M pairs) | scale 10 (100k/150k, 75k links, ~28.5M pairs) |
|---|---|---|
| profiled 13-sweep run | 67.5 s (5.19 s/sweep) | 258.7 s (19.9 s/sweep) |
| end-to-end (22 sweeps) | ≈ 114 s (derived) | ≈ 438 s (derived) |
| `simulateD` R-level accumulation | **84.6%** (self 72.5 + cbind 5.5 + as.vector 2.8 + S 2.4 + %*% 1.5) | **84.6%** (self 73.1 + cbind 5.1 + as.vector 2.4 + S 2.4 + %*% 1.6) |
| `sampleD` scan | 12.6% | 13.2% |
| everything else | ≤ 2% (buildCandidatePairs 0.9, sampleH 0.4, Delta/M-step ≈ 0) | ≤ 2% (idem) |

      Shares are scale-stable ⇒ the optimisation targets do not change with
      problem size. The R temporaries the C++ move would eliminate are
      ~10 × nPairs × 8 B (~2.3 GB at scale 10) — the move also extends the
      feasible problem size (61 GB machine RAM).
- [ ] Phase 3 — **revised order** (replaces §4 Phase 3; per §3 the batch
      scheme is not the first lever):
      1. [x] **3a — move the `simulateD` per-PIV likelihood accumulation to
         C++** — done in commit 0035346 (`pairLogLik()`): one call per sweep
         computes LLA/LLB/LLLq for all PIVs; unstable PIVs assemble
         covariates/times in C++ and call `model_dynamics$S` once per PIV
         (custom survival models keep working). Operation order matches the
         old R arithmetic left-to-right ⇒ bitwise-identical values
         (golden-parity suite passes). Measured:
         end-to-end 8.29 → 2.60 s (scale 1), ≈114 → 36.7 s (scale 5),
         ≈438 → 136.4 s (scale 10) — 3.1–3.2× at every scale.
         Profile at scale 10 now: `sampleD` **44.6%**, `pairLogLik` 28.2%,
         survival-S R-internal work (`as.vector`+`%*%`) ~22%,
         `buildCandidatePairs` 2.9%, `sampleH` 0.6%.
         Optional 3a+: compute S in C++ for the built-in models (saves the
         ~22% S-internal cost; custom user models keep the R callback).
      2. [x] **3b — OpenMP** — done in commits 2e59a0c + 5d3e924 (2026-10-21):
         `src/Makevars`/`Makevars.win` pick up `$(SHLIB_OPENMP_CXXFLAGS)`
         (empty where unsupported ⇒ still builds single-threaded); pragmas
         guarded `#ifdef _OPENMP`; `threads` argument (default 1) on
         `StEM()`/`simulateD()`/`pairLogLik()`, threaded through `gibbs_step`
         and the post-sample Delta loop. Every pair writes its own slot, no
         reductions ⇒ output bitwise identical for ANY thread count
         (`test-threads.R`: threads=1 vs 2 StEM-identical). The R survival
         callback stays outside the parallel regions. **Honest corrections
         vs the plan text:** `sampleNL`/`sampleL` were NOT threaded (0.6% of
         a sweep — pointless); measured scaling is DRAM-bandwidth-bound, not
         3–6×: micro-bench on a fixed 97M-pair state gives stable-PIV pass
         0.181 s (1 thread) → 0.096 s (32 threads) ≈ 50 GB/s — the machine's
         memory wall; `pairLogLik` end-to-end ≈ 2.5–3× at 8 threads, capped
         by the *sequential* `S()` R call (the scaling floor). End-to-end
         (pre-3c state): 138.7 → 116.1 s (scale 10, 22 sweeps).
      3. [x] **3c — `sampleD` restructure** — done in commits bc4e598 +
         ba61ce8 (2026-10-21): the global `std::map<int, std::set<int>>
         _DeltaMap` (+ `initDeltaMap()`/`Deltafind()`) is replaced by flat
         0-based partner arrays (`linkOfA[i]`/`linkOfB[j]`) initialised from
         a new `links` argument — O(1) linked-test/insert/erase, no global
         state, identical link ordering; the map itself turned out to be
         <1% of the scan (branch-2 pairs are only nlinks per sweep), the
         real cost was per-pair transdendentals. So part 2 hoists the
         loop-invariant logs (gamma is a scalar: `log γ`, `log(1-γ)` once;
         `log(nB - nlinkrec)` maintained across flips) and moves the
         per-pair non-LLL delta to a **parallel precompute pass** (the plan's
         Option C), leaving the serial scan ≈ 1 exp + 1 uniform + branch.
         Delta expressions keep the original association order ⇒ every
         delta/pswitch/loglik value bitwise identical; only the RNG stream
         position changes (uniform block per pair index; distributions
         unchanged — goldens pass unchanged). Measured (scale 10, 22 sweeps,
         threads=1): 136.4 → **95.1 s** (6.20 → 4.32 s/sweep, 1.43×);
         scale 1: 2.60 → **2.27 s**. Profile (13-sweep canonical, threads=1):
         `pairLogLik` **40.7%**, `sampleD` **19.8%** (was 44.6% — 3.1× on
         the scan), S-internal (`as.vector`+`<Anonymous>`+`%*%`) **30.7%**,
         `buildCandidatePairs` 4.2%, `sampleH` 0.9%.
         **Profiling caveat discovered:** Rprof + OpenMP over-attributes
         wall time — at threads=8 Rprof showed the run 1.6× SLOWER than
         threads=1 while the unprofiled wall-clock is 1.3× faster (SIGPROF
         lands in spin-waiting workers). Always profile OpenMP code at
         threads=1 and verify parallel gains by wall-clock only.
      4. [ ] **Re-run the Option A decision rule with the new numbers.**
         `sampleD`'s serial scan is now ~0.84 s/sweep (19.8%) — borderline
         against the 15–20% trigger — but its content is now mostly
         irreducible exp+uniform work, so batching buys at most ~0.6 s/sweep
         at scale 10 for the full exactness machinery. Decide AFTER 3a+
         (next), which shifts every share again.
      5. [ ] **3a+ — compute `S` in C++ for the built-in survival models**
         (exponential/Weibull/Gompertz/piecewise; custom models keep the R
         callback). Now the second-biggest block (30.7%) AND the scaling
         floor of `pairLogLik` under OpenMP: the sequential `as.vector`+
         `%*%`+`exp` on the nPairs×(1+ncov) matrix is thread-independent, so
         8 threads currently buy only ~1.3× end-to-end (95.1 → 73.3 s). A
         C++ `S` (naive dot + vectorised exp, parallel) removes both: expect
         S-internal ~1.3 → ~0.2 s/sweep and pairLogLik 1.73 → ~0.5-0.6 s/sweep
         at 8 threads. BLAS-vs-naive-dot differences ~1e-16, inside the
         golden tolerance.
      6. [ ] **Optional: fuse pairLogLik's per-PIV passes** into one pass
         over pairs (stable terms first, unstable after S returns) — same
         per-q addition order ⇒ bitwise identical; cuts streamed traffic
         ~2–3× on the DRAM-bound stable loops. Only if 3a+ leaves
         `pairLogLik` on top.
