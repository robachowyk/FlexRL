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

Known per-kernel costs that matter at scale:

- `sampleNL` / `sampleL` call `Rcpp::sample()` **once per record / per linked
  pair** — each call crosses into R's `sample()` machinery (SEXP allocation,
  argument checking). At ~1–2 µs per call and (nA+nB)·n_PIVs + n_links·n_PIVs
  calls per sweep this is plausibly the dominant cost.
- `sampleL` allocates 3+ `NumericVector`s of size nval² *per pair* and multiplies
  them through Rcpp sugar (temporaries + GC pressure); the sampling support is an
  nval² grid built by `ExpandGrid` every sweep.
- `simulateD` builds `LLL` as an nA×nB `Matrix::sparseMatrix` and does sparse
  sub-assignment `LLL[select] <- LLL[select] + …` **per PIV**, then extracts
  `LLL[select]` three more times (anyNA check, `log_lik`'s `sum(LLL[links])`,
  pass to `sampleD`). `Matrix` sparse sub-assignment has heavy constants.
- `simulateD` builds the candidate set through a string pipeline:
  `pasteIntoPattern → unique → factor → indexPatterns → pairPatterns`
  (per-record `ostringstream` + `factor()` over millions of strings per sweep).
- `sampleD`'s scan itself is cheap per pair but touches ~6 transcendentals and
  one `runif` per actionable pair, plus `std::map` insert/erase per accepted flip.

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

- [ ] Phase 1.1 — inverse-CDF sampler + `sampleNL` on pre-drawn uniforms
- [ ] Phase 1.2 — factorized `sampleL` (joint (H^A, H^B) return)
- [ ] Phase 1.3 — per-pair `LLLq`, `sampleD` loglik reconstruction + `linksQ`
- [ ] Phase 1.4 — `buildCandidatePairs()` replaces the string pipeline
- [ ] Phase 2 — bench pre/post numbers recorded here
- [ ] Phase 3 — decision per §3 rule
