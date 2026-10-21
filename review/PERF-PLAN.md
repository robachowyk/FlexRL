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
