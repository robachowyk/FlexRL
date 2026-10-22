# Changelog

## FlexRL 1.0.1 (release candidate)

- Bugs fixed.
- Performance: the Gibbs machinery is ~100x faster end-to-end on the
  benchmark problem (`review/bench-stem.R`, nA=10000, nB=15000,
  n_links=7500, 22 sweeps: 238.7 s -\> 2.29 s); on a 100k/150k-record
  problem (75k links, 22 sweeps) a serial run takes 88 s and
  `threads = 8` brings it to ~51 s. The candidate-pair log-likelihood is
  now a plain per-pair vector (the nA x nB sparse matrix and data-frame
  row subsetting are gone), candidate pairs are built by a single C++
  call (`buildCandidatePairs()`), the per-sweep likelihood pieces are
  computed by one fused C++ pass (`pairLogLik()`), `sampleL` samples the
  joint (H^A, H^B) through an exact O(nval) factorization, and all
  kernels draw from pre-drawn uniform blocks.
- **Breaking (internal kernels):** code calling the samplers directly
  must be adjusted:
  - `sampleD()` signature is now
    `(S, LLA, LLB, LLL, gamma, nlinkrec, links, threads)`. `gamma` is a
    single scalar; the linkage state is passed as `links`, a 0-based k x
    2 matrix of (A, B) record pairs (0x2 when nothing is linked),
    replacing the `sumRowD`/`sumColD` arguments — the flags are derived
    from `links` and still returned. The return list gains `links` (the
    updated 0-based pairs) next to `linksQ` (1-based rows of `S`);
    `loglik` was dropped from the inputs (the value cancels in every
    acceptance ratio and is reconstructed internally).
  - [`simulateD()`](https://robachowyk.github.io/FlexRL/reference/simulateD.md)
    likewise takes the current `linksR` and no `sumRowD`/`sumColD`; it
    gains `threads` and `model_dynamics`.
  - The global C++ link map — and with it `initDeltaMap()` and
    `Deltafind()` — is gone: a sweep carries no state between calls,
    everything it needs lives in its arguments.
  - `sampleNL()` takes a `unif` argument; `sampleL()` takes
    `(survivalpSameH, nval, phikA, phikB, eta, unif)` and returns the
    joint (H^A, H^B) matrix.
  - Kernels draw their Uniform(0,1) deviates from pre-drawn blocks
    consumed by pair position rather than calling the RNG per draw: for
    a given [`set.seed()`](https://rdrr.io/r/base/Random.html) the
    realised draws differ from \<= 1.0.0 (distributions do not).
- OpenMP: the per-pair likelihood loops in the Gibbs kernels are
  parallelised; new `threads` argument on
  [`StEM()`](https://robachowyk.github.io/FlexRL/reference/StEM.md)
  (default 4; on a build without OpenMP the code runs sequentially
  regardless). Results are bitwise identical for any thread count. The
  package still builds — and runs sequentially — on platforms without
  OpenMP.
- The built-in survival models (exponential, Weibull, Gompertz,
  piecewise) evaluate through a fused, parallel C++ kernel
  (`survivalS_cpp`) instead of R-level matrix algebra; custom survival
  functions are unaffected and keep working exactly as before. Values
  agree with the previous R implementation to ~1e-16 (BLAS rounding);
  [`survival_model()`](https://robachowyk.github.io/FlexRL/reference/survival_model.md)
  objects gain a `cuts` field.
- Parity verified against the pre-optimisation implementation and pinned
  by permanent regression tests: `tests/testthat/test-golden-parity.R`
  compares the current kernels to golden fixtures generated from the
  pre-optimisation build (`review/generate-goldens.R`, baseline
  `glm-notes` @ 4b09088) — likelihoods and candidate-pair sets exactly,
  StEM seed clouds and the tiny problem’s stationary frequencies at
  Monte-Carlo level. See `review/PARITY-REPORT.md`.

## FlexRL 1.0.0 (release candidate)

CRAN release: 2026-10-07

- Functions and arguments renamed, data configuration functions added.
- New
  [`survival_model()`](https://robachowyk.github.io/FlexRL/reference/survival_model.md):
  the change over time of a structured PIV can be modelled with an
  exponential (default), Weibull, Gompertz, piecewise-constant or
  user-defined survival function, passed to
  [`StEM()`](https://robachowyk.github.io/FlexRL/reference/StEM.md)
  through the `survival_model` argument. `alpha` chains are named after
  the model parameters. `SurvivalUnstable()` and `loglikSurvival()` are
  obsolete. `cond_hazard_params` are now given on the parameter scale of
  the model (log-hazards for the exponential model).
- New
  [`RL_diagnostics()`](https://robachowyk.github.io/FlexRL/reference/RL_diagnostics.md):
  provides record linkage diagnostics over linkage decision rules
  (thresholds from 0.50 to 0.99); internally computes FDP_measures and
  discrepancy_measures useful for downstream inference;
  `print(x, threshold = )` and `plot(x, type = )` display the
  diagnostics,
- New `link_with_<pkg>()` record linkage wrappers created for several
  packages: multilink, fastLink, BRL, reclin2, diyar, fedmatch, FlexRL.
  They return a named list (`idxA`, `idxB`, `LinkScore`) that can be
  passed directly to
  [`RL_diagnostics()`](https://robachowyk.github.io/FlexRL/reference/RL_diagnostics.md).
- Useless dependencies dropped; R \>= 4.1.0.
- Bugs corrected, warnings added.

## FlexRL 0.1.0 (development)

CRAN release: 2024-12-09

- Package creation. Initial CRAN submission.
