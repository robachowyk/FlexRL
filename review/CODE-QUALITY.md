# FlexRL review — CODE QUALITY NOTES

Observations from a full read of `R/FlexRecordLinkage.r` (~3k lines),
`src/FlexRLfunctions.cpp` (~700 lines), and the test suite. Tone: colleague
notes for the author; none of these are correctness bugs (see FINDINGS.md for
those).

## Structure

1. **Single-file layout.** Everything R-side lives in one 3k-line file:
   simulation, data preparation, StEM, FDP estimators, six package wrappers,
   diagnostics. Splitting into `simulate-data.r`, `prepare-data.r`, `stem.r`,
   `fdp.r`, `wrappers-*.r` would make navigation and code review much easier,
   and would let the wrappers' roxygen examples stay close to their code.
2. **C++ file is comparatively tidy** — kernels are small and documented with
   the corresponding paper equations, which made the exact-posterior review
   feasible. Keep that property when adding kernels.

## API and conventions

3. **Wrapper argument forwarding** —
   `arguments[intersect(names(arguments), names(formals(f::...)))]` is elegant
   and keeps wrappers future-proof, but it **silently drops mistyped argument
   names** (e.g. `niter` vs `nIter` changes behaviour with no signal). A
   one-line check — warn about `setdiff(names(arguments), names(formals(...)))`
   — would turn silent misconfiguration into a visible hint.
4. **Naming style is consistent enough** (`snake_case` functions, `PIVs_config`
   style parameters); the internal `pivsA/pivsB` vs `PIVs` vs `PIVs_stable`
   trio in StEM takes a moment to untangle but is fine.
5. **`same_mistakes` semantics** — default `TRUE` silently requires identical
   per-file mistake rates, and `prepare_data` rejects otherwise. Now covered by
   a test; consider a doc sentence on the StEM page since the constraint spans
   two functions.

## Validation and error messages

6. **Big improvement this review.** Config validation (`PIVs_config` shape,
   dynamics/mistakes coherence, cond_hazard_cov presence), the encoded-data
   contract check in StEM (#K), and the tightened `is_numeric_or_na` (#G) now
   cover all 15 deranged-probe inputs with clear messages.
7. **Remaining unguarded spots** (low priority, documented in FINDINGS §2):
   `nlminb` convergence flag unchecked; `_DeltaMap` persistence between
   externally driven chains; `FDP_synth`'s silent degenerate output when the
   synthetic path is handicapped (`date = 0`) — a warning there would save
   users from trusting a 0.000 FDP.

## C++ specifics

8. **Non-const reference parameters**: `sampleD` was the only kernel taking
   `LogicalVector&` and writing through it (#I-a). It now clones at entry.
   Consider making all kernel parameters const references (or by-value) so this
   class of bug cannot recur; Rcpp makes the intent hard to see at call sites.
9. **Sparse constructions**: always build general matrices via
   `Matrix::sparseMatrix(..., x = numeric(0))` rather than
   `Matrix(0, n, n, sparse = TRUE)` (#I-b). Worth a short comment convention in
   the source; both sites are fixed and commented now.
10. **Integer arithmetic on registered values**: the "0 = missing" sentinel is
    used consistently, but nothing in the C++ defends against `NA_INTEGER`
    (which overflows in `G(i) - 1`). The R-side contract check (#K) now gates
    this; if C++ should ever be called directly, a defensive
    `if (G(i) == NA_INTEGER) G(i) = 0;` at kernel entry would be cheap.

## Documentation

11. Roxygen is thorough and the examples run — genuinely good shape for CRAN.
    Two small gaps found: `proba_same_H` is `NA` for stable-only configs
    (undocumented contract, now tested); the fedmatch wrapper's returned
    indices are now documented as 1-based row positions (behaviour fixed).
12. `vignettes/` ships `OLDFlexRL-*` vignettes (duplicate index entries →
    check NOTE), `.orig` and `.html` artifacts. Moving the OLD pair out of
    `vignettes/` (e.g. `inst/doc/` or a `dev/` folder) clears the NOTE.

## Tests

13. 14 files / 446 expectations, ~66 s — CRAN-friendly. The exact-posterior
    test (~45 s) is the most valuable asset in the suite; if CRAN time limits
    ever bite, split it into a long-test skipped on CRAN rather than removing it.
14. Per FINDINGS §5: never trust exit-code-based suite runners; count
    `failed`/`error` from the result objects. The old harness masked six files
    with real failures.
15. Suggested additions (cheap, high value): a `synthesise`→`StEM` end-to-end
    round trip on the tiny problem; an FDP estimator test with a *known* FDP
    (the sims do this manually; a 10-line test would pin #C permanently).

## Reproducibility hygiene

16. `set.seed()` calls inside exported functions would be inappropriate, so the
    package correctly leaves RNG state to the user — but `simulate_data` and
    `StEM` consume different amounts of RNG per call, so results are only
    reproducible with an explicit seed immediately before each call. A "Reproducibility" note in the README/vignette would help.
