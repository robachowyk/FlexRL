# FlexRL 1.0.1 (release candidate)

- Bugs fixed.
- Performance: the Gibbs machinery is ~29x faster end-to-end on the benchmark
  problem (`review/bench-stem.R`, nA=10000, nB=15000, n_links=7500, 22 sweeps:
  238.7 s -> 8.3 s). The candidate-pair log-likelihood is now a plain per-pair
  vector (the nA x nB sparse matrix and data-frame row subsetting are gone),
  candidate pairs are built by a single C++ call (`buildCandidatePairs()`),
  `sampleL` samples the joint (H^A, H^B) through an exact O(nval) factorization,
  and all kernels draw from pre-drawn uniform blocks.
- **Breaking (internal kernels):** `sampleNL()` takes a `unif` argument;
  `sampleL()` takes `(survivalpSameH, nval, phikA, phikB, eta, unif)` and
  returns the joint (H^A, H^B) matrix; `sampleD()` no longer takes `loglik`
  (it reconstructs the value, which cancels in every acceptance ratio) and
  returns `linksQ`. Random streams are not comparable to <= 1.0.0 seeds:
  results for a given `set.seed()` differ, distributions do not.

# FlexRL 1.0.0 (release candidate)

- Functions and arguments renamed, data configuration functions added.
- New `survival_model()`: the change over time of a structured PIV can be
  modelled with an exponential (default), Weibull, Gompertz, piecewise-constant
  or user-defined survival function, passed to `StEM()` through the 
  `survival_model` argument. `alpha` chains are named after the model parameters. 
  `SurvivalUnstable()` and `loglikSurvival()` are obsolete. `cond_hazard_params` 
  are now given on the parameter scale of the model (log-hazards for the
  exponential model).
- New `RL_diagnostics()`: provides record linkage diagnostics over linkage 
  decision rules (thresholds from 0.50 to 0.99); internally computes 
  FDP_measures and discrepancy_measures useful for downstream inference;
  `print(x, threshold = )` and `plot(x, type = )` display the diagnostics,
- New `link_with_<pkg>()` record linkage wrappers created for several packages:
  multilink, fastLink, BRL, reclin2, diyar, fedmatch, FlexRL. They return a 
  named list (`idxA`, `idxB`, `LinkScore`) that can be passed directly 
  to `RL_diagnostics()`.
- Useless dependencies dropped; R >= 4.1.0.
- Bugs corrected, warnings added.

# FlexRL 0.1.0 (development)

- Package creation. Initial CRAN submission.
