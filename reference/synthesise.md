# Augment file B with synthetic records

Fits a generative model on `encodedB[, PIVs]` and draws `n_synth` new
synthetic records from it, appended to `encodedB` with
`source = "synthetic"`; used by
[`compute_augmRL_FDP_synth()`](compute_augmRL_FDP_synth.md) to estimate
the false discovery proportion of a record-linkage method without ground
truth.

## Usage

``` r
synthesise(
  method,
  encodedA,
  encodedB,
  PIVs,
  n_synth,
  restrict_support_intersection = TRUE
)
```

## Arguments

- method:

  One of `"arf"`
  ([`arf::adversarial_rf()`](https://bips-hb.github.io/arf/reference/adversarial_rf.html)),
  `"synthpop"`
  ([`synthpop::syn()`](https://rdrr.io/pkg/synthpop/man/syn.html)), or
  `"mice"`
  ([`mice::mice()`](https://amices.org/mice/reference/mice.html)).

- encodedA:

  The (already prepared / encoded) data source.

- encodedB:

  The (already prepared / encoded) data source.

- PIVs:

  Character vector, names of the PIVs to synthesise.

- n_synth:

  Integer, number of synthetic records to generate.

- restrict_support_intersection:

  Logical; drop synthetic records whose PIV values fall outside the
  support shared with `encodedA` (default `TRUE`).

## Value

List with `dataA` (unchanged `encodedA`, support-restricted) and `dataB`
(`encodedB` plus the synthetic records).
