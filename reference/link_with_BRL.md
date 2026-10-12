# Run BRL through the uniform wrapper interface

Links `dataA` and `dataB` with the package and returns the linked pairs
in the common format used by
[`compute_RL_FDP_score()`](https://robachowyk.github.io/FlexRL/reference/compute_RL_FDP_score.md)
and
[`compute_augmRL_FDP_synth()`](https://robachowyk.github.io/FlexRL/reference/compute_augmRL_FDP_synth.md).
`dataA` and `dataB` can be the outputs of
[`synthesise()`](https://robachowyk.github.io/FlexRL/reference/synthesise.md)
or the original encoded data sources.

## Usage

``` r
link_with_BRL(dataA, dataB, arguments, ...)
```

## Arguments

- dataA:

  Data source to link: the encoded source from
  [`prepare_data()`](https://robachowyk.github.io/FlexRL/reference/prepare_data.md),
  or the augmented sources from
  [`synthesise()`](https://robachowyk.github.io/FlexRL/reference/synthesise.md).

- dataB:

  Data source to link: the encoded source from
  [`prepare_data()`](https://robachowyk.github.io/FlexRL/reference/prepare_data.md),
  or the augmented sources from
  [`synthesise()`](https://robachowyk.github.io/FlexRL/reference/synthesise.md).

- arguments:

  List of extra arguments forwarded to
  [`BRL::compareRecords()`](https://rdrr.io/pkg/BRL/man/compareRecords.html)
  and
  [`BRL::bipartiteGibbs()`](https://rdrr.io/pkg/BRL/man/bipartiteGibbs.html)
  (e.g. `flds`, `types`, `nIter`).

- ...:

  Ignored; kept for compatibility with
  [`RL_diagnostics()`](https://robachowyk.github.io/FlexRL/reference/rl_diagnostics.md).

## Value

Named list with `idxA`, `idxB` and `LinkScore` (row indices in `dataA`,
`dataB` and linkage scores of potential linked records).

## Details

More details about BRL on
https://cran.r-project.org/web/packages/BRL/index.html
