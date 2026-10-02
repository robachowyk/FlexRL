# Run reclin2 through the uniform wrapper interface

Links `dataA` and `dataB` with the package and returns the linked pairs
in the common format used by
[`compute_RL_FDP_score()`](compute_RL_FDP_score.md) and
[`compute_augmRL_FDP_synth()`](compute_augmRL_FDP_synth.md). `dataA` and
`dataB` can be the outputs of [`synthesise()`](synthesise.md) or the
original encoded data sources.

## Usage

``` r
link_with_reclin2(dataA, dataB, arguments, ...)
```

## Arguments

- dataA:

  Data source to link: the encoded source from
  [`prepare_data()`](prepare_data.md), or the augmented sources from
  [`synthesise()`](synthesise.md).

- dataB:

  Data source to link: the encoded source from
  [`prepare_data()`](prepare_data.md), or the augmented sources from
  [`synthesise()`](synthesise.md).

- arguments:

  List of extra arguments forwarded across the `reclin2` pipeline
  (`pair()`, `compare_pairs()`, `problink_em()`,
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  `select_threshold()`): e.g. `on`, `formula`, `type`, `add`,
  `variable`, `score`, `threshold`.

- ...:

  Ignored; kept for compatibility with
  [`RL_diagnostics()`](rl_diagnostics.md).

## Value

Named list with `idxA`, `idxB` and `LinkScore` (row indices in `dataA`,
`dataB` and linkage scores of potential linked records for the
`select_threshold()` argument).

## Details

More details about reclin2 on
https://cran.r-project.org/web/packages/reclin2/index.html
