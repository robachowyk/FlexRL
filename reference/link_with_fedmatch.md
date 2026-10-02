# Run fedmatch through the uniform wrapper interface

Links `dataA` and `dataB` with the package and returns the linked pairs
in the common format used by
[`compute_RL_FDP_score()`](compute_RL_FDP_score.md) and
[`compute_augmRL_FDP_synth()`](compute_augmRL_FDP_synth.md). `dataA` and
`dataB` can be the outputs of [`synthesise()`](synthesise.md) or the
original encoded data sources.

## Usage

``` r
link_with_fedmatch(dataA, dataB, arguments, ...)
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

  List of extra arguments forwarded to
  [`fedmatch::merge_plus()`](https://rdrr.io/pkg/fedmatch/man/merge_plus.html)
  (e.g. `by`, `match_type`, `unique_key_1`, `unique_key_2`,
  `multivar_settings`).

- ...:

  Ignored; kept for compatibility with
  [`RL_diagnostics()`](rl_diagnostics.md).

## Value

Named list with `idxA`, `idxB` and `LinkScore` (row indices in `dataA`,
`dataB` and linkage scores of potential linked records).

## Details

More details about fedmatch on
https://cran.r-project.org/web/packages/fedmatch/index.html
