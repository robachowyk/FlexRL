# Run multilink through the uniform wrapper interface

Links `dataA` and `dataB` with the package and returns the linked pairs
in the common format used by
[`compute_RL_FDP_score()`](compute_RL_FDP_score.md) and
[`compute_augmRL_FDP_synth()`](compute_augmRL_FDP_synth.md). `dataA` and
`dataB` can be the outputs of [`synthesise()`](synthesise.md) or the
original encoded data sources.

## Usage

``` r
link_with_multilink(dataA, dataB, arguments, ...)
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

  List of extra arguments forwarded across the `multilink` pipeline
  (`create_comparison_data()`, `reduce_comparison_data()`,
  `specify_prior()`, `gibbs_sampler()`, `find_bayes_estimate()`,
  `relabel_bayes_estimate()`): e.g. `records`, `types`, `breaks`,
  `duplicates`, `n_iter`.

- ...:

  Ignored; kept for compatibility with
  [`RL_diagnostics()`](rl_diagnostics.md).

## Value

Named list with `idxA`, `idxB` and `LinkScore = NULL` (row indices of
default potential linked records in `dataA`, `dataB`, this package does
not return scores).

## Details

More details about multilink on
https://cran.r-project.org/web/packages/multilink/index.html
