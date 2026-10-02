# Plot post-linkage diagnostics

Plot post-linkage diagnostics

## Usage

``` r
# S3 method for class 'RL_diagnostics'
plot(x, type, threshold = NULL, ...)
```

## Arguments

- x:

  An `RL_diagnostics` object, see
  [`RL_diagnostics()`](rl_diagnostics.md).

- type:

  One of `"scores"` (linkage score histogram), `"distributions"` (linked
  subset vs. data sources, per variable, at `threshold`),
  `"convergence"` (StEM trace plots, see
  [`plot_StEM_convergence()`](plot_stem_convergence.md)), `"FDP"` (FDP
  estimates by threshold, only if `FDP_estimation = TRUE` was used) or
  `"discrepancy"` (MMD, SMD, IoU and agreement by threshold, see
  [`plot.discrepancy_curves()`](plot.discrepancy_curves.md)).

- threshold:

  Numeric, the linkage decision rule defining the linked subset for
  `type = "distributions"` (ignored for methods without scores).

- ...:

  Passed on to the underlying plotting helper.

## Value

`x`, invisibly; called for its plotting side effect.
