# Plot the distribution of linkage scores

Plot the distribution of linkage scores

## Usage

``` r
plot_linkage_scores(n_pairs, LinkScore)
```

## Arguments

- n_pairs:

  Integer, total number of candidate pairs considered
  (`nrow(A) * nrow(B)`).

- LinkScore:

  Numeric vector, linkage scores of the pairs above 0 (e.g. `Delta$x`).

## Value

`NULL`, invisibly; called for its plotting side effect.

## Examples

``` r
fit_Delta_x <- c(0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.5,0.6,0.7,0.7,0.7,0.7,0.7,0.7,0.8,0.8)
plot_linkage_scores(1000, fit_Delta_x)
```
