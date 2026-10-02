# Log number of possible linkage configurations

Computes `log(nB! / (nB - sumD)!)`, i.e. the log number of ways to
choose `sumD` ordered links among `n_records_B` records of the larger
file; used in [`log_lik()`](log_lik.md) to normalise the likelihood of
the linkage matrix.

## Usage

``` r
log_possible_config(n_records_B, sumD)
```

## Arguments

- n_records_B:

  Integer, number of records in the larger data source (B).

- sumD:

  Integer, number of currently linked records.

## Value

Numeric, `sum(log((n_records_B - sumD + 1):n_records_B))`, or `0` if
`sumD == 0`.

## Examples

``` r
log_possible_config(n_records_B = 15, sumD = 5)
#> [1] 12.79486
```
