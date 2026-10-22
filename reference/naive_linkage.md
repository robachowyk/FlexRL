# Naive (deterministic exact-match) record linkage

Links records that agree exactly on every non-missing PIV. Does not
enforce the one-to-one assignment constraint, so should be used only to
gauge the difficulty of the linkage task (amount of duplication, and
discriminative power of the PIVs together).

## Usage

``` r
naive_linkage(PIVs, encodedA, encodedB, na_match = TRUE, na_is_zero = TRUE)
```

## Arguments

- PIVs:

  Character vector, names of the PIVs (columns present in both files).

- encodedA:

  The data source (PIVs encoded to natural numbers).

- encodedB:

  The data source (PIVs encoded to natural numbers).

- na_match:

  Logical; if `TRUE`, a missing PIV value is treated as matching any
  value in the other file (default `TRUE`).

- na_is_zero:

  Logical; if `TRUE`, missing values are already coded as `0` (FlexRL's
  convention); if `FALSE`, `NA` will be recoded to `0` (default `TRUE`).

## Value

Data frame with columns `idxA`, `idxB`: pairs of record indices that
match.

## Examples

``` r
PIVs_config <- list( V1 = list(dynamics = "stable",
                               bound_mistakes = c(0.10,0.10),
                               fix_mistakes = c(NA,NA)),
                     V2 = list(dynamics = "stable",
                               bound_mistakes = c(0.10,0.10),
                               fix_mistakes = c(NA,NA)),
                     V3 = list(dynamics = "flexible",
                               bound_mistakes = c(NA,NA),
                               fix_mistakes = c(NA,NA)),
                     V4 = list(dynamics = "structured",
                               bound_mistakes = c(NA,NA),
                               fix_mistakes = c(0.03,0.03),
                               cond_hazard_cov = list(cov1=c("Xe", "Xf"),
                                                      cov2=c())) )
n_values  <- c( 5, 6, 7, 12 )
p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
cond_hazard_params <- list( V1 = c(), V2 = c(), 
                             V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
                           p_mistake, p_missing, cond_hazard_params, TRUE )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> '2' is the larger source, saved as B; '1' saved as A.
naive_linkage( names(prep_data$PIVs_config), 
               prep_data$encodedA, prep_data$encodedB )
#>     idxA idxB
#> 1      2    2
#> 2      3    3
#> 3      4    4
#> 4      5    5
#> 5      6    6
#> 6      7    7
#> 7      7  187
#> 8      9    9
#> 9     10   10
#> 10    11   11
#> 11    11   41
#> 12    11   74
#> 13    41   11
#> 14    41   41
#> 15    41   74
#> 16    74   11
#> 17    74   41
#> 18    74   74
#> 19    12   12
#> 20   100   12
#> 21    13   13
#> 22    14   14
#> 23    15   15
#> 24    16   16
#> 25    17   17
#> 26    18   18
#> 27    19   19
#> 28    21   21
#> 29    22   22
#> 30    24   24
#> 31    25   25
#> 32    26   26
#> 33    29   29
#> 34    30   30
#> 35    31   31
#> 36    32   32
#> 37    33   33
#> 38    34   34
#> 39    35   35
#> 40    37   37
#> 41    39  150
#> 42    43   43
#> 43    44   44
#> 44    47   44
#> 45    45   45
#> 46    46   46
#> 47    50   50
#> 48    52   52
#> 49    54   54
#> 50   122   54
#> 51    55   48
#> 52    57   57
#> 53    58   58
#> 54    61   61
#> 55    65   65
#> 56    67   67
#> 57    85   67
#> 58    69  177
#> 59    70  151
#> 60    75   75
#> 61    80   80
#> 62    82   82
#> 63    83   83
#> 64    86   86
#> 65    87   87
#> 66    89   89
#> 67    90   90
#> 68    90  139
#> 69    91   91
#> 70    93   93
#> 71   133   93
#> 72    96   96
#> 73    97   97
#> 74   109  193
#> 75   120   39
#> 76   129  160
#> 77   142  145
#> 155    1    1
#> 156  124   87
#> 157  124  147
#> 158   59   98
#> 159   98   98
#> 160  108   98
```
