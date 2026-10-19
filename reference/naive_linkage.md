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
#> 1      1    1
#> 2      2    2
#> 3      3    3
#> 4      3  104
#> 5      4    4
#> 6    124    4
#> 7    131    4
#> 8      5    5
#> 9      7    7
#> 10    11   11
#> 11    12   12
#> 12    13   13
#> 13    14   14
#> 14    15   15
#> 15    16   16
#> 16   148   16
#> 17    18   18
#> 18    18   46
#> 19    46   18
#> 20    46   46
#> 21    19   19
#> 22    21   21
#> 23    22   22
#> 24    23   23
#> 25    23   74
#> 26    24   24
#> 27    25   25
#> 28    26   26
#> 29    27   27
#> 30    28   28
#> 31    29   29
#> 32    30   30
#> 33    31   31
#> 34    32   32
#> 35    32  124
#> 36    32  147
#> 37    33   33
#> 38    34   34
#> 39    35   35
#> 40    36   36
#> 41    37   37
#> 42    40   40
#> 43    41   41
#> 44    42   42
#> 45   112   42
#> 46    43  138
#> 47    44   44
#> 48    44  194
#> 49    45   45
#> 50    45  131
#> 51    48   48
#> 52    49   49
#> 53    50   50
#> 54    50  189
#> 55    51   51
#> 56    53   53
#> 57    58   58
#> 58    59   59
#> 59    59  114
#> 60    59  173
#> 61    63   63
#> 62    65   38
#> 63    66  122
#> 64   116  122
#> 65    67   67
#> 66    78   78
#> 67    79   79
#> 68    80   80
#> 69    81   81
#> 70    83   83
#> 71    83  118
#> 72    84   84
#> 73    84  198
#> 74    87   87
#> 75    87  182
#> 76    90   90
#> 77    92   92
#> 78    94   94
#> 79    95   95
#> 80    97   97
#> 81    99   99
#> 82   111   66
#> 83   120  165
#> 84   121  185
#> 85   123  187
#> 86   130  150
#> 87   144  188
#> 175    6    6
#> 176   64   64
#> 177   93   93
#> 178   80  110
```
