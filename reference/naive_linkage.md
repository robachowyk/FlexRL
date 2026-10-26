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
#> Warning: PIVs 'V4' and 'V3' are associated (Cramer's V = 0.30 in data1, 0.23 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
#> '2' is the larger source, saved as B; '1' saved as A.
naive_linkage( names(prep_data$PIVs_config), 
               prep_data$encodedA, prep_data$encodedB )
#>     idxA idxB
#> 1      1    1
#> 2      2    2
#> 3      3    3
#> 4      3  116
#> 5      4    4
#> 6      5    5
#> 7      6    6
#> 8      7    7
#> 9      8    8
#> 10     9    9
#> 11    10   10
#> 12    12   12
#> 13    13   13
#> 14    15   15
#> 15    16   16
#> 16    17   17
#> 17    18   18
#> 18    21   21
#> 19    24   24
#> 20    25   25
#> 21    28   28
#> 22    29   29
#> 23    32   32
#> 24    33   33
#> 25    33  126
#> 26    34   34
#> 27    35   35
#> 28    36   36
#> 29    40   81
#> 30    81   81
#> 31    41   41
#> 32    42   42
#> 33    44   44
#> 34    45   45
#> 35    46   46
#> 36    46   95
#> 37    95   46
#> 38    95   95
#> 39    47   47
#> 40    48   48
#> 41    49   49
#> 42    59   49
#> 43    50   50
#> 44    54   54
#> 45    58   58
#> 46    65   65
#> 47    66   66
#> 48    69   69
#> 49    72   72
#> 50    77   77
#> 51    79   79
#> 52    80   80
#> 53    84   84
#> 54    90   90
#> 55    93   90
#> 56    97   97
#> 57    99   99
#> 58   100   76
#> 59   103   37
#> 60   114  135
#> 61   129   67
#> 62   129  121
#> 63   133  168
#> 64   136  190
#> 129   78   78
#> 130   78  148
#> 131   78  155
#> 132   68   29
#> 133   68   68
#> 134  146   76
#> 135   92   92
#> 136  111   50
#> 137  111   77
#> 138  111   99
```
