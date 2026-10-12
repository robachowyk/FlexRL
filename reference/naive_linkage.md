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
#> 2      3   79
#> 3     79   79
#> 4      4   99
#> 5     99   99
#> 6      5    5
#> 7      6    6
#> 8      7    7
#> 9      8    8
#> 10    11   11
#> 11    11   16
#> 12    16   11
#> 13    16   16
#> 14    12   12
#> 15    15   15
#> 16    17   17
#> 17    18   18
#> 18    19   19
#> 19    20   20
#> 20    21   21
#> 21    30   21
#> 22    22   22
#> 23    23   23
#> 24    24   24
#> 25    25   25
#> 26    26   26
#> 27    27   27
#> 28   103   27
#> 29    28   28
#> 30    29   29
#> 31    31   31
#> 32    32   32
#> 33    33   33
#> 34    34   34
#> 35    35   35
#> 36    36   36
#> 37    37   37
#> 38    39   39
#> 39    40   40
#> 40   112   40
#> 41    41   41
#> 42    42   42
#> 43   123   42
#> 44    44   44
#> 45    45   45
#> 46    46   46
#> 47    47   47
#> 48    49   49
#> 49    52   52
#> 50    53   91
#> 51    54   54
#> 52    56   56
#> 53    56   64
#> 54   127   56
#> 55   127   64
#> 56    58   58
#> 57    62   62
#> 58    63   63
#> 59    63   94
#> 60    94   63
#> 61    94   94
#> 62    65   65
#> 63    67   67
#> 64    69   69
#> 65    70   70
#> 66    71   71
#> 67    72   80
#> 68    74   74
#> 69   126   74
#> 70    76   76
#> 71    77   77
#> 72    81  154
#> 73    82  114
#> 74    84   51
#> 75    84   84
#> 76    87   87
#> 77    97  197
#> 78   100  100
#> 79   108  157
#> 80   113  153
#> 81   132  134
#> 163   86   10
#> 164   86   71
#> 165   86   86
#> 166   38   38
#> 167   55  186
#> 168   28    9
#> 169  143   66
```
