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
gen_data <- simulate_data( PIVs_config, n_values, c(250, 300), 200, 
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
#> 4      5    5
#> 5     63    5
#> 6      8    8
#> 7      8  163
#> 8    245    8
#> 9    245  163
#> 10     9    9
#> 11     9  197
#> 12   197    9
#> 13   197  197
#> 14    11   11
#> 15    11  206
#> 16    12   12
#> 17    13   13
#> 18    14   14
#> 19    14   67
#> 20    67   14
#> 21    67   67
#> 22    15   15
#> 23    15  150
#> 24    16   90
#> 25    90   90
#> 26    17   17
#> 27    18   18
#> 28    18  181
#> 29   181   18
#> 30   181  181
#> 31    19   19
#> 32    20   20
#> 33    21   21
#> 34    22   22
#> 35    71   22
#> 36   117   22
#> 37    24   24
#> 38    25   25
#> 39    26   26
#> 40    27   27
#> 41    28   28
#> 42    30   30
#> 43    30   72
#> 44    72   30
#> 45    72   72
#> 46    33   33
#> 47    34   34
#> 48    36   36
#> 49    36  151
#> 50   151   36
#> 51   151  151
#> 52    37   37
#> 53    38   38
#> 54    39   39
#> 55    40   40
#> 56    41   41
#> 57    42   42
#> 58    43   43
#> 59    44   44
#> 60    45   45
#> 61    46   46
#> 62    46  169
#> 63   119   46
#> 64   119  169
#> 65    47   47
#> 66    49   49
#> 67    50   50
#> 68    51   51
#> 69    52   52
#> 70    53   53
#> 71    54   54
#> 72    55  133
#> 73    55  235
#> 74    56   56
#> 75    57   57
#> 76   196   57
#> 77    58   58
#> 78    59   59
#> 79    59  175
#> 80   175   59
#> 81   175  175
#> 82   220   59
#> 83   220  175
#> 84    60   60
#> 85    60  202
#> 86    64   64
#> 87    64  109
#> 88    65   65
#> 89    65  203
#> 90    66   66
#> 91    69   69
#> 92    70   70
#> 93    73   73
#> 94    74   74
#> 95    74   91
#> 96    91   74
#> 97    91   91
#> 98   213   74
#> 99   213   91
#> 100   75   75
#> 101   75  153
#> 102   78   78
#> 103   79  156
#> 104   80   80
#> 105   81   81
#> 106   83   83
#> 107  210   83
#> 108   84   84
#> 109   85   85
#> 110   85  216
#> 111   86  255
#> 112   87   87
#> 113   89   89
#> 114   92   92
#> 115   93   93
#> 116   94   94
#> 117   95  293
#> 118   96   96
#> 119   99   99
#> 120  102  102
#> 121  104  104
#> 122  105   76
#> 123  172   76
#> 124  107  107
#> 125  110  269
#> 126  111  111
#> 127  112  112
#> 128  112  271
#> 129  113  113
#> 130  114  114
#> 131  115  115
#> 132  238  115
#> 133  116  116
#> 134  118  118
#> 135  123  248
#> 136  124  124
#> 137  124  178
#> 138  125  180
#> 139  126  126
#> 140  127  252
#> 141  128  128
#> 142  129  129
#> 143  130  130
#> 144  130  221
#> 145  131  131
#> 146  131  187
#> 147  187  131
#> 148  187  187
#> 149  132  132
#> 150  137  137
#> 151  143  134
#> 152  144  144
#> 153  145  145
#> 154  146  146
#> 155  154  154
#> 156  154  186
#> 157  186  154
#> 158  186  186
#> 159  157  157
#> 160  158  158
#> 161  161  161
#> 162  161  185
#> 163  162  162
#> 164  165  165
#> 165  167  165
#> 166  166  166
#> 167  166  220
#> 168  168  168
#> 169  170  170
#> 170  171   97
#> 171  171  171
#> 172  174  276
#> 173  209  276
#> 174  176  176
#> 175  177  177
#> 176  178  193
#> 177  178  280
#> 178  193  193
#> 179  193  280
#> 180  179  179
#> 181  182  182
#> 182  183  183
#> 183  184  184
#> 184  184  277
#> 185  195  195
#> 186  198  198
#> 187  218  225
#> 188  221  289
#> 189  237  204
#> 190  248  188
#> 381   82  194
#> 382  194  194
#> 383  124  232
#> 384  247  232
#> 385  135  135
#> 386  121  121
#> 387  121  276
#> 388  190   62
#> 389  190   71
#> 390  190  190
#> 391  246   39
#> 392  246   53
#> 393  246  149
#> 394  246  207
#> 395   75  274
#> 396  225    9
#> 397  225  195
#> 398  225  197
#> 399  225  252
#> 400   49  237
#> 401  178  237
#> 402  193  237
#> 403  226  237
#> 404   77   77
#> 405  159   42
#> 406  159  154
#> 407  159  159
#> 408  159  186
#> 409  204   32
#> 410  204  123
#> 411  204  248
#> 412  229  157
#> 413  229  295
#> 414   59  217
#> 415  175  217
#> 416  220  217
#> 417  233  217
#> 418   64  240
#> 419  109  240
#> 420  164  240
#> 421   87  297
#> 422   95  297
#> 423  106  297
```
