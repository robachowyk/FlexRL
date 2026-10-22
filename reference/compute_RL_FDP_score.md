# Estimate the false discovery proportion of a record-linkage method via model specific linkage scores

Estimate the false discovery proportion of a record-linkage method via
model specific linkage scores

## Usage

``` r
compute_RL_FDP_score(
  encodedA,
  encodedB,
  PIVs,
  maxIter4CV = 10,
  n_repeats = 10,
  RL_method,
  ...
)
```

## Arguments

- encodedA:

  Encoded data source as returned by
  [`prepare_data()`](https://robachowyk.github.io/FlexRL/reference/prepare_data.md)
  (`encodedA` is the smaller one).

- encodedB:

  Encoded data source as returned by
  [`prepare_data()`](https://robachowyk.github.io/FlexRL/reference/prepare_data.md)
  (`encodedB` is the larger one).

- PIVs:

  Character vector, names of the PIVs.

- maxIter4CV:

  Integer, max number of retries per iteration if no valid FDP estimate
  is obtained.

- n_repeats:

  Integer, number of augmentation iterations to average over.

- RL_method:

  One of `"multilink"`, `"fastLink"`, `"BRL"`, `"reclin2"`, `"diyar"`,
  `"fedmatch"`, `"FlexRL"`.

- ...:

  Extra arguments forwarded to the chosen link_with\_\* wrapper (i.e. to
  the underlying record-linkage package).

## Value

List with `FDP_score_estimator`, `Linked_pairs`: data frames
(`n_repeats` rows x 50 thresholds, 0.50 to 0.99).

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
#> Warning: PIV 'V4': values out of the common support in data2 [1].
#> Warning: PIVs 'V4' and 'V3' are associated (Cramer's V = 0.32 in data1, 0.19 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)                          
compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
                      "BRL", flds = PIVs, types = rep("bi",length(PIVs)) )
#> BRL results (average over 2 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.19    0.19    0.19    0.19    0.18    0.17    0.15    0.13    0.04    0
#> Linked pairs (RL)            58  57  57  56  54  48  34  18  1   0
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1922989 0.1922989 0.1922989 0.1922989 0.1922989 0.1877583 0.1877583
#> 2 0.1945785 0.1892398 0.1892398 0.1892398 0.1892398 0.1892398 0.1892398
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1877583 0.1877583 0.1877583 0.1877583 0.1877583 0.1877583 0.1877583
#> 2 0.1892398 0.1892398 0.1892398 0.1892398 0.1892398 0.1892398 0.1859921
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1877583 0.1877583 0.1850397 0.1823838 0.1823838 0.1798971 0.1774843
#> 2 0.1859921 0.1859921 0.1859921 0.1859921 0.1859921 0.1859921 0.1859921
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.1774843 0.1774843 0.1735512 0.1699773 0.1664066 0.1629136 0.1629136
#> 2 0.1819342 0.1819342 0.1819342 0.1769935 0.1754000 0.1740136 0.1713712
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.1615404 0.1564167 0.1479085 0.1395785 0.1362140 0.1331111 0.1273545
#> 2 0.1661240 0.1610541 0.1546732 0.1508602 0.1425778 0.1366138 0.1321605
#>        0.85      0.86      0.87      0.88       0.89       0.90       0.91 0.92
#> 1 0.1235185 0.1145370 0.1129293 0.1029630 0.09555556 0.08722222 0.08722222    0
#> 2 0.1305882 0.1266667 0.1166667 0.1133333 0.10666667 0.00000000 0.00000000    0
#>   0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   58   58   58   58   58   57   57   57   57   57   57   57   57   57   57
#> 2   58   57   57   57   57   57   57   57   57   57   57   57   57   56   56
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   57   56   55   55   54   53   53   53   51   49   47   45   45   44   40
#> 2   56   56   56   56   56   56   54   54   54   51   50   49   47   43   39
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   34   29   27   25   21   18   12   11    6    4    2    2    0    0    0
#> 2   34   31   25   21   18   17   13    5    3    1    0    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
                      "FlexRL", data = prep_data, StEM_iter = 5, 
                      StEM_burnin = 2, gibbs_iter = 5, gibbs_burnin = 2,
                      n_post_sample = 5 )
#> FlexRL
#> Warning: PIV 'V1', file A hit its mistake bound on 40% of StEM iterations.
#>   PIV 'V1', file B hit its mistake bound on 40% of StEM iterations.
#> FlexRL
#> FlexRL results (average over 2 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.14    0.14    0.05    0.05    0.05    0.05    0   0   0   0
#> Linked pairs (RL)            80  80  60  60  60  60  44  44  44  44
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1488372 0.1488372 0.1488372 0.1488372 0.1488372 0.1488372 0.1488372
#> 2 0.1260274 0.1260274 0.1260274 0.1260274 0.1260274 0.1260274 0.1260274
#>        0.57      0.58      0.59       0.60       0.61       0.62       0.63
#> 1 0.1488372 0.1488372 0.1488372 0.05714286 0.05714286 0.05714286 0.05714286
#> 2 0.1260274 0.1260274 0.1260274 0.04912281 0.04912281 0.04912281 0.04912281
#>         0.64       0.65       0.66       0.67       0.68       0.69       0.70
#> 1 0.05714286 0.05714286 0.05714286 0.05714286 0.05714286 0.05714286 0.05714286
#> 2 0.04912281 0.04912281 0.04912281 0.04912281 0.04912281 0.04912281 0.04912281
#>         0.71       0.72       0.73       0.74       0.75       0.76       0.77
#> 1 0.05714286 0.05714286 0.05714286 0.05714286 0.05714286 0.05714286 0.05714286
#> 2 0.04912281 0.04912281 0.04912281 0.04912281 0.04912281 0.04912281 0.04912281
#>         0.78       0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> 1 0.05714286 0.05714286    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.04912281 0.04912281    0    0    0    0    0    0    0    0    0    0    0
#>   0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   86   86   86   86   86   86   86   86   86   86   63   63   63   63   63
#> 2   73   73   73   73   73   73   73   73   73   73   57   57   57   57   57
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   63   63   63   63   63   63   63   63   63   63   63   63   63   63   63
#> 2   57   57   57   57   57   57   57   57   57   57   57   57   57   57   57
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   45   45   45   45   45   45   45   45   45   45   45   45   45   45   45
#> 2   43   43   43   43   43   43   43   43   43   43   43   43   43   43   43
#>   0.95 0.96 0.97 0.98 0.99
#> 1   45   45   45   45   45
#> 2   43   43   43   43   43
#> 
```
