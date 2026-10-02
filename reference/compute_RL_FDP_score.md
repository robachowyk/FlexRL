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

  Encoded data source as returned by [`prepare_data()`](prepare_data.md)
  (`encodedA` is the smaller one).

- encodedB:

  Encoded data source as returned by [`prepare_data()`](prepare_data.md)
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
gen_data <- simulate_data( PIVs_config, n_values, c(250, 300), 200, 
                           p_mistake, p_missing, cond_hazard_params, TRUE )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)                          
compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
                      "BRL", flds = PIVs, types = rep("bi",length(PIVs)) )
#> BRL results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.15   0.15   0.15   0.15   0.15   0.15   0.15
#> Linked pairs (RL)         107.00 107.00 107.00 107.00 106.00 106.00 106.00
#>                             0.57   0.58   0.59   0.60   0.61   0.62   0.63
#> FDP model score estimator   0.15   0.14   0.14   0.14   0.14   0.14   0.14
#> Linked pairs (RL)         105.00 104.00 104.00 104.00 104.00 104.00 104.00
#>                             0.64   0.65   0.66   0.67   0.68   0.69   0.70
#> FDP model score estimator   0.14   0.14   0.14   0.14   0.14   0.14   0.14
#> Linked pairs (RL)         104.00 103.00 103.00 103.00 103.00 103.00 103.00
#>                             0.71   0.72   0.73  0.74  0.75  0.76  0.77  0.78
#> FDP model score estimator   0.14   0.14   0.14  0.14  0.14  0.14  0.14  0.13
#> Linked pairs (RL)         103.00 102.00 100.00 99.00 99.00 98.00 97.00 96.00
#>                            0.79  0.80  0.81  0.82  0.83  0.84  0.85 0.86 0.87
#> FDP model score estimator  0.13  0.13  0.13  0.12  0.12  0.11  0.11  0.1  0.1
#> Linked pairs (RL)         95.00 92.00 91.00 80.00 76.00 65.00 56.00 51.0 46.0
#>                           0.88  0.89  0.90  0.91 0.92 0.93 0.94 0.95 0.96 0.97
#> FDP model score estimator  0.1  0.09  0.08  0.08 0.07 0.07    0    0    0    0
#> Linked pairs (RL)         41.0 31.00 21.00 15.00 5.00 2.00    0    0    0    0
#>                           0.98 0.99
#> FDP model score estimator    0    0
#> Linked pairs (RL)            0    0
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1530737 0.1530737 0.1530737 0.1530737 0.1500943 0.1500943 0.1500943
#> 2 0.1530737 0.1530737 0.1530737 0.1530737 0.1500943 0.1500943 0.1500943
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1473545 0.1447222 0.1447222 0.1447222 0.1447222 0.1447222 0.1447222
#> 2 0.1473545 0.1447222 0.1447222 0.1447222 0.1447222 0.1447222 0.1447222
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1447222 0.1427184 0.1427184 0.1427184 0.1427184 0.1427184 0.1427184
#> 2 0.1447222 0.1427184 0.1427184 0.1427184 0.1427184 0.1427184 0.1427184
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.1427184 0.1413617 0.1387444 0.1375084 0.1375084 0.1363832 0.1353837
#> 2 0.1427184 0.1413617 0.1387444 0.1375084 0.1375084 0.1363832 0.1353837
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.1345023 0.1336959 0.1313647 0.1307204 0.1234306 0.1207456 0.1132991
#> 2 0.1345023 0.1336959 0.1313647 0.1307204 0.1234306 0.1207456 0.1132991
#>        0.85     0.86      0.87       0.88       0.89       0.90       0.91
#> 1 0.1064286 0.102658 0.0992029 0.09650407 0.09032258 0.08391534 0.08059259
#> 2 0.1064286 0.102658 0.0992029 0.09650407 0.09032258 0.08391534 0.08059259
#>         0.92       0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1 0.07311111 0.06666667    0    0    0    0    0    0
#> 2 0.07311111 0.06666667    0    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1  107  107  107  107  106  106  106  105  104  104  104  104  104  104  104
#> 2  107  107  107  107  106  106  106  105  104  104  104  104  104  104  104
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1  103  103  103  103  103  103  103  102  100   99   99   98   97   96   95
#> 2  103  103  103  103  103  103  103  102  100   99   99   98   97   96   95
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   92   91   80   76   65   56   51   46   41   31   21   15    5    2    0
#> 2   92   91   80   76   65   56   51   46   41   31   21   15    5    2    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
                      "FlexRL", data = prep_data, StEM_iter = 5, 
                      StEM_burnin = 2, gibbs_iter = 5, gibbs_burnin = 2,
                      n_post_samp = 5 )
#> FlexRL
#> FlexRL
#> FlexRL results (average over 2 iterations):
#>                            0.50  0.51  0.52  0.53  0.54  0.55  0.56  0.57  0.58
#> FDP model score estimator   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1
#> Linked pairs (RL)         169.0 169.0 169.0 169.0 169.0 169.0 169.0 169.0 169.0
#>                            0.59   0.60   0.61   0.62   0.63   0.64   0.65
#> FDP model score estimator   0.1   0.05   0.05   0.05   0.05   0.05   0.05
#> Linked pairs (RL)         169.0 141.00 141.00 141.00 141.00 141.00 141.00
#>                             0.66   0.67   0.68   0.69   0.70   0.71   0.72
#> FDP model score estimator   0.05   0.05   0.05   0.05   0.05   0.05   0.05
#> Linked pairs (RL)         141.00 141.00 141.00 141.00 141.00 141.00 141.00
#>                             0.73   0.74   0.75   0.76   0.77   0.78   0.79 0.80
#> FDP model score estimator   0.05   0.05   0.05   0.05   0.05   0.05   0.05    0
#> Linked pairs (RL)         141.00 141.00 141.00 141.00 141.00 141.00 141.00  109
#>                           0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)          109  109  109  109  109  109  109  109  109  109
#>                           0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> FDP model score estimator    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)          109  109  109  109  109  109  109  109  109
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1000000 0.1000000 0.1000000 0.1000000 0.1000000 0.1000000 0.1000000
#> 2 0.1081395 0.1081395 0.1081395 0.1081395 0.1081395 0.1081395 0.1081395
#>        0.57      0.58      0.59       0.60       0.61       0.62       0.63
#> 1 0.1000000 0.1000000 0.1000000 0.04428571 0.04428571 0.04428571 0.04428571
#> 2 0.1081395 0.1081395 0.1081395 0.04647887 0.04647887 0.04647887 0.04647887
#>         0.64       0.65       0.66       0.67       0.68       0.69       0.70
#> 1 0.04428571 0.04428571 0.04428571 0.04428571 0.04428571 0.04428571 0.04428571
#> 2 0.04647887 0.04647887 0.04647887 0.04647887 0.04647887 0.04647887 0.04647887
#>         0.71       0.72       0.73       0.74       0.75       0.76       0.77
#> 1 0.04428571 0.04428571 0.04428571 0.04428571 0.04428571 0.04428571 0.04428571
#> 2 0.04647887 0.04647887 0.04647887 0.04647887 0.04647887 0.04647887 0.04647887
#>         0.78       0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> 1 0.04428571 0.04428571    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.04647887 0.04647887    0    0    0    0    0    0    0    0    0    0    0
#>   0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1  166  166  166  166  166  166  166  166  166  166  140  140  140  140  140
#> 2  172  172  172  172  172  172  172  172  172  172  142  142  142  142  142
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1  140  140  140  140  140  140  140  140  140  140  140  140  140  140  140
#> 2  142  142  142  142  142  142  142  142  142  142  142  142  142  142  142
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1  109  109  109  109  109  109  109  109  109  109  109  109  109  109  109
#> 2  109  109  109  109  109  109  109  109  109  109  109  109  109  109  109
#>   0.95 0.96 0.97 0.98 0.99
#> 1  109  109  109  109  109
#> 2  109  109  109  109  109
#> 
```
