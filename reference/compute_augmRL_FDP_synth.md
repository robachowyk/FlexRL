# Estimate the false discovery proportion of a record-linkage method via synthetic augmentation

Repeatedly augments file B with synthetic records
([`synthesise()`](synthesise.md)), runs the chosen record-linkage method
link_with\_\*, and compares the synthetic-vs-real proportion among
linked pairs to estimate the false discovery proportion, for a range of
score thresholds (0.50 to 0.99). See https://doi.org/10.1002/sim.70292
for the method.

## Usage

``` r
compute_augmRL_FDP_synth(
  synth_method,
  encodedA,
  encodedB,
  PIVs,
  n_synth = NULL,
  restrict_support_intersection = TRUE,
  maxIter4CV = 10,
  n_repeats = 10,
  RL_method,
  ...
)
```

## Arguments

- synth_method:

  Passed to [`synthesise()`](synthesise.md): `"arf"`, `"synthpop"`, or
  `"mice"`.

- encodedA, encodedB:

  The two encoded data sources as returned by
  [`prepare_data()`](prepare_data.md) (`encodedA` is the smaller one).

- PIVs:

  Character vector, names of the PIVs.

- n_synth:

  Integer, number of synthetic records to generate per iteration
  (default: 10% of `nrow(encodedB)`).

- restrict_support_intersection:

  Passed to [`synthesise()`](synthesise.md).

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

List with `FDP_score_estimator`, `FDP_synth_estimator`,
`Linked_pairs_augm`, `Linked_pairs`: data frames (`n_repeats` rows x 50
thresholds, 0.50 to 0.99).

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
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "BRL", flds = PIVs, 
                          types = rep("bi",length(PIVs)) )
#> Iteration: 0, Accuracy: 50.68%
#> Iteration: 1, Accuracy: 40.07%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> Iteration: 0, Accuracy: 49.32%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> BRL results (average over 2 iterations):
#>                            0.50  0.51  0.52  0.53  0.54  0.55  0.56  0.57  0.58
#> FDP model score estimator  0.16  0.16  0.16  0.16  0.16  0.16  0.15  0.15  0.15
#> FDP synth data estimator   0.16  0.10  0.10  0.10  0.10  0.11  0.11  0.11  0.11
#> Linked pairs (augm. RL)   96.00 96.00 96.00 96.00 96.00 95.00 94.00 94.00 94.00
#> Linked pairs (RL)         98.00 98.00 96.00 96.00 96.00 96.00 96.00 96.00 96.00
#>                            0.59  0.60  0.61  0.62  0.63  0.64  0.65  0.66  0.67
#> FDP model score estimator  0.15  0.15  0.15  0.15  0.15  0.15  0.15  0.15  0.15
#> FDP synth data estimator   0.11  0.11  0.11  0.11  0.11  0.11  0.11  0.11  0.11
#> Linked pairs (augm. RL)   94.00 94.00 94.00 94.00 94.00 94.00 94.00 94.00 94.00
#> Linked pairs (RL)         96.00 96.00 96.00 96.00 96.00 96.00 96.00 96.00 96.00
#>                            0.68  0.69  0.70  0.71  0.72  0.73  0.74  0.75  0.76
#> FDP model score estimator  0.15  0.15  0.15  0.15  0.15  0.15  0.15  0.15  0.15
#> FDP synth data estimator   0.11  0.11  0.11  0.11  0.11  0.11  0.11  0.11  0.11
#> Linked pairs (augm. RL)   94.00 94.00 94.00 94.00 94.00 94.00 94.00 94.00 93.00
#> Linked pairs (RL)         96.00 96.00 96.00 96.00 96.00 96.00 95.00 94.00 94.00
#>                            0.77  0.78  0.79  0.80  0.81  0.82  0.83  0.84  0.85
#> FDP model score estimator  0.15  0.15  0.15  0.14  0.14  0.14  0.14  0.13  0.13
#> FDP synth data estimator   0.11  0.11  0.11  0.12  0.13  0.14  0.15  0.00  0.00
#> Linked pairs (augm. RL)   92.00 92.00 89.00 82.00 78.00 72.00 66.00 60.00 48.00
#> Linked pairs (RL)         94.00 92.00 90.00 84.00 78.00 72.00 68.00 60.00 48.00
#>                            0.86  0.87  0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95
#> FDP model score estimator  0.12  0.11  0.11  0.1 0.09 0.04    0    0    0    0
#> FDP synth data estimator   0.00  0.00  0.00  0.0 0.00 0.00    0    0    0    0
#> Linked pairs (augm. RL)   35.00 24.00 16.00  8.0 3.00 0.00    0    0    0    0
#> Linked pairs (RL)         35.00 24.00 16.00  8.0 3.00 0.00    0    0    0    0
#>                           0.96 0.97 0.98 0.99
#> FDP model score estimator    0    0    0    0
#> FDP synth data estimator     0    0    0    0
#> Linked pairs (augm. RL)      0    0    0    0
#> Linked pairs (RL)            0    0    0    0
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1620635 0.1620635 0.1553241 0.1553241 0.1553241 0.1553241 0.1553241
#> 2 0.1620862 0.1586140 0.1586140 0.1586140 0.1586140 0.1555556 0.1524795
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1553241 0.1553241 0.1553241 0.1553241 0.1553241 0.1553241 0.1553241
#> 2 0.1524795 0.1524795 0.1524795 0.1524795 0.1524795 0.1524795 0.1524795
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1553241 0.1553241 0.1553241 0.1553241 0.1553241 0.1553241 0.1553241
#> 2 0.1524795 0.1524795 0.1524795 0.1524795 0.1524795 0.1524795 0.1524795
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.1553241 0.1553241 0.1553241 0.1553241 0.1542339 0.1532742 0.1532742
#> 2 0.1524795 0.1524795 0.1524795 0.1512530 0.1512530 0.1512530 0.1503345
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.1525448 0.1512088 0.1468122 0.1436850 0.1408559 0.1371973 0.1335593
#> 2 0.1495531 0.1473658 0.1431861 0.1398433 0.1354773 0.1337418 0.1303097
#>        0.85      0.86      0.87      0.88      0.89       0.90       0.91 0.92
#> 1 0.1276329 0.1178175 0.1132323 0.1090972 0.1005556 0.09500000 0.00000000    0
#> 2 0.1250000 0.1212169 0.1131687 0.1056944 0.1002222 0.09222222 0.08777778    0
#>   0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0
#> 
#> $FDP_synth_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1030928 0.1030928 0.1052632 0.1052632 0.1052632 0.1052632 0.1052632
#> 2 0.2083333 0.1041667 0.1041667 0.1041667 0.1041667 0.1052632 0.1063830
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1052632 0.1052632 0.1052632 0.1052632 0.1052632 0.1052632 0.1052632
#> 2 0.1063830 0.1063830 0.1063830 0.1063830 0.1063830 0.1063830 0.1063830
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1052632 0.1052632 0.1052632 0.1052632 0.1052632 0.1052632 0.1052632
#> 2 0.1063830 0.1063830 0.1063830 0.1063830 0.1063830 0.1063830 0.1063830
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.1052632 0.1052632 0.1052632 0.1052632 0.1063830 0.1075269 0.1075269
#> 2 0.1063830 0.1063830 0.1063830 0.1075269 0.1075269 0.1075269 0.1086957
#>        0.78      0.79      0.80      0.81      0.82      0.83 0.84 0.85 0.86
#> 1 0.1086957 0.1111111 0.1204819 0.1282051 0.1369863 0.1515152    0    0    0
#> 2 0.1098901 0.1136364 0.1219512 0.1298701 0.1428571 0.1492537    0    0    0
#>   0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   97   97   95   95   95   95   95   95   95   95   95   95   95   95   95
#> 2   96   96   96   96   96   95   94   94   94   94   94   94   94   94   94
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   95   95   95   95   95   95   95   95   95   95   94   93   93   92   90
#> 2   94   94   94   94   94   94   94   94   94   93   93   93   92   91   88
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   83   78   73   66   59   46   28   22   16    6    2    0    0    0    0
#> 2   82   77   70   67   61   50   42   27   16   10    4    1    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   98   98   96   96   96   96   96   96   96   96   96   96   96   96   96
#> 2   98   97   97   97   97   96   95   95   95   95   95   95   95   95   95
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   96   96   96   96   96   96   96   96   96   96   95   94   94   93   91
#> 2   95   95   95   95   95   95   95   95   95   94   94   94   93   92   89
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   84   79   74   67   59   46   28   22   16    6    2    0    0    0    0
#> 2   83   78   71   68   61   50   42   27   16   10    4    1    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "FlexRL", data = prep_data, 
                          StEM_iter = 5, StEM_burnin = 2, 
                          gibbs_iter = 5, gibbs_burnin = 2, n_post_samp = 5 )
#> Iteration: 0, Accuracy: 49.32%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> Iteration: 0, Accuracy: 47.14%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> FlexRL results (average over 2 iterations):
#>                             0.50   0.51   0.52   0.53   0.54   0.55   0.56
#> FDP model score estimator   0.11   0.11   0.11   0.11   0.11   0.11   0.11
#> FDP synth data estimator    0.31   0.31   0.31   0.31   0.31   0.31   0.31
#> Linked pairs (augm. RL)   175.00 175.00 175.00 175.00 175.00 175.00 175.00
#> Linked pairs (RL)         180.00 180.00 180.00 180.00 180.00 180.00 180.00
#>                             0.57   0.58   0.59   0.60   0.61   0.62   0.63
#> FDP model score estimator   0.11   0.11   0.11   0.05   0.05   0.05   0.05
#> FDP synth data estimator    0.31   0.31   0.31   0.34   0.34   0.34   0.34
#> Linked pairs (augm. RL)   175.00 175.00 175.00 144.00 144.00 144.00 144.00
#> Linked pairs (RL)         180.00 180.00 180.00 148.00 148.00 148.00 148.00
#>                             0.64   0.65   0.66   0.67   0.68   0.69   0.70
#> FDP model score estimator   0.05   0.05   0.05   0.05   0.05   0.05   0.05
#> FDP synth data estimator    0.34   0.34   0.34   0.34   0.34   0.34   0.34
#> Linked pairs (augm. RL)   144.00 144.00 144.00 144.00 144.00 144.00 144.00
#> Linked pairs (RL)         148.00 148.00 148.00 148.00 148.00 148.00 148.00
#>                             0.71   0.72   0.73   0.74   0.75   0.76   0.77
#> FDP model score estimator   0.05   0.05   0.05   0.05   0.05   0.05   0.05
#> FDP synth data estimator    0.34   0.34   0.34   0.34   0.34   0.34   0.34
#> Linked pairs (augm. RL)   144.00 144.00 144.00 144.00 144.00 144.00 144.00
#> Linked pairs (RL)         148.00 148.00 148.00 148.00 148.00 148.00 148.00
#>                             0.78   0.79   0.80   0.81   0.82   0.83   0.84
#> FDP model score estimator   0.05   0.05   0.00   0.00   0.00   0.00   0.00
#> FDP synth data estimator    0.34   0.34   0.27   0.27   0.27   0.27   0.27
#> Linked pairs (augm. RL)   144.00 144.00 109.00 109.00 109.00 109.00 109.00
#> Linked pairs (RL)         148.00 148.00 112.00 112.00 112.00 112.00 112.00
#>                             0.85   0.86   0.87   0.88   0.89   0.90   0.91
#> FDP model score estimator   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> FDP synth data estimator    0.27   0.27   0.27   0.27   0.27   0.27   0.27
#> Linked pairs (augm. RL)   109.00 109.00 109.00 109.00 109.00 109.00 109.00
#> Linked pairs (RL)         112.00 112.00 112.00 112.00 112.00 112.00 112.00
#>                             0.92   0.93   0.94   0.95   0.96   0.97   0.98
#> FDP model score estimator   0.00   0.00   0.00   0.00   0.00   0.00   0.00
#> FDP synth data estimator    0.27   0.27   0.27   0.27   0.27   0.27   0.27
#> Linked pairs (augm. RL)   109.00 109.00 109.00 109.00 109.00 109.00 109.00
#> Linked pairs (RL)         112.00 112.00 112.00 112.00 112.00 112.00 112.00
#>                             0.99
#> FDP model score estimator   0.00
#> FDP synth data estimator    0.27
#> Linked pairs (augm. RL)   109.00
#> Linked pairs (RL)         112.00
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1213483 0.1213483 0.1213483 0.1213483 0.1213483 0.1213483 0.1213483
#> 2 0.1016393 0.1016393 0.1016393 0.1016393 0.1016393 0.1016393 0.1016393
#>        0.57      0.58      0.59       0.60       0.61       0.62       0.63
#> 1 0.1213483 0.1213483 0.1213483 0.04822695 0.04822695 0.04822695 0.04822695
#> 2 0.1016393 0.1016393 0.1016393 0.05000000 0.05000000 0.05000000 0.05000000
#>         0.64       0.65       0.66       0.67       0.68       0.69       0.70
#> 1 0.04822695 0.04822695 0.04822695 0.04822695 0.04822695 0.04822695 0.04822695
#> 2 0.05000000 0.05000000 0.05000000 0.05000000 0.05000000 0.05000000 0.05000000
#>         0.71       0.72       0.73       0.74       0.75       0.76       0.77
#> 1 0.04822695 0.04822695 0.04822695 0.04822695 0.04822695 0.04822695 0.04822695
#> 2 0.05000000 0.05000000 0.05000000 0.05000000 0.05000000 0.05000000 0.05000000
#>         0.78       0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> 1 0.04822695 0.04822695    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.05000000 0.05000000    0    0    0    0    0    0    0    0    0    0    0
#>   0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0
#> 
#> $FDP_synth_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.2298851 0.2298851 0.2298851 0.2298851 0.2298851 0.2298851 0.2298851
#> 2 0.3977273 0.3977273 0.3977273 0.3977273 0.3977273 0.3977273 0.3977273
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.2298851 0.2298851 0.2298851 0.2173913 0.2173913 0.2173913 0.2173913
#> 2 0.3977273 0.3977273 0.3977273 0.4697987 0.4697987 0.4697987 0.4697987
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913
#> 2 0.4697987 0.4697987 0.4697987 0.4697987 0.4697987 0.4697987 0.4697987
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913
#> 2 0.4697987 0.4697987 0.4697987 0.4697987 0.4697987 0.4697987 0.4697987
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.2173913 0.2173913 0.1904762 0.1904762 0.1904762 0.1904762 0.1904762
#> 2 0.4697987 0.4697987 0.3539823 0.3539823 0.3539823 0.3539823 0.3539823
#>        0.85      0.86      0.87      0.88      0.89      0.90      0.91
#> 1 0.1904762 0.1904762 0.1904762 0.1904762 0.1904762 0.1904762 0.1904762
#> 2 0.3539823 0.3539823 0.3539823 0.3539823 0.3539823 0.3539823 0.3539823
#>        0.92      0.93      0.94      0.95      0.96      0.97      0.98
#> 1 0.1904762 0.1904762 0.1904762 0.1904762 0.1904762 0.1904762 0.1904762
#> 2 0.3539823 0.3539823 0.3539823 0.3539823 0.3539823 0.3539823 0.3539823
#>        0.99
#> 1 0.1904762
#> 2 0.3539823
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1  174  174  174  174  174  174  174  174  174  174  138  138  138  138  138
#> 2  176  176  176  176  176  176  176  176  176  176  149  149  149  149  149
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1  138  138  138  138  138  138  138  138  138  138  138  138  138  138  138
#> 2  149  149  149  149  149  149  149  149  149  149  149  149  149  149  149
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1  105  105  105  105  105  105  105  105  105  105  105  105  105  105  105
#> 2  113  113  113  113  113  113  113  113  113  113  113  113  113  113  113
#>   0.95 0.96 0.97 0.98 0.99
#> 1  105  105  105  105  105
#> 2  113  113  113  113  113
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1  178  178  178  178  178  178  178  178  178  178  141  141  141  141  141
#> 2  183  183  183  183  183  183  183  183  183  183  156  156  156  156  156
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1  141  141  141  141  141  141  141  141  141  141  141  141  141  141  141
#> 2  156  156  156  156  156  156  156  156  156  156  156  156  156  156  156
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1  107  107  107  107  107  107  107  107  107  107  107  107  107  107  107
#> 2  117  117  117  117  117  117  117  117  117  117  117  117  117  117  117
#>   0.95 0.96 0.97 0.98 0.99
#> 1  107  107  107  107  107
#> 2  117  117  117  117  117
#> 
```
