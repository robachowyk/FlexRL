# Estimate the false discovery proportion of a record-linkage method via synthetic augmentation

Repeatedly augments file B with synthetic records
([`synthesise()`](https://robachowyk.github.io/FlexRL/reference/synthesise.md)),
runs the chosen record-linkage method link_with\_\*, and compares the
synthetic-vs-real proportion among linked pairs to estimate the false
discovery proportion, for a range of score thresholds (0.50 to 0.99).
See https://doi.org/10.1002/sim.70292 for the method.

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

  Passed to
  [`synthesise()`](https://robachowyk.github.io/FlexRL/reference/synthesise.md):
  `"arf"`, `"synthpop"`, or `"mice"`.

- encodedA, encodedB:

  The two encoded data sources as returned by
  [`prepare_data()`](https://robachowyk.github.io/FlexRL/reference/prepare_data.md)
  (`encodedA` is the smaller one).

- PIVs:

  Character vector, names of the PIVs.

- n_synth:

  Integer, number of synthetic records to generate per iteration
  (default: 10% of `nrow(encodedB)`).

- restrict_support_intersection:

  Passed to
  [`synthesise()`](https://robachowyk.github.io/FlexRL/reference/synthesise.md).

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
gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
                           p_mistake, p_missing, cond_hazard_params, TRUE )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> Warning: PIVs 'V4' and 'V1' are associated (Cramer's V = 0.29 in data1, 0.33 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)                         
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "BRL", flds = PIVs, 
                          types = rep("bi",length(PIVs)) )
#> Iteration: 0, Accuracy: 45.94%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> Iteration: 0, Accuracy: 46.85%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> BRL results (average over 2 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.29    0.29    0.29    0.28    0.27    0.24    0   0   0   0
#> FDP synth data estimator     0.11    0.11    0.11    0   0   0   0   0   0   0
#> Linked pairs (augm. RL)      46  46  46  45  34  6   0   0   0   0
#> Linked pairs (RL)            47  47  47  45  34  6   0   0   0   0
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.2871739 0.2871739 0.2871739 0.2871739 0.2871739 0.2871739 0.2871739
#> 2 0.2830556 0.2830556 0.2830556 0.2830556 0.2830556 0.2830556 0.2830556
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.2871739 0.2871739 0.2871739 0.2871739 0.2871739 0.2871739 0.2852593
#> 2 0.2830556 0.2830556 0.2830556 0.2830556 0.2830556 0.2830556 0.2830556
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.2835606 0.2835606 0.2792683 0.2792683 0.2767521 0.2734259 0.2708081
#> 2 0.2813475 0.2797101 0.2781975 0.2757623 0.2745503 0.2693994 0.2665033
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77 0.78
#> 1 0.2641453 0.2598413 0.2547222 0.2438889 0.2371111 0.2340741 0.2255556    0
#> 2 0.2636559 0.2613889 0.2559444 0.2396825 0.2374074 0.2300000 0.2233333    0
#>   0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93
#> 1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0
#> 2    0    0    0    0    0    0
#> 
#> $FDP_synth_estimator
#>       0.50     0.51     0.52     0.53     0.54     0.55     0.56     0.57
#> 1 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
#> 2 0.212766 0.212766 0.212766 0.212766 0.212766 0.212766 0.212766 0.212766
#>       0.58     0.59     0.60     0.61     0.62     0.63 0.64 0.65 0.66 0.67
#> 1 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000    0    0    0    0
#> 2 0.212766 0.212766 0.212766 0.212766 0.212766 0.212766    0    0    0    0
#>   0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82
#> 1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97
#> 1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.98 0.99
#> 1    0    0
#> 2    0    0
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   46   46   46   46   46   46   46   46   46   46   46   46   46   45   44
#> 2   47   47   47   47   47   47   47   47   47   47   47   47   47   47   47
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   44   41   41   39   36   33   26   21   16    8    5    3    1    0    0
#> 2   46   45   43   42   37   34   31   28   20    7    6    3    1    0    0
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   46   46   46   46   46   46   46   46   46   46   46   46   46   45   44
#> 2   48   48   48   48   48   48   48   48   48   48   48   48   48   48   47
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   44   41   41   39   36   33   26   21   16    8    5    3    1    0    0
#> 2   46   45   43   42   37   34   31   28   20    7    6    3    1    0    0
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "FlexRL", data = prep_data, 
                          StEM_iter = 5, StEM_burnin = 2, 
                          gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5 )
#> Iteration: 0, Accuracy: 50.75%
#> Iteration: 1, Accuracy: 37.84%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> Iteration: 0, Accuracy: 52.02%
#> Iteration: 1, Accuracy: 38.58%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> FlexRL results (average over 2 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.16    0.16    0.06    0.06    0.06    0.06    0   0   0   0
#> FDP synth data estimator     0.51    0.51    0.4 0.4 0.4 0.4 0.23    0.23    0.23    0.23
#> Linked pairs (augm. RL)      88  88  64  64  64  64  44  44  44  44
#> Linked pairs (RL)            92  92  66  66  66  66  45  45  45  45
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1695652 0.1695652 0.1695652 0.1695652 0.1695652 0.1695652 0.1695652
#> 2 0.1483871 0.1483871 0.1483871 0.1483871 0.1483871 0.1483871 0.1483871
#>        0.57      0.58      0.59       0.60       0.61       0.62       0.63
#> 1 0.1695652 0.1695652 0.1695652 0.06349206 0.06349206 0.06349206 0.06349206
#> 2 0.1483871 0.1483871 0.1483871 0.06571429 0.06571429 0.06571429 0.06571429
#>         0.64       0.65       0.66       0.67       0.68       0.69       0.70
#> 1 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206
#> 2 0.06571429 0.06571429 0.06571429 0.06571429 0.06571429 0.06571429 0.06571429
#>         0.71       0.72       0.73       0.74       0.75       0.76       0.77
#> 1 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206
#> 2 0.06571429 0.06571429 0.06571429 0.06571429 0.06571429 0.06571429 0.06571429
#>         0.78       0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> 1 0.06349206 0.06349206    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.06571429 0.06571429    0    0    0    0    0    0    0    0    0    0    0
#>   0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0
#> 
#> $FDP_synth_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.5747126 0.5747126 0.5747126 0.5747126 0.5747126 0.5747126 0.5747126
#> 2 0.4494382 0.4494382 0.4494382 0.4494382 0.4494382 0.4494382 0.4494382
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.5747126 0.5747126 0.5747126 0.5000000 0.5000000 0.5000000 0.5000000
#> 2 0.4494382 0.4494382 0.4494382 0.2941176 0.2941176 0.2941176 0.2941176
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
#> 2 0.2941176 0.2941176 0.2941176 0.2941176 0.2941176 0.2941176 0.2941176
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
#> 2 0.2941176 0.2941176 0.2941176 0.2941176 0.2941176 0.2941176 0.2941176
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.5000000 0.5000000 0.2380952 0.2380952 0.2380952 0.2380952 0.2380952
#> 2 0.2941176 0.2941176 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913
#>        0.85      0.86      0.87      0.88      0.89      0.90      0.91
#> 1 0.2380952 0.2380952 0.2380952 0.2380952 0.2380952 0.2380952 0.2380952
#> 2 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913
#>        0.92      0.93      0.94      0.95      0.96      0.97      0.98
#> 1 0.2380952 0.2380952 0.2380952 0.2380952 0.2380952 0.2380952 0.2380952
#> 2 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913 0.2173913
#>        0.99
#> 1 0.2380952
#> 2 0.2173913
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   87   87   87   87   87   87   87   87   87   87   60   60   60   60   60
#> 2   89   89   89   89   89   89   89   89   89   89   68   68   68   68   68
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   60   60   60   60   60   60   60   60   60   60   60   60   60   60   60
#> 2   68   68   68   68   68   68   68   68   68   68   68   68   68   68   68
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   42   42   42   42   42   42   42   42   42   42   42   42   42   42   42
#> 2   46   46   46   46   46   46   46   46   46   46   46   46   46   46   46
#>   0.95 0.96 0.97 0.98 0.99
#> 1   42   42   42   42   42
#> 2   46   46   46   46   46
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   92   92   92   92   92   92   92   92   92   92   63   63   63   63   63
#> 2   93   93   93   93   93   93   93   93   93   93   70   70   70   70   70
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   63   63   63   63   63   63   63   63   63   63   63   63   63   63   63
#> 2   70   70   70   70   70   70   70   70   70   70   70   70   70   70   70
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   43   43   43   43   43   43   43   43   43   43   43   43   43   43   43
#> 2   47   47   47   47   47   47   47   47   47   47   47   47   47   47   47
#>   0.95 0.96 0.97 0.98 0.99
#> 1   43   43   43   43   43
#> 2   47   47   47   47   47
#> 
```
