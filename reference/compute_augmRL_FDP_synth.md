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
thresholds, 0.50 to 0.99), and `FDP_score_sd`, `FDP_synth_sd` Numeric,
standard deviations of the corresponding estimators over the `n_repeats`
iterations per threshold (`NA` when fewer than two valid repeats are
available).

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
#> Warning: PIVs 'V4' and 'V3' are associated (Cramer's V = 0.30 in data1, 0.25 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)                         
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "BRL", flds = PIVs, 
                          types = rep("bi",length(PIVs)) )
#> Iteration: 0, Accuracy: 48.61%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: 1 real records were dropped from file B during synthesis; the real/synthetic index split used by the FDP_synth estimator is invalid for data that was not support-restricted by prepare_data().
#> Iteration: 0, Accuracy: 49.24%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: 1 real records were dropped from file B during synthesis; the real/synthetic index split used by the FDP_synth estimator is invalid for data that was not support-restricted by prepare_data().
#> BRL results (average over 2 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold              0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator (mean) 0.17    0.16    0.16    0.16    0.16    0.16    0.15    0.13    0.1 0
#> FDP model score estimator (sd)   0   0   0   0   0   0   0   0   0   0
#> FDP synth data estimator (mean)  0.1 0   0   0   0   0   0   0   0   0
#> FDP synth data estimator (sd)    0.14    0   0   0   0   0   0   0   0   0
#> Linked pairs (augm. RL)          56  54  54  54  54  54  44  22  1   0
#> Linked pairs (RL)                56  54  54  54  54  54  44  22  1   0
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1674343 0.1674343 0.1674343 0.1616667 0.1616667 0.1616667 0.1616667
#> 2 0.1724951 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1616667 0.1616667 0.1616667 0.1616667 0.1616667 0.1616667 0.1616667
#> 2 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1616667 0.1616667 0.1616667 0.1616667 0.1616667 0.1616667 0.1616667
#> 2 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.1616667 0.1616667 0.1616667 0.1616667 0.1616667 0.1584829 0.1569935
#> 2 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677 0.1607677 0.1581132
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.1541270 0.1529398 0.1469251 0.1469251 0.1431339 0.1377778 0.1340613
#> 2 0.1555556 0.1519213 0.1495411 0.1473990 0.1438889 0.1385948 0.1306838
#>        0.85      0.86      0.87      0.88      0.89       0.90 0.91 0.92 0.93
#> 1 0.1285507 0.1208889 0.1142222 0.1113889 0.1027778 0.09888889    0    0    0
#> 2 0.1248148 0.1193056 0.1121212 0.1076389 0.1035556 0.09444444    0    0    0
#>   0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0
#> 2    0    0    0    0    0    0
#> 
#> $FDP_synth_estimator
#>        0.50      0.51      0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61
#> 1 0.1949318 0.1949318 0.1949318    0    0    0    0    0    0    0    0    0
#> 2 0.0000000 0.0000000 0.0000000    0    0    0    0    0    0    0    0    0
#>   0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76
#> 1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91
#> 1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0
#> 
#> $FDP_score_sd
#>         0.50         0.51         0.52         0.53         0.54         0.55 
#> 0.0035785142 0.0047140452 0.0047140452 0.0006356819 0.0006356819 0.0006356819 
#>         0.56         0.57         0.58         0.59         0.60         0.61 
#> 0.0006356819 0.0006356819 0.0006356819 0.0006356819 0.0006356819 0.0006356819 
#>         0.62         0.63         0.64         0.65         0.66         0.67 
#> 0.0006356819 0.0006356819 0.0006356819 0.0006356819 0.0006356819 0.0006356819 
#>         0.68         0.69         0.70         0.71         0.72         0.73 
#> 0.0006356819 0.0006356819 0.0006356819 0.0006356819 0.0006356819 0.0006356819 
#>         0.74         0.75         0.76         0.77         0.78         0.79 
#> 0.0006356819 0.0006356819 0.0016155769 0.0007917782 0.0010101525 0.0007202014 
#>         0.80         0.81         0.82         0.83         0.84         0.85 
#> 0.0018497901 0.0003351158 0.0005338555 0.0005777016 0.0023882829 0.0026416872 
#>         0.86         0.87         0.88         0.89         0.90         0.91 
#> 0.0011195857 0.0014856385 0.0026516504 0.0005499719 0.0031426968 0.0000000000 
#>         0.92         0.93         0.94         0.95         0.96         0.97 
#> 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
#>         0.98         0.99 
#> 0.0000000000 0.0000000000 
#> 
#> $FDP_synth_sd
#>      0.50      0.51      0.52      0.53      0.54      0.55      0.56      0.57 
#> 0.1378376 0.1378376 0.1378376 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
#>      0.58      0.59      0.60      0.61      0.62      0.63      0.64      0.65 
#> 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
#>      0.66      0.67      0.68      0.69      0.70      0.71      0.72      0.73 
#> 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
#>      0.74      0.75      0.76      0.77      0.78      0.79      0.80      0.81 
#> 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
#>      0.82      0.83      0.84      0.85      0.86      0.87      0.88      0.89 
#> 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
#>      0.90      0.91      0.92      0.93      0.94      0.95      0.96      0.97 
#> 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
#>      0.98      0.99 
#> 0.0000000 0.0000000 
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   54   54   54   54   54   54   54   54   54   54   54   54   54   54   54
#> 2   57   55   55   55   55   55   55   55   55   55   55   55   55   55   55
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   54   54   54   54   54   54   54   54   54   54   54   52   51   49   48
#> 2   55   55   55   55   55   55   55   55   55   55   55   55   53   51   48
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   43   43   39   33   29   23   15   10    8    2    1    0    0    0    0
#> 2   46   44   40   34   26   21   16   11    8    5    1    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   55   55   55   54   54   54   54   54   54   54   54   54   54   54   54
#> 2   57   55   55   55   55   55   55   55   55   55   55   55   55   55   55
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   54   54   54   54   54   54   54   54   54   54   54   52   51   49   48
#> 2   55   55   55   55   55   55   55   55   55   55   55   55   53   51   48
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   43   43   39   33   29   23   15   10    8    2    1    0    0    0    0
#> 2   46   44   40   34   26   21   16   11    8    5    1    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "FlexRL", data = prep_data, 
                          StEM_iter = 5, StEM_burnin = 2, 
                          gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5 )
#> Iteration: 0, Accuracy: 49.24%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: 1 real records were dropped from file B during synthesis; the real/synthetic index split used by the FDP_synth estimator is invalid for data that was not support-restricted by prepare_data().
#> FlexRL
#> Iteration: 0, Accuracy: 48.1%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> Warning: 1 real records were dropped from file B during synthesis; the real/synthetic index split used by the FDP_synth estimator is invalid for data that was not support-restricted by prepare_data().
#> FlexRL
#> FlexRL results (average over 2 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold              0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator (mean) 0.11    0.11    0.05    0.05    0.05    0.05    0   0   0   0
#> FDP model score estimator (sd)   0.02    0.02    0.02    0.02    0.02    0.02    0   0   0   0
#> FDP synth data estimator (mean)  0.52    0.52    0.63    0.63    0.63    0.63    0.88    0.88    0.88    0.88
#> FDP synth data estimator (sd)    0.08    0.08    0.12    0.12    0.12    0.12    0.28    0.28    0.28    0.28
#> Linked pairs (augm. RL)          88  88  72  72  72  72  54  54  54  54
#> Linked pairs (RL)                92  92  76  76  76  76  58  58  58  58
#> $FDP_score_estimator
#>         0.50       0.51       0.52       0.53       0.54       0.55       0.56
#> 1 0.12467532 0.12467532 0.12467532 0.12467532 0.12467532 0.12467532 0.12467532
#> 2 0.09245283 0.09245283 0.09245283 0.09245283 0.09245283 0.09245283 0.09245283
#>         0.57       0.58       0.59       0.60       0.61       0.62       0.63
#> 1 0.12467532 0.12467532 0.12467532 0.06349206 0.06349206 0.06349206 0.06349206
#> 2 0.09245283 0.09245283 0.09245283 0.03777778 0.03777778 0.03777778 0.03777778
#>         0.64       0.65       0.66       0.67       0.68       0.69       0.70
#> 1 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206
#> 2 0.03777778 0.03777778 0.03777778 0.03777778 0.03777778 0.03777778 0.03777778
#>         0.71       0.72       0.73       0.74       0.75       0.76       0.77
#> 1 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206 0.06349206
#> 2 0.03777778 0.03777778 0.03777778 0.03777778 0.03777778 0.03777778 0.03777778
#>         0.78       0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> 1 0.06349206 0.06349206    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.03777778 0.03777778    0    0    0    0    0    0    0    0    0    0    0
#>   0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0
#> 
#> $FDP_synth_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.5767844 0.5767844 0.5767844 0.5767844 0.5767844 0.5767844 0.5767844
#> 2 0.4613610 0.4613610 0.4613610 0.4613610 0.4613610 0.4613610 0.4613610
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.5767844 0.5767844 0.5767844 0.7136485 0.7136485 0.7136485 0.7136485
#> 2 0.4613610 0.4613610 0.4613610 0.5471956 0.5471956 0.5471956 0.5471956
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.7136485 0.7136485 0.7136485 0.7136485 0.7136485 0.7136485 0.7136485
#> 2 0.5471956 0.5471956 0.5471956 0.5471956 0.5471956 0.5471956 0.5471956
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.7136485 0.7136485 0.7136485 0.7136485 0.7136485 0.7136485 0.7136485
#> 2 0.5471956 0.5471956 0.5471956 0.5471956 0.5471956 0.5471956 0.5471956
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.7136485 0.7136485 1.0796221 1.0796221 1.0796221 1.0796221 1.0796221
#> 2 0.5471956 0.5471956 0.6820119 0.6820119 0.6820119 0.6820119 0.6820119
#>        0.85      0.86      0.87      0.88      0.89      0.90      0.91
#> 1 1.0796221 1.0796221 1.0796221 1.0796221 1.0796221 1.0796221 1.0796221
#> 2 0.6820119 0.6820119 0.6820119 0.6820119 0.6820119 0.6820119 0.6820119
#>        0.92      0.93      0.94      0.95      0.96      0.97      0.98
#> 1 1.0796221 1.0796221 1.0796221 1.0796221 1.0796221 1.0796221 1.0796221
#> 2 0.6820119 0.6820119 0.6820119 0.6820119 0.6820119 0.6820119 0.6820119
#>        0.99
#> 1 1.0796221
#> 2 0.6820119
#> 
#> $FDP_score_sd
#>       0.50       0.51       0.52       0.53       0.54       0.55       0.56 
#> 0.02278474 0.02278474 0.02278474 0.02278474 0.02278474 0.02278474 0.02278474 
#>       0.57       0.58       0.59       0.60       0.61       0.62       0.63 
#> 0.02278474 0.02278474 0.02278474 0.01818275 0.01818275 0.01818275 0.01818275 
#>       0.64       0.65       0.66       0.67       0.68       0.69       0.70 
#> 0.01818275 0.01818275 0.01818275 0.01818275 0.01818275 0.01818275 0.01818275 
#>       0.71       0.72       0.73       0.74       0.75       0.76       0.77 
#> 0.01818275 0.01818275 0.01818275 0.01818275 0.01818275 0.01818275 0.01818275 
#>       0.78       0.79       0.80       0.81       0.82       0.83       0.84 
#> 0.01818275 0.01818275 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 
#>       0.85       0.86       0.87       0.88       0.89       0.90       0.91 
#> 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 
#>       0.92       0.93       0.94       0.95       0.96       0.97       0.98 
#> 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 
#>       0.99 
#> 0.00000000 
#> 
#> $FDP_synth_sd
#>       0.50       0.51       0.52       0.53       0.54       0.55       0.56 
#> 0.08161668 0.08161668 0.08161668 0.08161668 0.08161668 0.08161668 0.08161668 
#>       0.57       0.58       0.59       0.60       0.61       0.62       0.63 
#> 0.08161668 0.08161668 0.08161668 0.11769998 0.11769998 0.11769998 0.11769998 
#>       0.64       0.65       0.66       0.67       0.68       0.69       0.70 
#> 0.11769998 0.11769998 0.11769998 0.11769998 0.11769998 0.11769998 0.11769998 
#>       0.71       0.72       0.73       0.74       0.75       0.76       0.77 
#> 0.11769998 0.11769998 0.11769998 0.11769998 0.11769998 0.11769998 0.11769998 
#>       0.78       0.79       0.80       0.81       0.82       0.83       0.84 
#> 0.11769998 0.11769998 0.28115287 0.28115287 0.28115287 0.28115287 0.28115287 
#>       0.85       0.86       0.87       0.88       0.89       0.90       0.91 
#> 0.28115287 0.28115287 0.28115287 0.28115287 0.28115287 0.28115287 0.28115287 
#>       0.92       0.93       0.94       0.95       0.96       0.97       0.98 
#> 0.28115287 0.28115287 0.28115287 0.28115287 0.28115287 0.28115287 0.28115287 
#>       0.99 
#> 0.28115287 
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   73   73   73   73   73   73   73   73   73   73   59   59   59   59   59
#> 2  102  102  102  102  102  102  102  102  102  102   86   86   86   86   86
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   59   59   59   59   59   59   59   59   59   59   59   59   59   59   59
#> 2   86   86   86   86   86   86   86   86   86   86   86   86   86   86   86
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   39   39   39   39   39   39   39   39   39   39   39   39   39   39   39
#> 2   69   69   69   69   69   69   69   69   69   69   69   69   69   69   69
#>   0.95 0.96 0.97 0.98 0.99
#> 1   39   39   39   39   39
#> 2   69   69   69   69   69
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   77   77   77   77   77   77   77   77   77   77   63   63   63   63   63
#> 2  106  106  106  106  106  106  106  106  106  106   90   90   90   90   90
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   63   63   63   63   63   63   63   63   63   63   63   63   63   63   63
#> 2   90   90   90   90   90   90   90   90   90   90   90   90   90   90   90
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   43   43   43   43   43   43   43   43   43   43   43   43   43   43   43
#> 2   73   73   73   73   73   73   73   73   73   73   73   73   73   73   73
#>   0.95 0.96 0.97 0.98 0.99
#> 1   43   43   43   43   43
#> 2   73   73   73   73   73
#> 
```
