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
#> Warning: PIVs 'V4' and 'V1' are associated (Cramer's V = 0.32 in data1, 0.22 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)                         
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "BRL", flds = PIVs, 
                          types = rep("bi",length(PIVs)) )
#> Iteration: 0, Accuracy: 50.25%
#> Iteration: 1, Accuracy: 43.04%
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Iteration: 0, Accuracy: 48.24%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> BRL results (average over 2 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.18    0.18    0.18    0.18    0.18    0.18    0.16    0.14    0   0
#> FDP synth data estimator     0.21    0.21    0.21    0.21    0.21    0.21    0.29    0.44    0   0
#> Linked pairs (augm. RL)      60  60  60  60  60  58  43  12  0   0
#> Linked pairs (RL)            60  60  60  60  60  59  44  13  0   0
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1813936 0.1813936 0.1813936 0.1813936 0.1813936 0.1813936 0.1813936
#> 2 0.1733513 0.1733513 0.1733513 0.1733513 0.1733513 0.1733513 0.1733513
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1813936 0.1813936 0.1813936 0.1813936 0.1813936 0.1813936 0.1813936
#> 2 0.1733513 0.1733513 0.1733513 0.1733513 0.1733513 0.1733513 0.1733513
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1813936 0.1813936 0.1813936 0.1813936 0.1813936 0.1813936 0.1813936
#> 2 0.1733513 0.1733513 0.1733513 0.1733513 0.1733513 0.1733513 0.1733513
#>        0.71      0.72      0.73      0.74      0.75      0.76     0.77
#> 1 0.1813936 0.1813936 0.1813936 0.1800192 0.1787135 0.1752675 0.169229
#> 2 0.1733513 0.1733513 0.1733513 0.1733513 0.1720219 0.1720219 0.171000
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.1657971 0.1622222 0.1589167 0.1569006 0.1543492 0.1500383 0.1432778
#> 2 0.1692337 0.1648008 0.1606944 0.1584198 0.1528228 0.1494097 0.1441667
#>        0.85      0.86      0.87      0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95
#> 1 0.1367521 0.1277778 0.1236111 0.1166667    0    0    0    0    0    0    0
#> 2 0.1345299 0.1286111 0.1235556 0.1144444    0    0    0    0    0    0    0
#>   0.96 0.97 0.98 0.99
#> 1    0    0    0    0
#> 2    0    0    0    0
#> 
#> $FDP_synth_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1814882 0.1814882 0.1814882 0.1814882 0.1814882 0.1814882 0.1814882
#> 2 0.2341920 0.2341920 0.2341920 0.2341920 0.2341920 0.2341920 0.2341920
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1814882 0.1814882 0.1814882 0.1814882 0.1814882 0.1814882 0.1814882
#> 2 0.2341920 0.2341920 0.2341920 0.2341920 0.2341920 0.2341920 0.2341920
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1814882 0.1814882 0.1814882 0.1814882 0.1814882 0.1814882 0.1814882
#> 2 0.2341920 0.2341920 0.2341920 0.2341920 0.2341920 0.2341920 0.2341920
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.1814882 0.1814882 0.1814882 0.1846722 0.1879699 0.1986097 0.2192982
#> 2 0.2341920 0.2341920 0.2341920 0.2341920 0.2380952 0.2380952 0.2421308
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.2339181 0.2506266 0.2699055 0.2844950 0.3095975 0.3759398 0.5540166
#> 2 0.2506266 0.2747253 0.3039514 0.3246753 0.0000000 0.0000000 0.0000000
#>       0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98
#> 1 0.877193    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.000000    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.99
#> 1    0
#> 2    0
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   58   58   58   58   58   58   58   58   58   58   58   58   58   58   58
#> 2   61   61   61   61   61   61   61   61   61   61   61   61   61   61   61
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   58   58   58   58   58   58   58   58   58   57   56   53   48   45   42
#> 2   61   61   61   61   61   61   61   61   61   61   60   60   59   57   52
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   39   37   34   28   19   12    7    4    1    0    0    0    0    0    0
#> 2   47   44   37   32   24   13    8    5    1    0    0    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   59   59   59   59   59   59   59   59   59   59   59   59   59   59   59
#> 2   62   62   62   62   62   62   62   62   62   62   62   62   62   62   62
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   59   59   59   59   59   59   59   59   59   58   57   54   49   46   43
#> 2   62   62   62   62   62   62   62   62   62   62   61   61   60   58   53
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   40   38   35   29   20   13    7    4    1    0    0    0    0    0    0
#> 2   48   45   37   32   24   13    8    5    1    0    0    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "FlexRL", data = prep_data, 
                          StEM_iter = 5, StEM_burnin = 2, 
                          gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5 )
#> Iteration: 0, Accuracy: 49.62%
#> FlexRL
#> Iteration: 0, Accuracy: 47.58%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> FlexRL
#> FlexRL results (average over 2 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.15    0.15    0.06    0.06    0.06    0.06    0   0   0   0
#> FDP synth data estimator     0.82    0.82    0.55    0.55    0.55    0.55    0.46    0.46    0.46    0.46
#> Linked pairs (augm. RL)      80  80  60  60  60  60  42  42  42  42
#> Linked pairs (RL)            86  86  64  64  64  64  44  44  44  44
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1769231 0.1769231 0.1769231 0.1769231 0.1769231 0.1769231 0.1769231
#> 2 0.1263158 0.1263158 0.1263158 0.1263158 0.1263158 0.1263158 0.1263158
#>        0.57      0.58      0.59       0.60       0.61       0.62       0.63
#> 1 0.1769231 0.1769231 0.1769231 0.07169811 0.07169811 0.07169811 0.07169811
#> 2 0.1263158 0.1263158 0.1263158 0.05333333 0.05333333 0.05333333 0.05333333
#>         0.64       0.65       0.66       0.67       0.68       0.69       0.70
#> 1 0.07169811 0.07169811 0.07169811 0.07169811 0.07169811 0.07169811 0.07169811
#> 2 0.05333333 0.05333333 0.05333333 0.05333333 0.05333333 0.05333333 0.05333333
#>         0.71       0.72       0.73       0.74       0.75       0.76       0.77
#> 1 0.07169811 0.07169811 0.07169811 0.07169811 0.07169811 0.07169811 0.07169811
#> 2 0.05333333 0.05333333 0.05333333 0.05333333 0.05333333 0.05333333 0.05333333
#>         0.78       0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> 1 0.07169811 0.07169811    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.05333333 0.05333333    0    0    0    0    0    0    0    0    0    0    0
#>   0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0
#> 
#> $FDP_synth_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.5405405 0.5405405 0.5405405 0.5405405 0.5405405 0.5405405 0.5405405
#> 2 1.1015912 1.1015912 1.1015912 1.1015912 1.1015912 1.1015912 1.1015912
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.5405405 0.5405405 0.5405405 0.1923077 0.1923077 0.1923077 0.1923077
#> 2 1.1015912 1.1015912 1.1015912 0.9153318 0.9153318 0.9153318 0.9153318
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077
#> 2 0.9153318 0.9153318 0.9153318 0.9153318 0.9153318 0.9153318 0.9153318
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077
#> 2 0.9153318 0.9153318 0.9153318 0.9153318 0.9153318 0.9153318 0.9153318
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.1923077 0.1923077 0.3030303 0.3030303 0.3030303 0.3030303 0.3030303
#> 2 0.9153318 0.9153318 0.6072874 0.6072874 0.6072874 0.6072874 0.6072874
#>        0.85      0.86      0.87      0.88      0.89      0.90      0.91
#> 1 0.3030303 0.3030303 0.3030303 0.3030303 0.3030303 0.3030303 0.3030303
#> 2 0.6072874 0.6072874 0.6072874 0.6072874 0.6072874 0.6072874 0.6072874
#>        0.92      0.93      0.94      0.95      0.96      0.97      0.98
#> 1 0.3030303 0.3030303 0.3030303 0.3030303 0.3030303 0.3030303 0.3030303
#> 2 0.6072874 0.6072874 0.6072874 0.6072874 0.6072874 0.6072874 0.6072874
#>        0.99
#> 1 0.3030303
#> 2 0.6072874
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   74   74   74   74   74   74   74   74   74   74   52   52   52   52   52
#> 2   86   86   86   86   86   86   86   86   86   86   69   69   69   69   69
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   52   52   52   52   52   52   52   52   52   52   52   52   52   52   52
#> 2   69   69   69   69   69   69   69   69   69   69   69   69   69   69   69
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   33   33   33   33   33   33   33   33   33   33   33   33   33   33   33
#> 2   52   52   52   52   52   52   52   52   52   52   52   52   52   52   52
#>   0.95 0.96 0.97 0.98 0.99
#> 1   33   33   33   33   33
#> 2   52   52   52   52   52
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   78   78   78   78   78   78   78   78   78   78   53   53   53   53   53
#> 2   95   95   95   95   95   95   95   95   95   95   75   75   75   75   75
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   53   53   53   53   53   53   53   53   53   53   53   53   53   53   53
#> 2   75   75   75   75   75   75   75   75   75   75   75   75   75   75   75
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   34   34   34   34   34   34   34   34   34   34   34   34   34   34   34
#> 2   55   55   55   55   55   55   55   55   55   55   55   55   55   55   55
#>   0.95 0.96 0.97 0.98 0.99
#> 1   34   34   34   34   34
#> 2   55   55   55   55   55
#> 
```
