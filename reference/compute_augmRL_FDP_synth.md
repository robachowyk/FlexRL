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
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)                         
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "BRL", flds = PIVs, 
                          types = rep("bi",length(PIVs)) )
#> Iteration: 0, Accuracy: 51.27%
#> Iteration: 1, Accuracy: 36.2%
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> Iteration: 0, Accuracy: 45.59%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> BRL results (average over 2 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.22    0.22    0.22    0.22    0.21    0.2 0.18    0.07    0   0
#> FDP synth data estimator     0.1 0.1 0.1 0.1 0.1 0   0   0   0   0
#> Linked pairs (augm. RL)      54  53  53  53  52  42  20  1   0   0
#> Linked pairs (RL)            54  54  54  54  52  42  20  1   0   0
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.2207757 0.2207757 0.2207757 0.2207757 0.2207757 0.2207757 0.2207757
#> 2 0.2180606 0.2180606 0.2180606 0.2180606 0.2180606 0.2137037 0.2137037
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.2207757 0.2207757 0.2207757 0.2207757 0.2207757 0.2207757 0.2207757
#> 2 0.2137037 0.2137037 0.2137037 0.2137037 0.2137037 0.2137037 0.2137037
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.2207757 0.2207757 0.2207757 0.2207757 0.2207757 0.2171024 0.2171024
#> 2 0.2137037 0.2137037 0.2137037 0.2137037 0.2137037 0.2118029 0.2118029
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.2137642 0.2106619 0.2079506 0.2066414 0.2032249 0.2000877 0.1951178
#> 2 0.2118029 0.2089760 0.2062358 0.2049769 0.1978571 0.1956111 0.1925526
#>       0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.190728 0.1889712 0.1825731 0.1754630 0.1713889 0.1616667 0.0000000
#> 2 0.186810 0.1816239 0.1762963 0.1729012 0.1636667 0.1564815 0.1477778
#>        0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98
#> 1 0.0000000    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.1444444    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.99
#> 1    0
#> 2    0
#> 
#> $FDP_synth_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077
#> 2 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077
#> 2 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#>        0.64      0.65      0.66      0.67      0.68 0.69 0.70      0.71
#> 1 0.1923077 0.1923077 0.1923077 0.1923077 0.1923077  0.2  0.2 0.2083333
#> 2 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000  0.0  0.0 0.0000000
#>        0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.84 0.85
#> 1 0.2173913    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.0000000    0    0    0    0    0    0    0    0    0    0    0    0    0
#>   0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   52   52   52   52   52   52   52   52   52   52   52   52   52   52   52
#> 2   55   55   55   55   55   54   54   54   54   54   54   54   54   54   54
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   52   52   52   52   50   50   48   46   45   44   41   38   33   29   27
#> 2   54   54   54   54   53   53   53   51   49   48   42   40   37   31   26
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   19   12    8    2    0    0    0    0    0    0    0    0    0    0    0
#> 2   21   18   10    6    3    2    0    0    0    0    0    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   53   53   53   53   53   53   53   53   53   53   53   53   53   53   53
#> 2   55   55   55   55   55   54   54   54   54   54   54   54   54   54   54
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   53   53   53   53   51   51   49   47   45   44   41   38   33   29   27
#> 2   54   54   54   54   53   53   53   51   49   48   42   40   37   31   26
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   19   12    8    2    0    0    0    0    0    0    0    0    0    0    0
#> 2   21   18   10    6    3    2    0    0    0    0    0    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "FlexRL", data = prep_data, 
                          StEM_iter = 5, StEM_burnin = 2, 
                          gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5 )
#> Iteration: 0, Accuracy: 45.2%
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> Iteration: 0, Accuracy: 48.61%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> FlexRL results (average over 2 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.14    0.14    0.07    0.07    0.07    0.07    0   0   0   0
#> FDP synth data estimator     0.26    0.26    0.27    0.27    0.27    0.27    0.2 0.2 0.2 0.2
#> Linked pairs (augm. RL)      96  96  76  76  76  76  50  50  50  50
#> Linked pairs (RL)            99  99  78  78  78  78  52  52  52  52
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1387755 0.1387755 0.1387755 0.1387755 0.1387755 0.1387755 0.1387755
#> 2 0.1380000 0.1380000 0.1380000 0.1380000 0.1380000 0.1380000 0.1380000
#>        0.57      0.58      0.59       0.60       0.61       0.62       0.63
#> 1 0.1387755 0.1387755 0.1387755 0.06753247 0.06753247 0.06753247 0.06753247
#> 2 0.1380000 0.1380000 0.1380000 0.06835443 0.06835443 0.06835443 0.06835443
#>         0.64       0.65       0.66       0.67       0.68       0.69       0.70
#> 1 0.06753247 0.06753247 0.06753247 0.06753247 0.06753247 0.06753247 0.06753247
#> 2 0.06835443 0.06835443 0.06835443 0.06835443 0.06835443 0.06835443 0.06835443
#>         0.71       0.72       0.73       0.74       0.75       0.76       0.77
#> 1 0.06753247 0.06753247 0.06753247 0.06753247 0.06753247 0.06753247 0.06753247
#> 2 0.06835443 0.06835443 0.06835443 0.06835443 0.06835443 0.06835443 0.06835443
#>         0.78       0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> 1 0.06753247 0.06753247    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.06835443 0.06835443    0    0    0    0    0    0    0    0    0    0    0
#>   0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0
#> 
#> $FDP_synth_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.3157895 0.3157895 0.3157895 0.3157895 0.3157895 0.3157895 0.3157895
#> 2 0.2040816 0.2040816 0.2040816 0.2040816 0.2040816 0.2040816 0.2040816
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.3157895 0.3157895 0.3157895 0.4054054 0.4054054 0.4054054 0.4054054
#> 2 0.2040816 0.2040816 0.2040816 0.1282051 0.1282051 0.1282051 0.1282051
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.4054054 0.4054054 0.4054054 0.4054054 0.4054054 0.4054054 0.4054054
#> 2 0.1282051 0.1282051 0.1282051 0.1282051 0.1282051 0.1282051 0.1282051
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.4054054 0.4054054 0.4054054 0.4054054 0.4054054 0.4054054 0.4054054
#> 2 0.1282051 0.1282051 0.1282051 0.1282051 0.1282051 0.1282051 0.1282051
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.4054054 0.4054054 0.4081633 0.4081633 0.4081633 0.4081633 0.4081633
#> 2 0.1282051 0.1282051 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#>        0.85      0.86      0.87      0.88      0.89      0.90      0.91
#> 1 0.4081633 0.4081633 0.4081633 0.4081633 0.4081633 0.4081633 0.4081633
#> 2 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#>        0.92      0.93      0.94      0.95      0.96      0.97      0.98
#> 1 0.4081633 0.4081633 0.4081633 0.4081633 0.4081633 0.4081633 0.4081633
#> 2 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#>        0.99
#> 1 0.4081633
#> 2 0.0000000
#> 
#> $Linked_pairs_augm
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   95   95   95   95   95   95   95   95   95   95   74   74   74   74   74
#> 2   98   98   98   98   98   98   98   98   98   98   78   78   78   78   78
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   74   74   74   74   74   74   74   74   74   74   74   74   74   74   74
#> 2   78   78   78   78   78   78   78   78   78   78   78   78   78   78   78
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   49   49   49   49   49   49   49   49   49   49   49   49   49   49   49
#> 2   52   52   52   52   52   52   52   52   52   52   52   52   52   52   52
#>   0.95 0.96 0.97 0.98 0.99
#> 1   49   49   49   49   49
#> 2   52   52   52   52   52
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   98   98   98   98   98   98   98   98   98   98   77   77   77   77   77
#> 2  100  100  100  100  100  100  100  100  100  100   79   79   79   79   79
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   77   77   77   77   77   77   77   77   77   77   77   77   77   77   77
#> 2   79   79   79   79   79   79   79   79   79   79   79   79   79   79   79
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   51   51   51   51   51   51   51   51   51   51   51   51   51   51   51
#> 2   52   52   52   52   52   52   52   52   52   52   52   52   52   52   52
#>   0.95 0.96 0.97 0.98 0.99
#> 1   51   51   51   51   51
#> 2   52   52   52   52   52
#> 
```
