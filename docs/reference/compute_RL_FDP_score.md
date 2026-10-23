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
#> Warning: PIVs 'V4' and 'V2' are associated (Cramer's V = 0.33 in data1, 0.25 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)                          
compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
                      "BRL", flds = PIVs, types = rep("bi",length(PIVs)) )
#> BRL results (average over 2 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.17    0.16    0.16    0.16    0.16    0.16    0.16    0.14    0   0
#> Linked pairs (RL)            60  58  58  58  58  58  54  22  0   0
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1699058 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678
#> 2 0.1662593 0.1662593 0.1608475 0.1554023 0.1554023 0.1554023 0.1554023
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678
#> 2 0.1554023 0.1554023 0.1554023 0.1554023 0.1554023 0.1554023 0.1554023
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678
#> 2 0.1554023 0.1554023 0.1554023 0.1554023 0.1554023 0.1554023 0.1554023
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678 0.1643678 0.1619246
#> 2 0.1554023 0.1554023 0.1554023 0.1554023 0.1554023 0.1554023 0.1554023
#>        0.78      0.79      0.80      0.81      0.82      0.83      0.84
#> 1 0.1619246 0.1600823 0.1583547 0.1577124 0.1496296 0.1462418 0.1431034
#> 2 0.1554023 0.1532143 0.1532143 0.1524646 0.1511740 0.1476329 0.1405208
#>        0.85      0.86      0.87      0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95
#> 1 0.1372515 0.1295833 0.1247222 0.0000000    0    0    0    0    0    0    0
#> 2 0.1357870 0.1276068 0.1201587 0.1122222    0    0    0    0    0    0    0
#>   0.96 0.97 0.98 0.99
#> 1    0    0    0    0
#> 2    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   59   58   58   58   58   58   58   58   58   58   58   58   58   58   58
#> 2   60   60   59   58   58   58   58   58   58   58   58   58   58   58   58
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   58   58   58   58   58   58   58   58   58   58   58   58   56   56   54
#> 2   58   58   58   58   58   58   58   58   58   58   58   58   58   58   56
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   52   51   39   34   29   19    8    4    0    0    0    0    0    0    0
#> 2   56   55   53   46   32   24   13    7    2    0    0    0    0    0    0
#>   0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0
#> 2    0    0    0    0    0
#> 
compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
                      "FlexRL", data = prep_data, StEM_iter = 5, 
                      StEM_burnin = 2, gibbs_iter = 5, gibbs_burnin = 2,
                      n_post_sample = 5 )
#> FlexRL
#> FlexRL
#> FlexRL results (average over 2 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.13    0.13    0.05    0.05    0.05    0.05    0   0   0   0
#> Linked pairs (RL)            95  95  76  76  76  76  55  55  55  55
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1360825 0.1360825 0.1360825 0.1360825 0.1360825 0.1360825 0.1360825
#> 2 0.1139785 0.1139785 0.1139785 0.1139785 0.1139785 0.1139785 0.1139785
#>        0.57      0.58      0.59       0.60       0.61       0.62       0.63
#> 1 0.1360825 0.1360825 0.1360825 0.06315789 0.06315789 0.06315789 0.06315789
#> 2 0.1139785 0.1139785 0.1139785 0.04533333 0.04533333 0.04533333 0.04533333
#>         0.64       0.65       0.66       0.67       0.68       0.69       0.70
#> 1 0.06315789 0.06315789 0.06315789 0.06315789 0.06315789 0.06315789 0.06315789
#> 2 0.04533333 0.04533333 0.04533333 0.04533333 0.04533333 0.04533333 0.04533333
#>         0.71       0.72       0.73       0.74       0.75       0.76       0.77
#> 1 0.06315789 0.06315789 0.06315789 0.06315789 0.06315789 0.06315789 0.06315789
#> 2 0.04533333 0.04533333 0.04533333 0.04533333 0.04533333 0.04533333 0.04533333
#>         0.78       0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> 1 0.06315789 0.06315789    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.04533333 0.04533333    0    0    0    0    0    0    0    0    0    0    0
#>   0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   97   97   97   97   97   97   97   97   97   97   76   76   76   76   76
#> 2   93   93   93   93   93   93   93   93   93   93   75   75   75   75   75
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   76   76   76   76   76   76   76   76   76   76   76   76   76   76   76
#> 2   75   75   75   75   75   75   75   75   75   75   75   75   75   75   75
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   52   52   52   52   52   52   52   52   52   52   52   52   52   52   52
#> 2   58   58   58   58   58   58   58   58   58   58   58   58   58   58   58
#>   0.95 0.96 0.97 0.98 0.99
#> 1   52   52   52   52   52
#> 2   58   58   58   58   58
#> 
```
