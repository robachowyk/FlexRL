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
(`n_repeats` rows x 50 thresholds, 0.50 to 0.99), and `FDP_score_sd`
Numeric, standard deviation of `FDP_score_estimator` over the
`n_repeats` iterations per threshold (`NA` when fewer than two valid
repeats are available).

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
#> Warning: PIVs 'V4' and 'V2' are associated (Cramer's V = 0.30 in data1, 0.26 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)                          
compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
                      "BRL", flds = PIVs, types = rep("bi",length(PIVs)) )
#> BRL results (average over 2 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold              0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator (mean) 0.17    0.15    0.15    0.15    0.15    0.15    0.15    0.12    0.09    0
#> FDP model score estimator (sd)   0   0   0   0   0   0   0   0   0   0
#> Linked pairs (RL)                61  59  59  59  59  59  51  28  4   0
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1659199 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269
#> 2 0.1659199 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269
#>        0.57      0.58      0.59      0.60      0.61      0.62      0.63
#> 1 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269
#> 2 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269
#>        0.64      0.65      0.66      0.67      0.68      0.69      0.70
#> 1 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269
#> 2 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269
#>        0.71      0.72      0.73      0.74      0.75      0.76      0.77
#> 1 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1533525
#> 2 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1547269 0.1533525
#>        0.78      0.79      0.80     0.81      0.82      0.83      0.84
#> 1 0.1521442 0.1498788 0.1457516 0.141513 0.1352033 0.1300309 0.1231852
#> 2 0.1521442 0.1498788 0.1457516 0.141513 0.1352033 0.1300309 0.1231852
#>        0.85      0.86      0.87      0.88       0.89       0.90       0.91 0.92
#> 1 0.1211508 0.1100585 0.1052083 0.1037778 0.09833333 0.08666667 0.08407407    0
#> 2 0.1211508 0.1100585 0.1052083 0.1037778 0.09833333 0.08666667 0.08407407    0
#>   0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0
#> 
#> $FDP_score_sd
#> 0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64 0.65 
#>    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
#> 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.81 
#>    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
#> 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 
#>    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
#> 0.98 0.99 
#>    0    0 
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   61   59   59   59   59   59   59   59   59   59   59   59   59   59   59
#> 2   61   59   59   59   59   59   59   59   59   59   59   59   59   59   59
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   59   59   59   59   59   59   59   59   59   59   59   59   58   57   55
#> 2   59   59   59   59   59   59   59   59   59   59   59   59   58   57   55
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   51   47   41   36   30   28   19   16   15   10    4    3    0    0    0
#> 2   51   47   41   36   30   28   19   16   15   10    4    3    0    0    0
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
#> Decision rule threshold              0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator (mean) 0.13    0.13    0.05    0.05    0.05    0.05    0   0   0   0
#> FDP model score estimator (sd)   0.01    0.01    0.01    0.01    0.01    0.01    0   0   0   0
#> Linked pairs (RL)                81  81  62  62  62  62  47  47  47  47
#> $FDP_score_estimator
#>        0.50      0.51      0.52      0.53      0.54      0.55      0.56
#> 1 0.1358025 0.1358025 0.1358025 0.1358025 0.1358025 0.1358025 0.1358025
#> 2 0.1234568 0.1234568 0.1234568 0.1234568 0.1234568 0.1234568 0.1234568
#>        0.57      0.58      0.59       0.60       0.61       0.62       0.63
#> 1 0.1358025 0.1358025 0.1358025 0.05483871 0.05483871 0.05483871 0.05483871
#> 2 0.1234568 0.1234568 0.1234568 0.04444444 0.04444444 0.04444444 0.04444444
#>         0.64       0.65       0.66       0.67       0.68       0.69       0.70
#> 1 0.05483871 0.05483871 0.05483871 0.05483871 0.05483871 0.05483871 0.05483871
#> 2 0.04444444 0.04444444 0.04444444 0.04444444 0.04444444 0.04444444 0.04444444
#>         0.71       0.72       0.73       0.74       0.75       0.76       0.77
#> 1 0.05483871 0.05483871 0.05483871 0.05483871 0.05483871 0.05483871 0.05483871
#> 2 0.04444444 0.04444444 0.04444444 0.04444444 0.04444444 0.04444444 0.04444444
#>         0.78       0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90
#> 1 0.05483871 0.05483871    0    0    0    0    0    0    0    0    0    0    0
#> 2 0.04444444 0.04444444    0    0    0    0    0    0    0    0    0    0    0
#>   0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99
#> 1    0    0    0    0    0    0    0    0    0
#> 2    0    0    0    0    0    0    0    0    0
#> 
#> $FDP_score_sd
#>        0.50        0.51        0.52        0.53        0.54        0.55 
#> 0.008729713 0.008729713 0.008729713 0.008729713 0.008729713 0.008729713 
#>        0.56        0.57        0.58        0.59        0.60        0.61 
#> 0.008729713 0.008729713 0.008729713 0.008729713 0.007349855 0.007349855 
#>        0.62        0.63        0.64        0.65        0.66        0.67 
#> 0.007349855 0.007349855 0.007349855 0.007349855 0.007349855 0.007349855 
#>        0.68        0.69        0.70        0.71        0.72        0.73 
#> 0.007349855 0.007349855 0.007349855 0.007349855 0.007349855 0.007349855 
#>        0.74        0.75        0.76        0.77        0.78        0.79 
#> 0.007349855 0.007349855 0.007349855 0.007349855 0.007349855 0.007349855 
#>        0.80        0.81        0.82        0.83        0.84        0.85 
#> 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 
#>        0.86        0.87        0.88        0.89        0.90        0.91 
#> 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 
#>        0.92        0.93        0.94        0.95        0.96        0.97 
#> 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 
#>        0.98        0.99 
#> 0.000000000 0.000000000 
#> 
#> $Linked_pairs
#>   0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64
#> 1   81   81   81   81   81   81   81   81   81   81   62   62   62   62   62
#> 2   81   81   81   81   81   81   81   81   81   81   63   63   63   63   63
#>   0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
#> 1   62   62   62   62   62   62   62   62   62   62   62   62   62   62   62
#> 2   63   63   63   63   63   63   63   63   63   63   63   63   63   63   63
#>   0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94
#> 1   45   45   45   45   45   45   45   45   45   45   45   45   45   45   45
#> 2   49   49   49   49   49   49   49   49   49   49   49   49   49   49   49
#>   0.95 0.96 0.97 0.98 0.99
#> 1   45   45   45   45   45
#> 2   49   49   49   49   49
#> 
```
