# Post-linkage diagnostics

Gathers, from a fitted linkage, the diagnostics needed to judge the
linked data before using it for inference: false discovery proportion
estimate, agreement among PIVs of linked pairs, comparison of the linked
subset with each source file (standardised mean differences, support
overlap, maximum mean discrepancy).

## Usage

``` r
RL_diagnostics(
  fit,
  encodedA,
  encodedB,
  compare_vars,
  vars_type_cont,
  true_pairs = NULL,
  FDP_estimation = TRUE,
  ...
)
```

## Arguments

- fit:

  List with either `Delta` (a data frame with columns `i`, `j`, `x`, as
  returned by [`StEM()`](stEM.md)) or the three elements `idxA`, `idxB`,
  `LinkScore` (e.g. the output of a link_with\_\* wrapper).

- encodedA:

  Data source used for linkage (same encoding and row order as passed to
  the linkage method).

- encodedB:

  Data source used for linkage (same encoding and row order as passed to
  the linkage method).

- compare_vars:

  Character vector, variables (PIVs or others) to compare between the
  linked subset and each source file.

- vars_type_cont:

  Named list or logical vector, one entry per `compare_vars`: `TRUE` if
  the variable is treated as continuous, `FALSE` if categorical (one SMD
  per level).

- true_pairs:

  Optional data frame with 2 columns of true (A, B) indices, when known,
  to report the realised FDP and sensitivity.

- FDP_estimation:

  Logical; if `TRUE`, run
  [`compute_RL_FDP_score()`](compute_RL_FDP_score.md) and
  [`compute_augmRL_FDP_synth()`](compute_augmRL_FDP_synth.md) to
  estimate the FDP over thresholds 0.50 to 0.95. `RL_method` and other
  arguments must then be given in `...`.

- ...:

  Arguments passed to
  [`compute_RL_FDP_score()`](compute_RL_FDP_score.md) and
  [`compute_augmRL_FDP_synth()`](compute_augmRL_FDP_synth.md) when
  `FDP_estimation = TRUE`: `RL_method`, `n_repeats`, `maxIter4CV`,
  optionally `synth_method` (default `"arf"`), `n_synth`, `PIVs`
  (default `compare_vars`), and the arguments of the chosen
  link_with\_\* wrapper.

## Value

An object of class `"RL_diagnostics"`, a list with:

- discrepancy_measures:

  data frame (class `"discrepancy_curves"`), one row per threshold `xi`:
  number of linked pairs, agreement rate per variable (see
  [`RL_agreement()`](rl_agreement.md)), standardised mean differences
  (one per continuous variable or per level, see [`smd()`](SMD.md)),
  support overlap per variable (see [`support_iou()`](support_iou.md))
  and multivariate discrepancy (see [`mmd()`](mmd.md)), each for the
  linked subset of A vs. A and of B vs. B

- FDP_measures:

  FDP estimates by threshold (class `"FDP_curves"`) if
  `FDP_estimation = TRUE`, else `NULL`

- true_agreement:

  only if `true_pairs` is supplied: agreement rate of the true pairs on
  each of `compare_vars`, to compare with the agreement of the linked
  pairs

- true_performance:

  only if `true_pairs` is supplied: data frame with the realised FDP and
  sensitivity by threshold

- idxA, idxB, LinkScore, RL_method, n_pairs, compare_vars, A, B:

  the linkage and data needed by the
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods

- gamma, eta, alpha, phi:

  the StEM chains from `fit`, if present

## Details

All measures are computed for every threshold from 0.50 to 0.99 (by
0.01) on the linkage scores; for a method without scores (or with a
single score value) they are computed once, for all returned pairs.
`print(x, threshold = )` summarises diagnostics at one threshold,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows the
curves.

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
#> Warning: PIV 'V4': values out of the common support in data2 [6].
#> Warning: PIVs 'V4' and 'V1' are associated (Cramer's V = 0.33 in data1, 0.23 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)  
PIVs_type <- list(V1=FALSE, V2=FALSE, V3=FALSE, V4=TRUE)

fit_flexrl <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2,
                    gibbs_iter = 5, gibbs_burnin = 2, n_post_samp = 10 )
#> FlexRL
diag_flexrl <- RL_diagnostics(fit_flexrl, prep_data$encodedA, prep_data$encodedB,
                              PIVs, PIVs_type, true_pairs = prep_data$true_pairs,
                              FDP_estimation = TRUE, RL_method = "FlexRL", 
                              data = prep_data,
                              StEM_iter = 5, StEM_burnin = 2, 
                              gibbs_iter = 5, gibbs_burnin = 2,
                              n_post_samp = 10,
                              maxIter4CV = 1, n_repeats = 1)
#> Iteration: 0, Accuracy: 50.65%
#> Iteration: 1, Accuracy: 39.84%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> FlexRL
#> FlexRL results (average over 1 iterations):
#>                            0.50  0.51  0.52  0.53  0.54  0.55  0.56  0.57  0.58
#> FDP model score estimator  0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14
#> FDP synth data estimator   0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14
#> Linked pairs (augm. RL)   73.00 73.00 73.00 73.00 73.00 73.00 73.00 73.00 73.00
#> Linked pairs (RL)         74.00 74.00 74.00 74.00 74.00 74.00 74.00 74.00 74.00
#>                            0.59  0.60  0.61  0.62  0.63  0.64  0.65  0.66  0.67
#> FDP model score estimator  0.14  0.10  0.10  0.10  0.10  0.10  0.10  0.10  0.10
#> FDP synth data estimator   0.14  0.16  0.16  0.16  0.16  0.16  0.16  0.16  0.16
#> Linked pairs (augm. RL)   73.00 63.00 63.00 63.00 63.00 63.00 63.00 63.00 63.00
#> Linked pairs (RL)         74.00 64.00 64.00 64.00 64.00 64.00 64.00 64.00 64.00
#>                            0.68  0.69  0.70  0.71  0.72  0.73  0.74  0.75  0.76
#> FDP model score estimator  0.10  0.10  0.06  0.06  0.06  0.06  0.06  0.06  0.06
#> FDP synth data estimator   0.16  0.16  0.00  0.00  0.00  0.00  0.00  0.00  0.00
#> Linked pairs (augm. RL)   63.00 63.00 52.00 52.00 52.00 52.00 52.00 52.00 52.00
#> Linked pairs (RL)         64.00 64.00 52.00 52.00 52.00 52.00 52.00 52.00 52.00
#>                            0.77  0.78  0.79  0.80  0.81  0.82  0.83  0.84  0.85
#> FDP model score estimator  0.06  0.06  0.06  0.03  0.03  0.03  0.03  0.03  0.03
#> FDP synth data estimator   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
#> Linked pairs (augm. RL)   52.00 52.00 52.00 44.00 44.00 44.00 44.00 44.00 44.00
#> Linked pairs (RL)         52.00 52.00 52.00 44.00 44.00 44.00 44.00 44.00 44.00
#>                            0.86  0.87  0.88  0.89 0.90 0.91 0.92 0.93 0.94 0.95
#> FDP model score estimator  0.03  0.03  0.03  0.03    0    0    0    0    0    0
#> FDP synth data estimator   0.00  0.00  0.00  0.00    0    0    0    0    0    0
#> Linked pairs (augm. RL)   44.00 44.00 44.00 44.00   31   31   31   31   31   31
#> Linked pairs (RL)         44.00 44.00 44.00 44.00   31   31   31   31   31   31
#>                           0.96 0.97 0.98 0.99
#> FDP model score estimator    0    0    0    0
#> FDP synth data estimator     0    0    0    0
#> Linked pairs (augm. RL)     31   31   31   31
#> Linked pairs (RL)           31   31   31   31
#> FlexRL
#> FlexRL results (average over 1 iterations):
#>                            0.50  0.51  0.52  0.53  0.54  0.55  0.56  0.57  0.58
#> FDP model score estimator  0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14
#> Linked pairs (RL)         87.00 87.00 87.00 87.00 87.00 87.00 87.00 87.00 87.00
#>                            0.59 0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68
#> FDP model score estimator  0.14  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
#> Linked pairs (RL)         87.00 76.0 76.0 76.0 76.0 76.0 76.0 76.0 76.0 76.0
#>                           0.69  0.70  0.71  0.72  0.73  0.74  0.75  0.76  0.77
#> FDP model score estimator  0.1  0.05  0.05  0.05  0.05  0.05  0.05  0.05  0.05
#> Linked pairs (RL)         76.0 62.00 62.00 62.00 62.00 62.00 62.00 62.00 62.00
#>                            0.78  0.79  0.80  0.81  0.82  0.83  0.84  0.85  0.86
#> FDP model score estimator  0.05  0.05  0.02  0.02  0.02  0.02  0.02  0.02  0.02
#> Linked pairs (RL)         62.00 62.00 52.00 52.00 52.00 52.00 52.00 52.00 52.00
#>                            0.87  0.88  0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96
#> FDP model score estimator  0.02  0.02  0.02    0    0    0    0    0    0    0
#> Linked pairs (RL)         52.00 52.00 52.00   40   40   40   40   40   40   40
#>                           0.97 0.98 0.99
#> FDP model score estimator    0    0    0
#> Linked pairs (RL)           40   40   40
diag_flexrl # print(diag_flexrl)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.50):  72
#> 
#>   FDP           (true pairs):    0.375
#>   Sensitivity   (true pairs):    0.464
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.136
#>                                                threshold 0.50: FDP ~ 0.136
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.142
#>                                                threshold 0.50: FDP ~ 0.142
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.140
#>                                                threshold 0.50: FDP ~ 0.140
#>                                            max threshold 0.99: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0115
#>   Multivariate MMD  (linked B vs. B):  0.0141
#> 
#>   IoU support V1  (linked A vs. A):  0.8000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  0.8571
#>   IoU support V4  (linked A vs. A):  1.0000
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  0.8182
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9722
#>   Agreement V2  (linked A vs. linked B):  0.9583
#>   Agreement V3  (linked A vs. linked B):  0.9444
#>   Agreement V4  (linked A vs. linked B):  0.8194
#>   Agreement V1  (true pairs):  0.9485
#>   Agreement V2  (true pairs):  0.9278
#>   Agreement V3  (true pairs):  0.8969
#>   Agreement V4  (true pairs):  0.7732
#> 
#>   SMD V1 values: {0, 2, 4, ...}  (linked A vs. A):  -0.1159, 0.0922, -0.0565, ...
#>   SMD V2 values: {2, 5, 6, ...}  (linked A vs. A):  0.1677, 0.0909, -0.2645, ...
#>   SMD V3 values: {0, 1, 6, ...}  (linked A vs. A):  -0.1159, 0.1562, -0.1222, ...
#>   SMD V4                         (linked A vs. A):  -0.0319
#>   SMD V1 values: {1, 2, 5, ...}  (linked B vs. B):  0.1038, 0.1699, -0.1400, ...
#>   SMD V2 values: {0, 1, 6, ...}  (linked B vs. B):  -0.0718, 0.1781, -0.2074, ...
#>   SMD V3 values: {1, 6, 7, ...}  (linked B vs. B):  0.2112, -0.1316, 0.1968, ...
#>   SMD V4                         (linked B vs. B):  0.0931
print(diag_flexrl, threshold = 0.75)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.75):  55
#> 
#>   FDP           (true pairs):    0.291
#>   Sensitivity   (true pairs):    0.402
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.136
#>                                                threshold 0.75: FDP ~ 0.052
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.142
#>                                                threshold 0.75: FDP ~ 0.056
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.140
#>                                                threshold 0.75: FDP ~ 0.000
#>                                            max threshold 0.99: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0173
#>   Multivariate MMD  (linked B vs. B):  0.0209
#> 
#>   IoU support V1  (linked A vs. A):  0.8000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  0.8571
#>   IoU support V4  (linked A vs. A):  0.8182
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  0.7273
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9818
#>   Agreement V2  (linked A vs. linked B):  0.9636
#>   Agreement V3  (linked A vs. linked B):  0.9455
#>   Agreement V4  (linked A vs. linked B):  0.8909
#>   Agreement V1  (true pairs):  0.9485
#>   Agreement V2  (true pairs):  0.9278
#>   Agreement V3  (true pairs):  0.8969
#>   Agreement V4  (true pairs):  0.7732
#> 
#>   SMD V1 values: {0, 2, 3, ...}  (linked A vs. A):  -0.1159, 0.1660, -0.0734, ...
#>   SMD V2 values: {1, 2, 6, ...}  (linked A vs. A):  0.1438, 0.1582, -0.2311, ...
#>   SMD V3 values: {0, 1, 7, ...}  (linked A vs. A):  -0.1159, 0.2060, -0.0760, ...
#>   SMD V4                         (linked A vs. A):  -0.0632
#>   SMD V1 values: {2, 3, 5, ...}  (linked B vs. B):  0.2084, -0.1830, -0.1984, ...
#>   SMD V2 values: {1, 2, 6, ...}  (linked B vs. B):  0.3223, 0.0867, -0.2045, ...
#>   SMD V3 values: {0, 1, 3, ...}  (linked B vs. B):  -0.1250, 0.2492, -0.0853, ...
#>   SMD V4                         (linked B vs. B):  -0.0643
plot(diag_flexrl, "scores")

plot(diag_flexrl, "distributions", threshold = 0.75)


plot(diag_flexrl, "convergence")










plot(diag_flexrl, "FDP")

plot(diag_flexrl, "discrepancy")





fit_brl <- link_with_BRL( prep_data$encodedA, prep_data$encodedB, 
                          list( flds = PIVs, 
                                types = rep("bi",length(PIVs)) ) )
diag_brl <- RL_diagnostics(fit_brl, prep_data$encodedA, prep_data$encodedB,
                           PIVs, PIVs_type, true_pairs = prep_data$true_pairs,
                           FDP_estimation = TRUE, RL_method = "BRL", 
                           flds = PIVs, types = rep("bi",length(PIVs)),
                           maxIter4CV = 1, n_repeats = 1)
#> Iteration: 0, Accuracy: 49.74%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> BRL results (average over 1 iterations):
#>                            0.50  0.51  0.52  0.53  0.54  0.55  0.56  0.57  0.58
#> FDP model score estimator  0.21  0.21  0.21  0.21  0.21  0.21  0.21  0.21  0.21
#> FDP synth data estimator   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
#> Linked pairs (augm. RL)   55.00 54.00 54.00 54.00 54.00 54.00 54.00 54.00 54.00
#> Linked pairs (RL)         55.00 54.00 54.00 54.00 54.00 54.00 54.00 54.00 54.00
#>                            0.59  0.60  0.61  0.62  0.63  0.64  0.65  0.66  0.67
#> FDP model score estimator  0.21  0.21  0.21  0.21  0.21  0.21  0.21  0.21  0.21
#> FDP synth data estimator   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
#> Linked pairs (augm. RL)   54.00 54.00 54.00 54.00 54.00 54.00 54.00 54.00 54.00
#> Linked pairs (RL)         54.00 54.00 54.00 54.00 54.00 54.00 54.00 54.00 54.00
#>                           0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76  0.77
#> FDP model score estimator  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.19
#> FDP synth data estimator   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.00
#> Linked pairs (augm. RL)   52.0 52.0 52.0 52.0 52.0 52.0 49.0 49.0 48.0 43.00
#> Linked pairs (RL)         52.0 52.0 52.0 52.0 52.0 52.0 49.0 49.0 48.0 43.00
#>                            0.78  0.79  0.80  0.81  0.82 0.83 0.84 0.85 0.86
#> FDP model score estimator  0.19  0.18  0.18  0.17  0.16 0.15 0.15 0.15    0
#> FDP synth data estimator   0.00  0.00  0.00  0.00  0.00 0.00 0.00 0.00    0
#> Linked pairs (augm. RL)   38.00 34.00 26.00 17.00 12.00 8.00 8.00 1.00    0
#> Linked pairs (RL)         38.00 34.00 26.00 17.00 12.00 8.00 8.00 1.00    0
#>                           0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> FDP synth data estimator     0    0    0    0    0    0    0    0    0    0
#> Linked pairs (augm. RL)      0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)            0    0    0    0    0    0    0    0    0    0
#>                           0.97 0.98 0.99
#> FDP model score estimator    0    0    0
#> FDP synth data estimator     0    0    0
#> Linked pairs (augm. RL)      0    0    0
#> Linked pairs (RL)            0    0    0
#> BRL results (average over 1 iterations):
#>                           0.50  0.51  0.52  0.53  0.54  0.55  0.56  0.57  0.58
#> FDP model score estimator  0.2  0.19  0.19  0.19  0.19  0.18  0.18  0.18  0.18
#> Linked pairs (RL)         55.0 54.00 54.00 54.00 54.00 53.00 53.00 53.00 53.00
#>                            0.59  0.60  0.61  0.62  0.63  0.64  0.65  0.66  0.67
#> FDP model score estimator  0.18  0.18  0.18  0.18  0.18  0.18  0.18  0.18  0.18
#> Linked pairs (RL)         53.00 53.00 53.00 53.00 53.00 53.00 53.00 53.00 53.00
#>                            0.68  0.69  0.70  0.71  0.72  0.73  0.74  0.75  0.76
#> FDP model score estimator  0.18  0.18  0.18  0.18  0.18  0.18  0.18  0.18  0.18
#> Linked pairs (RL)         53.00 53.00 52.00 51.00 51.00 51.00 51.00 51.00 50.00
#>                            0.77  0.78  0.79  0.80  0.81  0.82  0.83  0.84 0.85
#> FDP model score estimator  0.18  0.18  0.17  0.17  0.16  0.16  0.15  0.15 0.14
#> Linked pairs (RL)         50.00 46.00 43.00 38.00 31.00 26.00 18.00 13.00 8.00
#>                           0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95
#> FDP model score estimator 0.14    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)         3.00    0    0    0    0    0    0    0    0    0
#>                           0.96 0.97 0.98 0.99
#> FDP model score estimator    0    0    0    0
#> Linked pairs (RL)            0    0    0    0
diag_brl # print(diag_brl)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.50):  55
#> 
#>   FDP           (true pairs):    0.127
#>   Sensitivity   (true pairs):    0.495
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.195
#>                                                threshold 0.50: FDP ~ 0.195
#>                                            max threshold 0.85: FDP ~ 0.140
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.211
#>                                                threshold 0.50: FDP ~ 0.211
#>                                            max threshold 0.85: FDP ~ 0.149
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.000
#>                                                threshold 0.50: FDP ~ 0.000
#>                                            max threshold 0.85: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0183
#>   Multivariate MMD  (linked B vs. B):  0.0187
#> 
#>   IoU support V1  (linked A vs. A):  0.8000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  0.8571
#>   IoU support V4  (linked A vs. A):  0.7273
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  0.7273
#> 
#>   Agreement V1  (linked A vs. linked B):  1.0000
#>   Agreement V2  (linked A vs. linked B):  1.0000
#>   Agreement V3  (linked A vs. linked B):  1.0000
#>   Agreement V4  (linked A vs. linked B):  1.0000
#>   Agreement V1  (true pairs):  0.9485
#>   Agreement V2  (true pairs):  0.9278
#>   Agreement V3  (true pairs):  0.8969
#>   Agreement V4  (true pairs):  0.7732
#> 
#>   SMD V1 values: {3, 4, 5, ...}  (linked A vs. A):  0.1181, -0.1499, 0.1613, ...
#>   SMD V2 values: {1, 2, 5, ...}  (linked A vs. A):  0.0150, 0.1010, -0.1054, ...
#>   SMD V3 values: {0, 3, 5, ...}  (linked A vs. A):  -0.1159, 0.0969, 0.0665, ...
#>   SMD V4                         (linked A vs. A):  0.1380
#>   SMD V1 values: {0, 4, 5, ...}  (linked B vs. B):  -0.0718, -0.0478, 0.0330, ...
#>   SMD V2 values: {0, 1, 5, ...}  (linked B vs. B):  -0.0718, 0.0507, -0.1497, ...
#>   SMD V3 values: {0, 2, 7, ...}  (linked B vs. B):  -0.1250, -0.1013, 0.1210, ...
#>   SMD V4                         (linked B vs. B):  0.0403
print(diag_brl, threshold = 0.75)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.75):  51
#> 
#>   FDP           (true pairs):    0.137
#>   Sensitivity   (true pairs):    0.454
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.195
#>                                                threshold 0.75: FDP ~ 0.180
#>                                            max threshold 0.85: FDP ~ 0.140
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.211
#>                                                threshold 0.75: FDP ~ 0.197
#>                                            max threshold 0.85: FDP ~ 0.149
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.000
#>                                                threshold 0.75: FDP ~ 0.000
#>                                            max threshold 0.85: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0210
#>   Multivariate MMD  (linked B vs. B):  0.0216
#> 
#>   IoU support V1  (linked A vs. A):  0.8000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  0.8571
#>   IoU support V4  (linked A vs. A):  0.7273
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  0.7273
#> 
#>   Agreement V1  (linked A vs. linked B):  1.0000
#>   Agreement V2  (linked A vs. linked B):  1.0000
#>   Agreement V3  (linked A vs. linked B):  1.0000
#>   Agreement V4  (linked A vs. linked B):  1.0000
#>   Agreement V1  (true pairs):  0.9485
#>   Agreement V2  (true pairs):  0.9278
#>   Agreement V3  (true pairs):  0.8969
#>   Agreement V4  (true pairs):  0.7732
#> 
#>   SMD V1 values: {1, 3, 4, ...}  (linked A vs. A):  -0.1275, 0.1631, -0.1170, ...
#>   SMD V2 values: {2, 5, 6, ...}  (linked A vs. A):  0.1369, -0.0762, -0.1132, ...
#>   SMD V3 values: {0, 3, 5, ...}  (linked A vs. A):  -0.1159, 0.1254, 0.0988, ...
#>   SMD V4                         (linked A vs. A):  0.1527
#>   SMD V1 values: {0, 2, 5, ...}  (linked B vs. B):  -0.0718, 0.0454, -0.0343, ...
#>   SMD V2 values: {1, 5, 6, ...}  (linked B vs. B):  0.0774, -0.1216, -0.0848, ...
#>   SMD V3 values: {0, 1, 7, ...}  (linked B vs. B):  -0.1250, 0.0882, 0.1346, ...
#>   SMD V4                         (linked B vs. B):  0.0545
plot(diag_brl, "scores")

plot(diag_brl, "distributions", threshold = 0.75)


plot(diag_brl, "FDP")

plot(diag_brl, "discrepancy")



```
