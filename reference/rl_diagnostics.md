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
  PIVs,
  true_pairs = NULL,
  FDP_estimation = TRUE,
  ...
)
```

## Arguments

- fit:

  List with either `Delta` (a data frame with columns `i`, `j`, `x`, as
  returned by
  [`StEM()`](https://robachowyk.github.io/FlexRL/reference/StEM.md)) or
  the three elements `idxA`, `idxB`, `LinkScore` (e.g. the output of a
  link_with\_\* wrapper).

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

- PIVs:

  Character vector, partially identifying variables to compare between
  the linked records from A and from B.

- true_pairs:

  Optional data frame with 2 columns of true (A, B) indices, when known,
  to report the realised FDP and sensitivity.

- FDP_estimation:

  Logical; if `TRUE`, run
  [`compute_RL_FDP_score()`](https://robachowyk.github.io/FlexRL/reference/compute_RL_FDP_score.md)
  and
  [`compute_augmRL_FDP_synth()`](https://robachowyk.github.io/FlexRL/reference/compute_augmRL_FDP_synth.md)
  to estimate the FDP over thresholds 0.50 to 0.95. `RL_method` and
  other arguments must then be given in `...`.

- ...:

  Arguments passed to
  [`compute_RL_FDP_score()`](https://robachowyk.github.io/FlexRL/reference/compute_RL_FDP_score.md)
  and
  [`compute_augmRL_FDP_synth()`](https://robachowyk.github.io/FlexRL/reference/compute_augmRL_FDP_synth.md)
  when `FDP_estimation = TRUE`: `RL_method`, `n_repeats`, `maxIter4CV`,
  optionally `synth_method` (default `"arf"`), `n_synth`, `PIVs`
  (default `compare_vars`), and the arguments of the chosen
  link_with\_\* wrapper.

## Value

An object of class `"RL_diagnostics"`, a list with:

- discrepancy_measures:

  data frame (class `"discrepancy_curves"`), one row per threshold `xi`:
  number of linked pairs, agreement rate per variable (see
  [`RL_agreement()`](https://robachowyk.github.io/FlexRL/reference/RL_agreement.md)),
  standardised mean differences (one per continuous variable or per
  level, see
  [`smd()`](https://robachowyk.github.io/FlexRL/reference/smd.md)),
  support overlap per variable (see
  [`support_iou()`](https://robachowyk.github.io/FlexRL/reference/support_iou.md))
  and multivariate discrepancy (see
  [`mmd()`](https://robachowyk.github.io/FlexRL/reference/mmd.md)), each
  for the linked subset of A vs. A and of B vs. B

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
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)  
PIVs_type <- list(V1=FALSE, V2=FALSE, V3=FALSE, V4=TRUE)
fit_brl <- link_with_BRL( prep_data$encodedA, prep_data$encodedB, 
                          list( flds = PIVs, 
                                types = rep("bi",length(PIVs)),
                                nIter = 500 ) )
diag_brl <- RL_diagnostics(fit_brl, prep_data$encodedA, prep_data$encodedB,
                           PIVs, PIVs_type, PIVs, true_pairs = prep_data$true_pairs,
                           FDP_estimation = TRUE, RL_method = "BRL", 
                           flds = PIVs, types = rep("bi",length(PIVs)), nIter = 500,
                           maxIter4CV = 1, n_repeats = 1)
#> 
#> FDP estimation on augmented RL taskIteration: 0, Accuracy: 47.69%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> Warning: 2 real records were dropped from file B during synthesis; the real/synthetic index split used by the FDP_synth estimator is invalid for data that was not support-restricted by prepare_data().
#> BRL results (average over 1 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.17    0.16    0.15    0.15    0.14    0.14    0.13    0.11    0.08    0
#> FDP synth data estimator     0.35    0.18    0.18    0.18    0.19    0.2 0.22    0   0   0
#> Linked pairs (augm. RL)      60  58  57  57  54  53  47  32  11  0
#> Linked pairs (RL)            62  59  58  58  55  54  48  32  11  0
#> 
#> FDP estimation on original RL task
#> BRL results (average over 1 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.16    0.15    0.14    0.14    0.14    0.13    0.12    0.1 0.07    0.05
#> Linked pairs (RL)            64  62  59  59  59  55  49  38  18  1
diag_brl # print(diag_brl)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.50):  61
#> 
#>   FDP           (true pairs):    0.115
#>   Sensitivity   (true pairs):    0.540
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.161
#>                                                threshold 0.50: FDP ~ 0.161
#>                                            max threshold 0.94: FDP ~ 0.053
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.174
#>                                                threshold 0.50: FDP ~ 0.174
#>                                            max threshold 0.94: FDP ~ 0.054
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.351
#>                                                threshold 0.50: FDP ~ 0.351
#>                                            max threshold 0.94: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0155
#>   Multivariate MMD  (linked B vs. B):  0.0176
#> 
#>   IoU support V1  (linked A vs. A):  0.8000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  0.8571
#>   IoU support V4  (linked A vs. A):  0.9167
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  1.0000
#>   IoU support V4  (linked B vs. B):  1.0000
#> 
#>   Agreement V1  (linked A vs. linked B):  1.0000
#>   Agreement V2  (linked A vs. linked B):  1.0000
#>   Agreement V3  (linked A vs. linked B):  1.0000
#>   Agreement V4  (linked A vs. linked B):  1.0000
#>   Agreement V1  (true pairs):  0.9300
#>   Agreement V2  (true pairs):  0.9700
#>   Agreement V3  (true pairs):  0.8900
#>   Agreement V4  (true pairs):  0.7500
#> 
#>   SMD V1 values: {1, 2, 4, ...}  (linked A vs. A):  -0.1166, -0.1827, 0.1961, ...
#>   SMD V2 values: {2, 4, 6, ...}  (linked A vs. A):  0.0995, -0.0950, -0.0613, ...
#>   SMD V3 values: {1, 3, 6, ...}  (linked A vs. A):  0.1576, 0.1035, -0.1569, ...
#>   SMD V4                         (linked A vs. A):  -0.0347
#>   SMD V1 values: {1, 2, 4, ...}  (linked B vs. B):  -0.2123, -0.1902, 0.1682, ...
#>   SMD V2 values: {2, 4, 5, ...}  (linked B vs. B):  0.1064, -0.1268, 0.0750, ...
#>   SMD V3 values: {1, 6, 7, ...}  (linked B vs. B):  0.1881, -0.0910, -0.1227, ...
#>   SMD V4                         (linked B vs. B):  -0.0562
print(diag_brl, threshold = 0.75)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.75):  58
#> 
#>   FDP           (true pairs):    0.121
#>   Sensitivity   (true pairs):    0.510
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.161
#>                                                threshold 0.75: FDP ~ 0.127
#>                                            max threshold 0.94: FDP ~ 0.053
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.174
#>                                                threshold 0.75: FDP ~ 0.142
#>                                            max threshold 0.94: FDP ~ 0.054
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.351
#>                                                threshold 0.75: FDP ~ 0.199
#>                                            max threshold 0.94: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0154
#>   Multivariate MMD  (linked B vs. B):  0.0173
#> 
#>   IoU support V1  (linked A vs. A):  0.8000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  0.8571
#>   IoU support V4  (linked A vs. A):  0.9167
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  1.0000
#>   IoU support V4  (linked B vs. B):  1.0000
#> 
#>   Agreement V1  (linked A vs. linked B):  1.0000
#>   Agreement V2  (linked A vs. linked B):  1.0000
#>   Agreement V3  (linked A vs. linked B):  1.0000
#>   Agreement V4  (linked A vs. linked B):  1.0000
#>   Agreement V1  (true pairs):  0.9300
#>   Agreement V2  (true pairs):  0.9700
#>   Agreement V3  (true pairs):  0.8900
#>   Agreement V4  (true pairs):  0.7500
#> 
#>   SMD V1 values: {0, 2, 4, ...}  (linked A vs. A):  -0.1159, -0.1691, 0.1890, ...
#>   SMD V2 values: {2, 4, 6, ...}  (linked A vs. A):  0.1205, -0.0768, -0.0715, ...
#>   SMD V3 values: {1, 3, 6, ...}  (linked A vs. A):  0.1816, 0.1260, -0.1833, ...
#>   SMD V4                         (linked A vs. A):  -0.0701
#>   SMD V1 values: {1, 2, 4, ...}  (linked B vs. B):  -0.2008, -0.1767, 0.1613, ...
#>   SMD V2 values: {0, 2, 4, ...}  (linked B vs. B):  -0.0707, 0.1277, -0.1092, ...
#>   SMD V3 values: {1, 6, 7, ...}  (linked B vs. B):  0.2130, -0.1191, -0.1390, ...
#>   SMD V4                         (linked B vs. B):  -0.0907
plot(diag_brl, "scores")

plot(diag_brl, "distributions", threshold = 0.75)


plot(diag_brl, "FDP")

plot(diag_brl, "discrepancy")



```
