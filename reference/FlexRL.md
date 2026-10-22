# FlexRL: A Flexible Model for Record Linkage

Links records that refer to the same entities across two data sources
without a unique identifier, using partially identifying variables.
These are stable (not changing over time, probability of mistakes could
be bounded), flexible (dynamic but no information to model changes over
time) or structured (dynamic and information to model changes over time,
probability of mistakes could be fixed). The main function
[`StEM()`](https://robachowyk.github.io/FlexRL/reference/StEM.md) fits a
latent-variable model by stochastic expectation-maximisation; it models
registration errors (missing values and mistakes) and changes over time.
[`prepare_data()`](https://robachowyk.github.io/FlexRL/reference/prepare_data.md)
prepares the data sources for record linkage,
[`RL_diagnostics()`](https://robachowyk.github.io/FlexRL/reference/RL_diagnostics.md)
gathers diagnostics (FDP estimation and discrepancy metrics for
inference on the linked data.

## Details

Methodological paper:
[doi:10.1093/jrsssc/qlaf016](https://doi.org/10.1093/jrsssc/qlaf016) .
False discovery proportion estimation:
[doi:10.1002/sim.70292](https://doi.org/10.1002/sim.70292) . Experiments
repository: <https://github.com/robachowyk/FlexRL-experiments>.

## See also

Useful links:

- <https://github.com/robachowyk/FlexRL>

- <https://robachowyk.github.io/FlexRL/>

- Report bugs at <https://github.com/robachowyk/FlexRL/issues>

## Author

Kayané Robach

## Examples

``` r
# Link two simulated sources with 4 PIVs: two stable, one flexible, one
# structured. The true links are known, so performance can be computed.
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
cond_hazard_params <- list(V1 = c(), V2 = c(),
                           V3 = c(), V4 = log(c(0.7, 0.6, 0.5)))
gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100,
                           p_mistake, p_missing, cond_hazard_params,
                           TRUE, survival_model("exponential") )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> Warning: PIV 'V4': values out of the common support in data2 [6].
#> '2' is the larger source, saved as B; '1' saved as A.
fit <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2,
             gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5 )
#> FlexRL

# linked pairs and performance against the true pairs
linked <- fit$Delta[fit$Delta$x > 0.5, ]
linked_pairs <- paste(linked$i, linked$j, sep = "_")
true_pairs   <- paste(prep_data$true_pairs[[1]], prep_data$true_pairs[[2]], sep = "_")
tp <- length(intersect(linked_pairs, true_pairs))
fp <- length(setdiff(linked_pairs, true_pairs))
fn <- length(setdiff(true_pairs, linked_pairs))
c(LinkageDecisionRule = 0.5, FDP = fp / (tp + fp), Sensitivity = tp / (tp + fn))
#> LinkageDecisionRule                 FDP         Sensitivity 
#>           0.5000000           0.3580247           0.5306122 

# diagnostics for inference on the linked data
diag <- RL_diagnostics(fit, prep_data$encodedA, prep_data$encodedB,
                       names(PIVs_config),
                       list(V1 = FALSE, V2 = FALSE, V3 = FALSE, V4 = TRUE),
                       names(PIVs_config),
                       true_pairs = prep_data$true_pairs, FDP_estimation = TRUE,
                       RL_method = "FlexRL", data = prep_data,
                       StEM_iter = 5, StEM_burnin = 2,
                       gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5,
                       maxIter4CV = 1, n_repeats = 1)
#> 
#> FDP estimation on augmented RL taskIteration: 0, Accuracy: 45.22%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> Warning: 3 real records were dropped from file B during synthesis; the real/synthetic index split used by the FDP_synth estimator is invalid for data that was not support-restricted by prepare_data().
#> FlexRL
#> FlexRL results (average over 1 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.16    0.16    0.06    0.06    0.06    0.06    0   0   0   0
#> FDP synth data estimator     0.58    0.58    0.63    0.63    0.63    0.63    0.59    0.59    0.59    0.59
#> Linked pairs (augm. RL)      80  80  55  55  55  55  39  39  39  39
#> Linked pairs (RL)            84  84  58  58  58  58  41  41  41  41
#> 
#> FDP estimation on original RL task
#> FlexRL
#> FlexRL results (average over 1 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.1 0.1 0.05    0.05    0.05    0.05    0   0   0   0
#> Linked pairs (RL)            95  95  81  81  81  81  61  61  61  61
diag # print(diag)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.50):  81
#> 
#>   FDP           (true pairs):    0.358
#>   Sensitivity   (true pairs):    0.531
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.101
#>                                                threshold 0.50: FDP ~ 0.101
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.164
#>                                                threshold 0.50: FDP ~ 0.164
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.576
#>                                                threshold 0.50: FDP ~ 0.576
#>                                            max threshold 0.99: FDP ~ 0.591
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0093
#>   Multivariate MMD  (linked B vs. B):  0.0102
#> 
#>   IoU support V1  (linked A vs. A):  1.0000
#>   IoU support V2  (linked A vs. A):  0.8333
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  1.0000
#>   IoU support V1  (linked B vs. B):  1.0000
#>   IoU support V2  (linked B vs. B):  1.0000
#>   IoU support V3  (linked B vs. B):  1.0000
#>   IoU support V4  (linked B vs. B):  0.9091
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9259
#>   Agreement V2  (linked A vs. linked B):  0.9630
#>   Agreement V3  (linked A vs. linked B):  0.9259
#>   Agreement V4  (linked A vs. linked B):  0.6420
#>   Agreement V1  (true pairs):  0.9694
#>   Agreement V2  (true pairs):  0.9592
#>   Agreement V3  (true pairs):  0.8776
#>   Agreement V4  (true pairs):  0.7449
#> 
#>   SMD V1 values: {1, 2, 4, ...}  (linked A vs. A):  0.1483, -0.0400, -0.0844, ...
#>   SMD V2 values: {0, 2, 3, ...}  (linked A vs. A):  -0.1159, 0.1303, -0.0836, ...
#>   SMD V3 values: {2, 4, 7, ...}  (linked A vs. A):  -0.0671, 0.1298, -0.0829, ...
#>   SMD V4                         (linked A vs. A):  -0.0951
#>   SMD V1 values: {1, 2, 3, ...}  (linked B vs. B):  0.0658, 0.0592, -0.0539, ...
#>   SMD V2 values: {1, 4, 5, ...}  (linked B vs. B):  0.0535, 0.0494, -0.0468, ...
#>   SMD V3 values: {1, 3, 7, ...}  (linked B vs. B):  -0.0886, 0.2224, -0.1493, ...
#>   SMD V4                         (linked B vs. B):  -0.0101
print(diag, threshold = 0.75)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.75):  67
#> 
#>   FDP           (true pairs):    0.269
#>   Sensitivity   (true pairs):    0.500
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.101
#>                                                threshold 0.75: FDP ~ 0.049
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.164
#>                                                threshold 0.75: FDP ~ 0.059
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.576
#>                                                threshold 0.75: FDP ~ 0.629
#>                                            max threshold 0.99: FDP ~ 0.591
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0140
#>   Multivariate MMD  (linked B vs. B):  0.0140
#> 
#>   IoU support V1  (linked A vs. A):  1.0000
#>   IoU support V2  (linked A vs. A):  0.8333
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  1.0000
#>   IoU support V1  (linked B vs. B):  1.0000
#>   IoU support V2  (linked B vs. B):  1.0000
#>   IoU support V3  (linked B vs. B):  1.0000
#>   IoU support V4  (linked B vs. B):  0.9091
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9403
#>   Agreement V2  (linked A vs. linked B):  0.9701
#>   Agreement V3  (linked A vs. linked B):  0.9254
#>   Agreement V4  (linked A vs. linked B):  0.6866
#>   Agreement V1  (true pairs):  0.9694
#>   Agreement V2  (true pairs):  0.9592
#>   Agreement V3  (true pairs):  0.8776
#>   Agreement V4  (true pairs):  0.7449
#> 
#>   SMD V1 values: {1, 2, 5, ...}  (linked A vs. A):  0.1129, -0.0422, -0.0482, ...
#>   SMD V2 values: {0, 3, 4, ...}  (linked A vs. A):  -0.1159, -0.0665, 0.0925, ...
#>   SMD V3 values: {4, 6, 7, ...}  (linked A vs. A):  0.1552, -0.0682, -0.1112, ...
#>   SMD V4                         (linked A vs. A):  -0.1770
#>   SMD V1 values: {1, 4, 5, ...}  (linked B vs. B):  0.0765, 0.0573, -0.1142, ...
#>   SMD V2 values: {1, 4, 6, ...}  (linked B vs. B):  0.1179, 0.0791, -0.0808, ...
#>   SMD V3 values: {1, 3, 7, ...}  (linked B vs. B):  -0.0958, 0.2729, -0.1977, ...
#>   SMD V4                         (linked B vs. B):  0.0178
plot(diag, "scores")

plot(diag, "distributions", threshold = 0.75)


plot(diag, "convergence")










plot(diag, "FDP")

plot(diag, "discrepancy")



```
