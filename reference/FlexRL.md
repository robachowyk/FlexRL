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
#> Warning: PIVs 'V4' and 'V3' are associated (Cramer's V = 0.31 in data1, 0.27 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
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
#>           0.5000000           0.2727273           0.6400000 

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
#> FDP estimation on augmented RL taskIteration: 0, Accuracy: 53.52%
#> Iteration: 1, Accuracy: 40%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> FlexRL results (average over 1 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.11    0.11    0.03    0.03    0.03    0.03    0   0   0   0
#> FDP synth data estimator     0.44    0.44    0.42    0.42    0.42    0.42    0.49    0.49    0.49    0.49
#> Linked pairs (augm. RL)      90  90  71  71  71  71  61  61  61  61
#> Linked pairs (RL)            94  94  74  74  74  74  64  64  64  64
#> 
#> FDP estimation on original RL task
#> FlexRL
#> FlexRL results (average over 1 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.12    0.12    0.05    0.05    0.05    0.05    0   0   0   0
#> Linked pairs (RL)            105 105 84  84  84  84  65  65  65  65
diag # print(diag)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.50):  88
#> 
#>   FDP           (true pairs):    0.273
#>   Sensitivity   (true pairs):    0.640
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.116
#>                                                threshold 0.50: FDP ~ 0.116
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.106
#>                                                threshold 0.50: FDP ~ 0.106
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.444
#>                                                threshold 0.50: FDP ~ 0.444
#>                                            max threshold 0.99: FDP ~ 0.492
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0069
#>   Multivariate MMD  (linked B vs. B):  0.0081
#> 
#>   IoU support V1  (linked A vs. A):  1.0000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  0.9167
#>   IoU support V1  (linked B vs. B):  1.0000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  1.0000
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9205
#>   Agreement V2  (linked A vs. linked B):  0.9886
#>   Agreement V3  (linked A vs. linked B):  1.0000
#>   Agreement V4  (linked A vs. linked B):  0.8409
#>   Agreement V1  (true pairs):  0.9500
#>   Agreement V2  (true pairs):  0.9700
#>   Agreement V3  (true pairs):  0.9300
#>   Agreement V4  (true pairs):  0.7200
#> 
#>   SMD V1 values: {3, 4, 5, ...}  (linked A vs. A):  0.0709, 0.1003, -0.1197, ...
#>   SMD V2 values: {0, 2, 5, ...}  (linked A vs. A):  0.0575, 0.0592, -0.1132, ...
#>   SMD V3 values: {4, 6, 7, ...}  (linked A vs. A):  -0.0731, 0.1013, -0.0242, ...
#>   SMD V4                         (linked A vs. A):  -0.0526
#>   SMD V1 values: {3, 4, 5, ...}  (linked B vs. B):  -0.0208, 0.0912, -0.0812, ...
#>   SMD V2 values: {0, 2, 5, ...}  (linked B vs. B):  -0.1003, 0.0831, -0.1476, ...
#>   SMD V3 values: {0, 1, 7, ...}  (linked B vs. B):  -0.0707, -0.0760, 0.1551, ...
#>   SMD V4                         (linked B vs. B):  0.0250
print(diag, threshold = 0.75)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.75):  74
#> 
#>   FDP           (true pairs):    0.243
#>   Sensitivity   (true pairs):    0.560
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.116
#>                                                threshold 0.75: FDP ~ 0.045
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.106
#>                                                threshold 0.75: FDP ~ 0.027
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.444
#>                                                threshold 0.75: FDP ~ 0.423
#>                                            max threshold 0.99: FDP ~ 0.492
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0088
#>   Multivariate MMD  (linked B vs. B):  0.0084
#> 
#>   IoU support V1  (linked A vs. A):  1.0000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  0.9167
#>   IoU support V1  (linked B vs. B):  1.0000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  1.0000
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9189
#>   Agreement V2  (linked A vs. linked B):  0.9865
#>   Agreement V3  (linked A vs. linked B):  1.0000
#>   Agreement V4  (linked A vs. linked B):  0.8919
#>   Agreement V1  (true pairs):  0.9500
#>   Agreement V2  (true pairs):  0.9700
#>   Agreement V3  (true pairs):  0.9300
#>   Agreement V4  (true pairs):  0.7200
#> 
#>   SMD V1 values: {0, 3, 5, ...}  (linked A vs. A):  0.0839, 0.1325, -0.1392, ...
#>   SMD V2 values: {0, 4, 5, ...}  (linked A vs. A):  0.0839, 0.0559, -0.0813, ...
#>   SMD V3 values: {3, 5, 6, ...}  (linked A vs. A):  0.0879, -0.1625, 0.1687, ...
#>   SMD V4                         (linked A vs. A):  -0.0617
#>   SMD V1 values: {2, 4, 5, ...}  (linked B vs. B):  -0.0365, 0.0600, -0.0709, ...
#>   SMD V2 values: {0, 4, 5, ...}  (linked B vs. B):  -0.1003, 0.0945, -0.1167, ...
#>   SMD V3 values: {2, 5, 7, ...}  (linked B vs. B):  -0.0990, -0.2002, 0.1915, ...
#>   SMD V4                         (linked B vs. B):  -0.0338
plot(diag, "scores")

plot(diag, "distributions", threshold = 0.75)


plot(diag, "convergence")










plot(diag, "FDP")

plot(diag, "discrepancy")



```
