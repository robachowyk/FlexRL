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
#>              0.5000              0.3375              0.5300 

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
#> FDP estimation on augmented RL task
#> Iteration: 0, Accuracy: 45.45%
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: 2 real records were dropped from file B during synthesis; the real/synthetic index split used by the FDP_synth estimator is invalid for data that was not support-restricted by prepare_data().
#> FlexRL
#> FlexRL results (average over 1 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold              0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator (mean) 0.16    0.16    0.08    0.08    0.08    0.08    0   0   0   0
#> FDP model score estimator (sd)   NA  NA  NA  NA  NA  NA  NA  NA  NA  NA
#> FDP synth data estimator (mean)  0.14    0.14    0.19    0.19    0.19    0.19    0   0   0   0
#> FDP synth data estimator (sd)    NA  NA  NA  NA  NA  NA  NA  NA  NA  NA
#> Linked pairs (augm. RL)          78  78  57  57  57  57  35  35  35  35
#> Linked pairs (RL)                79  79  58  58  58  58  35  35  35  35
#> 
#> FDP estimation on original RL task
#> FlexRL
#> FlexRL results (average over 1 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold              0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator (mean) 0.15    0.15    0.04    0.04    0.04    0.04    0   0   0   0
#> FDP model score estimator (sd)   NA  NA  NA  NA  NA  NA  NA  NA  NA  NA
#> Linked pairs (RL)                85  85  60  60  60  60  48  48  48  48
diag # print(diag)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.50):  80
#> 
#>   FDP           (true pairs):    0.338
#>   Sensitivity   (true pairs):    0.530
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.146 (sd NA)
#>                                                threshold 0.50: FDP ~ 0.146 (sd NA)
#>                                            max threshold 0.99: FDP ~ 0.000 (sd NA)
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.165 (sd NA)
#>                                                threshold 0.50: FDP ~ 0.165 (sd NA)
#>                                            max threshold 0.99: FDP ~ 0.000 (sd NA)
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.142 (sd NA)
#>                                                threshold 0.50: FDP ~ 0.142 (sd NA)
#>                                            max threshold 0.99: FDP ~ 0.000 (sd NA)
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0073
#>   Multivariate MMD  (linked B vs. B):  0.0098
#> 
#>   IoU support V1  (linked A vs. A):  0.8000
#>   IoU support V2  (linked A vs. A):  0.8333
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  1.0000
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  1.0000
#>   IoU support V4  (linked B vs. B):  1.0000
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9625
#>   Agreement V2  (linked A vs. linked B):  0.9125
#>   Agreement V3  (linked A vs. linked B):  0.9000
#>   Agreement V4  (linked A vs. linked B):  0.8750
#>   Agreement V1  (true pairs):  0.9200
#>   Agreement V2  (true pairs):  0.9400
#>   Agreement V3  (true pairs):  0.9100
#>   Agreement V4  (true pairs):  0.7900
#> 
#>   SMD V1 values: {0, 1, 2, ...}  (linked A vs. A):  -0.0816, 0.0699, 0.0446, ...
#>   SMD V2 values: {0, 4, 5, ...}  (linked A vs. A):  -0.1159, 0.1226, -0.1087, ...
#>   SMD V3 values: {2, 4, 6, ...}  (linked A vs. A):  0.1194, -0.0997, -0.1079, ...
#>   SMD V4                         (linked A vs. A):  0.0324
#>   SMD V1 values: {2, 4, 5, ...}  (linked B vs. B):  0.1485, 0.0755, -0.1731, ...
#>   SMD V2 values: {0, 1, 5, ...}  (linked B vs. B):  -0.0707, 0.0984, -0.1133, ...
#>   SMD V3 values: {3, 4, 7, ...}  (linked B vs. B):  0.1021, 0.0651, -0.0666, ...
#>   SMD V4                         (linked B vs. B):  -0.1509
print(diag, threshold = 0.75)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.75):  66
#> 
#>   FDP           (true pairs):    0.288
#>   Sensitivity   (true pairs):    0.470
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.146 (sd NA)
#>                                                threshold 0.75: FDP ~ 0.040 (sd NA)
#>                                            max threshold 0.99: FDP ~ 0.000 (sd NA)
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.165 (sd NA)
#>                                                threshold 0.75: FDP ~ 0.079 (sd NA)
#>                                            max threshold 0.99: FDP ~ 0.000 (sd NA)
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.142 (sd NA)
#>                                                threshold 0.75: FDP ~ 0.195 (sd NA)
#>                                            max threshold 0.99: FDP ~ 0.000 (sd NA)
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0098
#>   Multivariate MMD  (linked B vs. B):  0.0124
#> 
#>   IoU support V1  (linked A vs. A):  0.8000
#>   IoU support V2  (linked A vs. A):  0.8333
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  1.0000
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  1.0000
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9848
#>   Agreement V2  (linked A vs. linked B):  0.9242
#>   Agreement V3  (linked A vs. linked B):  0.9242
#>   Agreement V4  (linked A vs. linked B):  0.8939
#>   Agreement V1  (true pairs):  0.9200
#>   Agreement V2  (true pairs):  0.9400
#>   Agreement V3  (true pairs):  0.9100
#>   Agreement V4  (true pairs):  0.7900
#> 
#>   SMD V1 values: {0, 1, 5, ...}  (linked A vs. A):  -0.0816, 0.1199, -0.0691, ...
#>   SMD V2 values: {0, 2, 3, ...}  (linked A vs. A):  -0.1159, 0.0969, -0.1308, ...
#>   SMD V3 values: {1, 2, 6, ...}  (linked A vs. A):  -0.0923, 0.1514, -0.0799, ...
#>   SMD V4                         (linked A vs. A):  0.0794
#>   SMD V1 values: {1, 2, 5, ...}  (linked B vs. B):  0.0924, 0.1615, -0.2152, ...
#>   SMD V2 values: {2, 3, 6, ...}  (linked B vs. B):  0.0753, -0.1055, 0.0838, ...
#>   SMD V3 values: {0, 4, 7, ...}  (linked B vs. B):  -0.1003, 0.1021, -0.1171, ...
#>   SMD V4                         (linked B vs. B):  -0.0711
plot(diag, "scores")

plot(diag, "distributions", threshold = 0.75)


plot(diag, "convergence")










plot(diag, "FDP")

plot(diag, "discrepancy")



```
