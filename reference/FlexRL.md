# FlexRL: A Flexible Model for Record Linkage

Links records that refer to the same entities across two data sources
without a unique identifier, using partially identifying variables.
These are stable (not changing over time, probability of mistakes could
be bounded), flexible (dynamic but no information to model changes over
time) or structured (dynamic and information to model changes over time,
probability of mistakes could be fixed). The main function
[`StEM()`](https://robachowyk.github.io/FlexRL/reference/stEM.md) fits a
latent-variable model by stochastic expectation-maximisation; it models
registration errors (missing values and mistakes) and changes over time.
[`prepare_data()`](https://robachowyk.github.io/FlexRL/reference/prepare_data.md)
prepares the data sources for record linkage,
[`RL_diagnostics()`](https://robachowyk.github.io/FlexRL/reference/rl_diagnostics.md)
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
#> Warning: PIVs 'V4' and 'V3' are associated (Cramer's V = 0.32 in data1, 0.23 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
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
#>           0.5000000           0.4556962           0.4300000 

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
#> FDP estimation on augmented RL taskIteration: 0, Accuracy: 48.49%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V3': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> FlexRL results (average over 1 iterations) (augmented RL task) `compute_augmRL_FDP_synth`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.13    0.13    0.05    0.05    0.05    0.05    0   0   0   0
#> FDP synth data estimator     0.34    0.34    0.14    0.14    0.14    0.14    0   0   0   0
#> Linked pairs (augm. RL)      88  88  69  69  69  69  54  54  54  54
#> Linked pairs (RL)            91  91  70  70  70  70  54  54  54  54
#> 
#> FDP estimation on original RL task
#> FlexRL
#> FlexRL results (average over 1 iterations) (original RL task) `compute_RL_FDP_score`
#> Decision rule threshold      0.5 0.55    0.6 0.65    0.7 0.75    0.8 0.85    0.9 0.95
#> FDP model score estimator    0.12    0.12    0.05    0.05    0.05    0.05    0   0   0   0
#> Linked pairs (RL)            89  89  72  72  72  72  54  54  54  54
diag # print(diag)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.50):  79
#> 
#>   FDP           (true pairs):    0.456
#>   Sensitivity   (true pairs):    0.430
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.117
#>                                                threshold 0.50: FDP ~ 0.117
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.127
#>                                                threshold 0.50: FDP ~ 0.127
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.341
#>                                                threshold 0.50: FDP ~ 0.341
#>                                            max threshold 0.99: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0069
#>   Multivariate MMD  (linked B vs. B):  0.0124
#> 
#>   IoU support V1  (linked A vs. A):  1.0000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  0.9091
#>   IoU support V1  (linked B vs. B):  1.0000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  1.0000
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9241
#>   Agreement V2  (linked A vs. linked B):  0.9494
#>   Agreement V3  (linked A vs. linked B):  0.8987
#>   Agreement V4  (linked A vs. linked B):  0.6582
#>   Agreement V1  (true pairs):  0.9600
#>   Agreement V2  (true pairs):  0.9500
#>   Agreement V3  (true pairs):  0.8700
#>   Agreement V4  (true pairs):  0.7300
#> 
#>   SMD V1 values: {0, 3, 4, ...}  (linked A vs. A):  0.0734, -0.2141, 0.0708, ...
#>   SMD V2 values: {1, 2, 5, ...}  (linked A vs. A):  -0.0963, 0.0915, 0.0769, ...
#>   SMD V3 values: {2, 4, 5, ...}  (linked A vs. A):  -0.0747, -0.1090, 0.1298, ...
#>   SMD V4                         (linked A vs. A):  -0.0090
#>   SMD V1 values: {2, 3, 4, ...}  (linked B vs. B):  0.1055, -0.1598, 0.0364, ...
#>   SMD V2 values: {0, 2, 4, ...}  (linked B vs. B):  -0.0707, 0.0591, -0.1008, ...
#>   SMD V3 values: {1, 6, 7, ...}  (linked B vs. B):  0.1186, 0.1326, -0.1088, ...
#>   SMD V4                         (linked B vs. B):  -0.0748
print(diag, threshold = 0.75)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.75):  54
#> 
#>   FDP           (true pairs):    0.389
#>   Sensitivity   (true pairs):    0.330
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.117
#>                                                threshold 0.75: FDP ~ 0.050
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.127
#>                                                threshold 0.75: FDP ~ 0.046
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.341
#>                                                threshold 0.75: FDP ~ 0.145
#>                                            max threshold 0.99: FDP ~ 0.000
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0122
#>   Multivariate MMD  (linked B vs. B):  0.0213
#> 
#>   IoU support V1  (linked A vs. A):  1.0000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  0.9091
#>   IoU support V1  (linked B vs. B):  1.0000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  0.9167
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9074
#>   Agreement V2  (linked A vs. linked B):  0.9444
#>   Agreement V3  (linked A vs. linked B):  0.9074
#>   Agreement V4  (linked A vs. linked B):  0.7222
#>   Agreement V1  (true pairs):  0.9600
#>   Agreement V2  (true pairs):  0.9500
#>   Agreement V3  (true pairs):  0.8700
#>   Agreement V4  (true pairs):  0.7300
#> 
#>   SMD V1 values: {0, 3, 4, ...}  (linked A vs. A):  0.1452, -0.2192, 0.1245, ...
#>   SMD V2 values: {3, 5, 6, ...}  (linked A vs. A):  0.0881, 0.0800, -0.1425, ...
#>   SMD V3 values: {3, 4, 5, ...}  (linked A vs. A):  0.1339, -0.2192, 0.0990, ...
#>   SMD V4                         (linked A vs. A):  -0.0869
#>   SMD V1 values: {2, 3, 4, ...}  (linked B vs. B):  0.0754, -0.2009, 0.0640, ...
#>   SMD V2 values: {1, 3, 6, ...}  (linked B vs. B):  0.1432, 0.1637, -0.1392, ...
#>   SMD V3 values: {1, 3, 4, ...}  (linked B vs. B):  0.1734, 0.1404, -0.2650, ...
#>   SMD V4                         (linked B vs. B):  -0.0738
plot(diag, "scores")

plot(diag, "distributions", threshold = 0.75)


plot(diag, "convergence")










plot(diag, "FDP")

plot(diag, "discrepancy")



```
