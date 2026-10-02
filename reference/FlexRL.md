# FlexRL: A Flexible Model for Record Linkage

Links records that refer to the same entities across two data sources
without a unique identifier, using partially identifying variables.
These are stable (not changing over time, probability of mistakes could
be bounded), flexible (dynamic but no information to model changes over
time) or structured (dynamic and information to model changes over time,
probability of mistakes could be fixed). The main function
[`StEM()`](stEM.md) fits a latent-variable model by stochastic
expectation-maximisation; it models registration errors (missing values
and mistakes) and changes over time. [`prepare_data()`](prepare_data.md)
prepares the data sources for record linkage,
[`RL_diagnostics()`](rl_diagnostics.md) gathers diagnostics (FDP
estimation and discrepancy metrics for inference on the linked data.

## Details

Methodological paper:
[doi:10.1093/jrsssc/qlaf016](https://doi.org/10.1093/jrsssc/qlaf016) .
False discovery proportion estimation:
[doi:10.1002/sim.70292](https://doi.org/10.1002/sim.70292) . Experiments
repository: <https://github.com/robachowyk/FlexRL-experiments>.

## See also

Useful links:

- <https://github.com/robachowyk/FlexRL>

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
fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 5,
             gibbs_iter = 10, gibbs_burnin = 5, n_post_samp = 10 )
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
#>           0.5000000           0.3012048           0.5800000 

# diagnostics for inference on the linked data
diag <- RL_diagnostics(fit, prep_data$encodedA, prep_data$encodedB,
                       names(PIVs_config),
                       list(V1 = FALSE, V2 = FALSE, V3 = FALSE, V4 = TRUE),
                       true_pairs = prep_data$true_pairs, FDP_estimation = TRUE, 
                       RL_method = "FlexRL", data = prep_data,
                       StEM_iter = 5, StEM_burnin = 2, 
                       gibbs_iter = 5, gibbs_burnin = 2, n_post_samp = 5,
                       maxIter4CV = 1, n_repeats = 2)
#> Iteration: 0, Accuracy: 46.21%
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
#> FlexRL
#> Iteration: 0, Accuracy: 51.9%
#> Iteration: 1, Accuracy: 46.43%
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> FlexRL
#> Warning: PIV 'V2', file A hit its mistake bound on 40% of StEM iterations.
#>   PIV 'V2', file B hit its mistake bound on 40% of StEM iterations.
#> FlexRL results (average over 2 iterations):
#>                            0.50  0.51  0.52  0.53  0.54  0.55  0.56  0.57  0.58
#> FDP model score estimator  0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14
#> FDP synth data estimator   0.43  0.43  0.43  0.43  0.43  0.43  0.43  0.43  0.43
#> Linked pairs (augm. RL)   70.00 70.00 70.00 70.00 70.00 70.00 70.00 70.00 70.00
#> Linked pairs (RL)         72.00 72.00 72.00 72.00 72.00 72.00 72.00 72.00 72.00
#>                            0.59  0.60  0.61  0.62  0.63  0.64  0.65  0.66  0.67
#> FDP model score estimator  0.14  0.06  0.06  0.06  0.06  0.06  0.06  0.06  0.06
#> FDP synth data estimator   0.43  0.27  0.27  0.27  0.27  0.27  0.27  0.27  0.27
#> Linked pairs (augm. RL)   70.00 55.00 55.00 55.00 55.00 55.00 55.00 55.00 55.00
#> Linked pairs (RL)         72.00 56.00 56.00 56.00 56.00 56.00 56.00 56.00 56.00
#>                            0.68  0.69  0.70  0.71  0.72  0.73  0.74  0.75  0.76
#> FDP model score estimator  0.06  0.06  0.06  0.06  0.06  0.06  0.06  0.06  0.06
#> FDP synth data estimator   0.27  0.27  0.27  0.27  0.27  0.27  0.27  0.27  0.27
#> Linked pairs (augm. RL)   55.00 55.00 55.00 55.00 55.00 55.00 55.00 55.00 55.00
#> Linked pairs (RL)         56.00 56.00 56.00 56.00 56.00 56.00 56.00 56.00 56.00
#>                            0.77  0.78  0.79  0.80  0.81  0.82  0.83  0.84  0.85
#> FDP model score estimator  0.06  0.06  0.06  0.00  0.00  0.00  0.00  0.00  0.00
#> FDP synth data estimator   0.27  0.27  0.27  0.29  0.29  0.29  0.29  0.29  0.29
#> Linked pairs (augm. RL)   55.00 55.00 55.00 38.00 38.00 38.00 38.00 38.00 38.00
#> Linked pairs (RL)         56.00 56.00 56.00 38.00 38.00 38.00 38.00 38.00 38.00
#>                            0.86  0.87  0.88  0.89  0.90  0.91  0.92  0.93  0.94
#> FDP model score estimator  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
#> FDP synth data estimator   0.29  0.29  0.29  0.29  0.29  0.29  0.29  0.29  0.29
#> Linked pairs (augm. RL)   38.00 38.00 38.00 38.00 38.00 38.00 38.00 38.00 38.00
#> Linked pairs (RL)         38.00 38.00 38.00 38.00 38.00 38.00 38.00 38.00 38.00
#>                            0.95  0.96  0.97  0.98  0.99
#> FDP model score estimator  0.00  0.00  0.00  0.00  0.00
#> FDP synth data estimator   0.29  0.29  0.29  0.29  0.29
#> Linked pairs (augm. RL)   38.00 38.00 38.00 38.00 38.00
#> Linked pairs (RL)         38.00 38.00 38.00 38.00 38.00
#> FlexRL
#> FlexRL
#> FlexRL results (average over 2 iterations):
#>                            0.50  0.51  0.52  0.53  0.54  0.55  0.56  0.57  0.58
#> FDP model score estimator  0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14  0.14
#> Linked pairs (RL)         76.00 76.00 76.00 76.00 76.00 76.00 76.00 76.00 76.00
#>                            0.59  0.60  0.61  0.62  0.63  0.64  0.65  0.66  0.67
#> FDP model score estimator  0.14  0.06  0.06  0.06  0.06  0.06  0.06  0.06  0.06
#> Linked pairs (RL)         76.00 59.00 59.00 59.00 59.00 59.00 59.00 59.00 59.00
#>                            0.68  0.69  0.70  0.71  0.72  0.73  0.74  0.75  0.76
#> FDP model score estimator  0.06  0.06  0.06  0.06  0.06  0.06  0.06  0.06  0.06
#> Linked pairs (RL)         59.00 59.00 59.00 59.00 59.00 59.00 59.00 59.00 59.00
#>                            0.77  0.78  0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86
#> FDP model score estimator  0.06  0.06  0.06    0    0    0    0    0    0    0
#> Linked pairs (RL)         59.00 59.00 59.00   41   41   41   41   41   41   41
#>                           0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96
#> FDP model score estimator    0    0    0    0    0    0    0    0    0    0
#> Linked pairs (RL)           41   41   41   41   41   41   41   41   41   41
#>                           0.97 0.98 0.99
#> FDP model score estimator    0    0    0
#> Linked pairs (RL)           41   41   41
diag # print(diag_flexrl)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.50):  83
#> 
#>   FDP           (true pairs):    0.301
#>   Sensitivity   (true pairs):    0.580
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.139
#>                                                threshold 0.50: FDP ~ 0.139
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.138
#>                                                threshold 0.50: FDP ~ 0.138
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.432
#>                                                threshold 0.50: FDP ~ 0.432
#>                                            max threshold 0.99: FDP ~ 0.286
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0078
#>   Multivariate MMD  (linked B vs. B):  0.0104
#> 
#>   IoU support V1  (linked A vs. A):  1.0000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  1.0000
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  0.9167
#> 
#>   Agreement V1  (linked A vs. linked B):  0.9880
#>   Agreement V2  (linked A vs. linked B):  0.9880
#>   Agreement V3  (linked A vs. linked B):  0.9398
#>   Agreement V4  (linked A vs. linked B):  0.6988
#>   Agreement V1  (true pairs):  0.9500
#>   Agreement V2  (true pairs):  0.9500
#>   Agreement V3  (true pairs):  0.9100
#>   Agreement V4  (true pairs):  0.7100
#> 
#>   SMD V1 values: {2, 3, 4, ...}  (linked A vs. A):  0.0956, 0.0758, -0.1394, ...
#>   SMD V2 values: {1, 2, 5, ...}  (linked A vs. A):  0.0641, 0.0517, -0.0519, ...
#>   SMD V3 values: {3, 5, 7, ...}  (linked A vs. A):  0.0824, 0.0891, -0.1134, ...
#>   SMD V4                         (linked A vs. A):  -0.0087
#>   SMD V1 values: {0, 1, 2, ...}  (linked B vs. B):  -0.1003, -0.1095, 0.0957, ...
#>   SMD V2 values: {0, 1, 2, ...}  (linked B vs. B):  -0.0707, 0.0469, 0.0457, ...
#>   SMD V3 values: {2, 3, 7, ...}  (linked B vs. B):  -0.1037, 0.1353, -0.0918, ...
#>   SMD V4                         (linked B vs. B):  -0.0274
print(diag, threshold = 0.75)
#> <RL_diagnostics>
#> 
#>   Linked pairs  (score > 0.75):  61
#> 
#>   FDP           (true pairs):    0.213
#>   Sensitivity   (true pairs):    0.480
#> 
#>   FDP score estimate (RL task):            min threshold 0.50: FDP ~ 0.139
#>                                                threshold 0.75: FDP ~ 0.061
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP score estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.138
#>                                                threshold 0.75: FDP ~ 0.063
#>                                            max threshold 0.99: FDP ~ 0.000
#>   FDP synth estimate (augmented RL task):  min threshold 0.50: FDP ~ 0.432
#>                                                threshold 0.75: FDP ~ 0.271
#>                                            max threshold 0.99: FDP ~ 0.286
#>   The score-based estimator is valid when the linkage model is well calibrated to the data.
#>   The synthetic-data estimator is valid when the augmented task is equivalent to the
#>   original one, which requires links and non-links to have similar distributions; similar
#>   score-based FDP estimates on the original and augmented tasks support this prerequisite.
#> 
#>   Multivariate MMD  (linked A vs. A):  0.0109
#>   Multivariate MMD  (linked B vs. B):  0.0110
#> 
#>   IoU support V1  (linked A vs. A):  0.8000
#>   IoU support V2  (linked A vs. A):  1.0000
#>   IoU support V3  (linked A vs. A):  1.0000
#>   IoU support V4  (linked A vs. A):  1.0000
#>   IoU support V1  (linked B vs. B):  0.8000
#>   IoU support V2  (linked B vs. B):  0.8333
#>   IoU support V3  (linked B vs. B):  0.8571
#>   IoU support V4  (linked B vs. B):  1.0000
#> 
#>   Agreement V1  (linked A vs. linked B):  1.0000
#>   Agreement V2  (linked A vs. linked B):  1.0000
#>   Agreement V3  (linked A vs. linked B):  0.9672
#>   Agreement V4  (linked A vs. linked B):  0.8033
#>   Agreement V1  (true pairs):  0.9500
#>   Agreement V2  (true pairs):  0.9500
#>   Agreement V3  (true pairs):  0.9100
#>   Agreement V4  (true pairs):  0.7100
#> 
#>   SMD V1 values: {0, 3, 4, ...}  (linked A vs. A):  -0.0816, 0.1200, -0.1207, ...
#>   SMD V2 values: {2, 3, 5, ...}  (linked A vs. A):  0.0734, 0.0626, -0.1097, ...
#>   SMD V3 values: {0, 3, 7, ...}  (linked A vs. A):  0.1191, 0.1158, -0.1121, ...
#>   SMD V4                         (linked A vs. A):  -0.0786
#>   SMD V1 values: {0, 2, 5, ...}  (linked B vs. B):  -0.1003, 0.0434, -0.0357, ...
#>   SMD V2 values: {0, 2, 4, ...}  (linked B vs. B):  -0.0707, 0.0672, -0.0610, ...
#>   SMD V3 values: {2, 3, 7, ...}  (linked B vs. B):  -0.1423, 0.1465, -0.1074, ...
#>   SMD V4                         (linked B vs. B):  -0.0279
plot(diag, "scores")

plot(diag, "distributions", threshold = 0.75)


plot(diag, "convergence")










plot(diag, "FDP")

plot(diag, "discrepancy")



```
