pkgname <- "FlexRL"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "FlexRL-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('FlexRL')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("FDP_score")
### * FDP_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FDP_score
### Title: Model specific, linkage score based, false discovery proportion
###   at a given threshold
### Aliases: FDP_score

### ** Examples

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
fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 5,
             gibbs_iter = 10, gibbs_burnin = 5, n_post_sample = 10 )
FDP_score(fit$Delta$x, 0.5)
FDP_score(fit$Delta$x, 0.75)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FDP_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("FDP_synth")
### * FDP_synth

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FDP_synth
### Title: Model agnostic, based on synthetic data, false discovery
###   proportion at a given threshold
### Aliases: FDP_synth

### ** Examples

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
PIVs <- names(PIVs_config)
n_synth <-as.integer(0.10 * nrow(prep_data$encodedB))
new_data <- synthesise("arf", prep_data$encodedA[, c(PIVs,"local_id","source",
                      "date",prep_data$PIVs_config$V4$cond_hazard_cov$covA)],
                      prep_data$encodedB[, c(PIVs,"local_id","source","date",
                      prep_data$PIVs_config$V4$cond_hazard_cov$covB)], PIVs, 
                      n_synth, TRUE)
new_data$dataA[PIVs][is.na(new_data$dataA[PIVs])] <- 0
new_data$dataB[PIVs][is.na(new_data$dataB[PIVs])] <- 0
# cannot model dynamics for synthetic data, set dates to 0
new_data$dataA$date[is.na(new_data$dataA$date)] <- 0
new_data$dataB$date[is.na(new_data$dataB$date)] <- 0
arguments <- list(data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
                  gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5)
fit_flexrl <- link_with_FlexRL(new_data$dataA, new_data$dataB, arguments)  
FDP_synth(fit_flexrl$idxA, fit_flexrl$idxB, fit_flexrl$LinkScore, 0.5,
          nrow(prep_data$encodedA), nrow(prep_data$encodedB), 
          n_synth)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FDP_synth", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("FlexRL")
### * FlexRL

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FlexRL
### Title: FlexRL: A Flexible Model for Record Linkage
### Aliases: FlexRL-package FlexRL

### ** Examples

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
fit <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2,
             gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5 )

# linked pairs and performance against the true pairs
linked <- fit$Delta[fit$Delta$x > 0.5, ]
linked_pairs <- paste(linked$i, linked$j, sep = "_")
true_pairs   <- paste(prep_data$true_pairs[[1]], prep_data$true_pairs[[2]], sep = "_")
tp <- length(intersect(linked_pairs, true_pairs))
fp <- length(setdiff(linked_pairs, true_pairs))
fn <- length(setdiff(true_pairs, linked_pairs))
c(LinkageDecisionRule = 0.5, FDP = fp / (tp + fp), Sensitivity = tp / (tp + fn))

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
diag # print(diag)
print(diag, threshold = 0.75)
plot(diag, "scores")
plot(diag, "distributions", threshold = 0.75)
plot(diag, "convergence")
plot(diag, "FDP")
plot(diag, "discrepancy")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FlexRL", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("RL_agreement")
### * RL_agreement

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: RL_agreement
### Title: Agreement rate of pairs on shared variables
### Aliases: RL_agreement

### ** Examples

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
fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 5, 
             gibbs_iter = 10, gibbs_burnin = 5, n_post_sample = 10 )                         
linked_pairs <- fit$Delta[fit$Delta$x > 0.5, ]
RL_agreement( prep_data$encodedA, prep_data$encodedB,
              names(PIVs_config), linked_pairs )
RL_agreement( prep_data$encodedA, prep_data$encodedB,
              names(PIVs_config), prep_data$true_pairs )



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("RL_agreement", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("RL_diagnostics")
### * RL_diagnostics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: RL_diagnostics
### Title: Post-linkage diagnostics
### Aliases: RL_diagnostics

### ** Examples

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
diag_brl # print(diag_brl)
print(diag_brl, threshold = 0.75)
plot(diag_brl, "scores")
plot(diag_brl, "distributions", threshold = 0.75)
plot(diag_brl, "FDP")
plot(diag_brl, "discrepancy")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("RL_diagnostics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("StEM")
### * StEM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: StEM
### Title: Stochastic Expectation-Maximisation for record linkage
### Aliases: StEM

### ** Examples

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
                   
cond_hazard_params1 <- list( V1 = c(), V2 = c(), 
                             V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
                           p_mistake, p_missing, cond_hazard_params1, TRUE )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
fit1 <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
              gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5,
              model_dynamics = survival_model("exponential") )
apply(fit1$alpha$V4[2:5,], 2, mean)
cond_hazard_params1$V4
head(fit1$Delta[fit1$Delta$x > 0.5, ])

cond_hazard_params2 <- list( V1 = c(), V2 = c(), 
                             V3 = c(), V4 = log(c(0.3, 0.7, 0.6, 0.5)) )
gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
                           p_mistake, p_missing, cond_hazard_params2, TRUE,
                           model_dynamics = survival_model("weibull") )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
fit2 <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
              gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5,
              model_dynamics = survival_model("weibull"))
apply(fit2$alpha$V4[2:5,], 2, mean)
cond_hazard_params2$V4
head(fit2$Delta[fit2$Delta$x > 0.5, ])

cond_hazard_params3 <- list( V1 = c(), V2 = c(), 
                             V3 = c(), V4 = log(c(0.3, 0.7, 0.6, 0.5)) )
gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
                           p_mistake, p_missing, cond_hazard_params3, TRUE,
                           model_dynamics = survival_model("gompertz") )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
fit3 <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
              gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5,
              model_dynamics = survival_model("gompertz"))
apply(fit3$alpha$V4[2:5,], 2, mean)
cond_hazard_params3$V4
head(fit3$Delta[fit3$Delta$x > 0.5, ])

cond_hazard_params4 <- list( V1 = c(), V2 = c(), 
                             V3 = c(), V4 = log(c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1)) )
gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
                           p_mistake, p_missing, cond_hazard_params4, TRUE,
                           model_dynamics = survival_model("piecewise", cuts = c(1,2,3)) )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
fit4 <- StEM( data = prep_data, StEM_iter = 2, StEM_burnin = 1, 
              gibbs_iter = 2, gibbs_burnin = 1, n_post_sample = 2,
              model_dynamics = survival_model("piecewise", cuts = c(1,2,3)))
apply(fit4$alpha$V4, 2, mean)
cond_hazard_params4$V4
head(fit4$Delta[fit4$Delta$x > 0.5, ])



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("StEM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_RL_FDP_score")
### * compute_RL_FDP_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_RL_FDP_score
### Title: Estimate the false discovery proportion of a record-linkage
###   method via model specific linkage scores
### Aliases: compute_RL_FDP_score

### ** Examples

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
PIVs <- names(PIVs_config)                          
compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
                      "BRL", flds = PIVs, types = rep("bi",length(PIVs)) )
compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
                      "FlexRL", data = prep_data, StEM_iter = 5, 
                      StEM_burnin = 2, gibbs_iter = 5, gibbs_burnin = 2,
                      n_post_sample = 5 )



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_RL_FDP_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_augmRL_FDP_synth")
### * compute_augmRL_FDP_synth

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_augmRL_FDP_synth
### Title: Estimate the false discovery proportion of a record-linkage
###   method via synthetic augmentation
### Aliases: compute_augmRL_FDP_synth

### ** Examples

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
PIVs <- names(PIVs_config)                         
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "BRL", flds = PIVs, 
                          types = rep("bi",length(PIVs)) )
compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
                          NULL, TRUE, 1, 2, "FlexRL", data = prep_data, 
                          StEM_iter = 5, StEM_burnin = 2, 
                          gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5 )



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_augmRL_FDP_synth", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("create_data_alpha")
### * create_data_alpha

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create_data_alpha
### Title: Book-keeping data frame for parameters of PIVs dynamics
### Aliases: create_data_alpha

### ** Examples

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
PIVs_stable <- sapply(PIVs_config, function(x) x$dynamics != "structured")
n_coef_unstable = c(0,0,0,3)
Valpha <- mapply(create_data_alpha, n_coef_unstable = n_coef_unstable,
                 stable = PIVs_stable, SIMPLIFY = FALSE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create_data_alpha", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_lik")
### * log_lik

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_lik
### Title: Log-likelihood of the linkage matrix
### Aliases: log_lik

### ** Examples

LLL <- Matrix::Matrix(0, nrow = 13, ncol = 15, sparse = TRUE)
LLA <- stats::runif(13, 0, 2)
LLB <- stats::runif(15, 0, 2)
links <- as.matrix(data.frame(idxA = c(5, 9, 11, 12, 13),
                              idxB = c(5, 9, 11, 13, 15)))
LLL[links] <- 0.67
sumRowD <- (seq_len(13) %in% links[, 1])
sumColD <- (seq_len(15) %in% links[, 2])
gamma <- 0.5
log_lik(LLL, LLA, LLB, links, sumRowD, sumColD, gamma)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_lik", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("log_possible_config")
### * log_possible_config

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: log_possible_config
### Title: Log number of possible linkage configurations
### Aliases: log_possible_config

### ** Examples

log_possible_config(n_records_B = 15, sumD = 5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("log_possible_config", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mmd")
### * mmd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mmd
### Title: Maximum Mean Discrepancy between two sets of variables
### Aliases: mmd

### ** Examples

x <- matrix(rnorm(50), ncol = 5)
y <- matrix(rnorm(30) + 0.5, ncol = 5)
mmd(x, y)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mmd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("naive_linkage")
### * naive_linkage

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: naive_linkage
### Title: Naive (deterministic exact-match) record linkage
### Aliases: naive_linkage

### ** Examples

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
naive_linkage( names(prep_data$PIVs_config), 
               prep_data$encodedA, prep_data$encodedB )



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("naive_linkage", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_StEM_convergence")
### * plot_StEM_convergence

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_StEM_convergence
### Title: Monte Carlo convergence (trace) plots for a fitted StEM model
### Aliases: plot_StEM_convergence

### ** Examples

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
fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 3,
             gibbs_iter = 10, gibbs_burnin = 3, n_post_sample = 10 )
plot_StEM_convergence(fit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_StEM_convergence", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_distributions")
### * plot_distributions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_distributions
### Title: Compare distributions of shared variables across data sets
###   (histograms/barplots)
### Aliases: plot_distributions

### ** Examples

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
fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 5,
             gibbs_iter = 10, gibbs_burnin = 5, n_post_sample = 10 )
threshold_strict <- stats::quantile(fit$Delta$x, 0.75)
data_list = list(data_baseline = prep_data$encodedA, 
data_select = prep_data$encodedA[fit$Delta[fit$Delta$x > threshold_strict, "i"],])
common_vars = names(PIVs_config)
plot_distributions(data_list, common_vars)
common_vars = c("Xe", "Xf")
plot_distributions(data_list, common_vars)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_distributions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_linkage_scores")
### * plot_linkage_scores

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_linkage_scores
### Title: Plot the distribution of linkage scores
### Aliases: plot_linkage_scores

### ** Examples

fit_Delta_x <- c(0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.5,0.6,0.7,0.7,0.7,0.7,0.7,0.7,0.8,0.8)
plot_linkage_scores(1000, fit_Delta_x)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_linkage_scores", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("prepare_data")
### * prepare_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: prepare_data
### Title: Prepare two data sources for 'StEM()'
### Aliases: prepare_data

### ** Examples

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
str(prep_data, max.level = 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("prepare_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("prop_level")
### * prop_level

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: prop_level
### Title: Proportion of a data set with a given variable at a given level
### Aliases: prop_level

### ** Examples

df <- data.frame(colour = sample(c("orange", "purple"), 100, replace = TRUE))
prop_level(df, "colour", "purple")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("prop_level", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulateD")
### * simulateD

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulateD
### Title: Simulate the linkage matrix D
### Aliases: simulateD

### ** Examples

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
data_StEM <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
PIVs_stable <- sapply(data_StEM$PIVs_config, function(x)
                        x$dynamics != "structured")
linksR = base::matrix(0,0,2)
linksCpp = linksR
sumRowD = rep(0, nrow(data_StEM$encodedA))
sumColD = rep(0, nrow(data_StEM$encodedB))
nlinkrec = 0
survivalpSameH = base::matrix(1, nrow(linksR), length(data_StEM$n_values))
gamma = 0.5
eta = lapply(data_StEM$n_values, function(x) rep(1/x,x))
phi = lapply(data_StEM$n_values, function(x)  c(0.9,0.9,0.1,0.1))
n_coef_unstable = lapply( seq_along(PIVs_stable), function(idx)
 if(PIVs_stable[idx]){ 0 }else{
   ncol(data_StEM$encodedA[, data_StEM$PIVs_config[[idx]]$cond_hazard_cov$covA,
                             drop=FALSE]) +
   ncol(data_StEM$encodedB[, data_StEM$PIVs_config[[idx]]$cond_hazard_cov$covB,
                             drop=FALSE]) + 1 } )
alpha = lapply( seq_along(PIVs_stable),
                function(idx) if(PIVs_stable[idx]){ c(-Inf) }
                else{ rep(log(0.05), n_coef_unstable[[idx]]) }
              )
newTruePivs = simulateH(data=data_StEM, links=linksCpp,
                        survivalpSameH=survivalpSameH,
                        sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
truepivsA = newTruePivs$truepivsA
truepivsB = newTruePivs$truepivsB
Dsample = simulateD(data=data_StEM, linksR=linksR,
                   truepivsA=truepivsA, truepivsB=truepivsB,
                   gamma=gamma, eta=eta, alpha=alpha, phi=phi)
linksCpp = Dsample$links
linksR = linksCpp + 1



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulateD", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulateH")
### * simulateH

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulateH
### Title: Simulate the true PIV values underlying the registered records
### Aliases: simulateH

### ** Examples

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
data_StEM <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
PIVs_stable <- sapply(data_StEM$PIVs_config, function(x)
                        x$dynamics != "structured")
linksR = base::matrix(0,0,2)
linksCpp = linksR
sumRowD = rep(0, nrow(data_StEM$encodedA))
sumColD = rep(0, nrow(data_StEM$encodedB))
nlinkrec = 0
survivalpSameH = base::matrix(1, nrow(linksR), length(data_StEM$n_values))
gamma = 0.5
eta = lapply(data_StEM$n_values, function(x) rep(1/x,x))
phi = lapply(data_StEM$n_values, function(x)  c(0.9,0.9,0.1,0.1))
n_coef_unstable = lapply( seq_along(PIVs_stable), function(idx)
 if(PIVs_stable[idx]){ 0 }else{
   ncol(data_StEM$encodedA[, data_StEM$PIVs_config[[idx]]$cond_hazard_cov$covA,
                             drop=FALSE]) +
   ncol(data_StEM$encodedB[, data_StEM$PIVs_config[[idx]]$cond_hazard_cov$covB,
                             drop=FALSE]) + 1 } )
alpha = lapply( seq_along(PIVs_stable),
                function(idx) if(PIVs_stable[idx]){ c(-Inf) }
                else{ rep(log(0.05), n_coef_unstable[[idx]]) }
              )
newTruePivs = simulateH(data=data_StEM, links=linksCpp,
                        survivalpSameH=survivalpSameH,
                        sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
truepivsA = newTruePivs$truepivsA
truepivsB = newTruePivs$truepivsB



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulateH", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_data")
### * simulate_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_data
### Title: Simulate two linked data sources for record linkage benchmarking
### Aliases: simulate_data

### ** Examples

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
str(gen_data, max.level = 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("smd")
### * smd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: smd
### Title: Standardised mean difference between a selected set and a
###   baseline
### Aliases: smd

### ** Examples

base <- data.frame(age = rnorm(200, 40, 10), sex = sample(c("M", "F"),
                    200, replace = TRUE))
select <- data.frame(age = rnorm(80, 43, 10), sex = sample(c("M", "F"), 80,
                    replace = TRUE, prob = c(0.6, 0.4)))
smd(select, base, "age")
smd(select, base, "sex", continuous = FALSE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("smd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("support_iou")
### * support_iou

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: support_iou
### Title: Intersection-over-union of two histogram supports
### Aliases: support_iou

### ** Examples

h1 <- graphics::hist(rnorm(200), plot = FALSE)
h2 <- graphics::hist(rnorm(200) + 1, plot = FALSE)
support_iou(h1, h2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("support_iou", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("survival_model")
### * survival_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: survival_model
### Title: Survival models for the dynamics of a PIV
### Aliases: survival_model

### ** Examples

X <- cbind(intercept = rep(1, 5))
times <- c(0.001, 0.2, 1.3, 1.5, 2)
Hequal <- c(TRUE, TRUE, TRUE, FALSE, FALSE)

expo <- survival_model("exponential")
expo$S(X, alpha = log(0.3), times)
stats::nlminb(expo$init(ncol(X)), expo$negloglik, X = X, times = times,
              Hequal = Hequal)$par

weib <- survival_model("weibull")
weib$S(X, alpha = c(0.5, -1), times)
stats::nlminb(weib$init(ncol(X)), weib$negloglik, X = X, times = times,
              Hequal = Hequal)$par

gomp <- survival_model("gompertz")
gomp$S(X, alpha = c(0.5, -1), times)
stats::nlminb(gomp$init(ncol(X)), gomp$negloglik, X = X, times = times,
              Hequal = Hequal)$par
              
pw <- survival_model("piecewise", cuts = 1)
pw$S(X, alpha = c(-2, -1), times)
stats::nlminb(pw$init(ncol(X)), pw$negloglik, X = X, times = times,
              Hequal = Hequal)$par

# a custom model: log-logistic with scale exp(X alpha[-1]) 
# and shape exp(alpha[1])
loglogistic <- survival_model("custom",
  S = function(X, alpha, times) {
    lambda <- exp(as.matrix(X) %*% alpha[-1])
    1 / ( 1 + (lambda * times)^exp(alpha[1]) )
  },
  n_par = function(n_cov) n_cov + 1,
  init = function(n_cov) c(stats::rnorm(1, 0, 0.1), stats::runif(n_cov, log(0.01), log(1))))
loglogistic$S(X, alpha = c(0, -2), times)
stats::nlminb(loglogistic$init(ncol(X)), loglogistic$negloglik, X = X, times = times,
              Hequal = Hequal)$par



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("survival_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
