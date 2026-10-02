# Stochastic Expectation-Maximisation for record linkage

Fits the FlexRL model with a Stochastic EM algorithm: each iteration
runs a Gibbs sampler (alternating between simulating the true PIV values
and the linkage matrix D) and then updates the model parameters
(`gamma`, `eta`, `alpha`, `phi`) from the post-burn-in Gibbs draws. See
the methodology paper https://10.1093/jrsssc/qlaf016 for details.

## Usage

``` r
StEM(
  data,
  StEM_iter = 30,
  StEM_burnin = 15,
  gibbs_iter = 20,
  gibbs_burnin = 10,
  music_on = FALSE,
  new_directory = NULL,
  save_info_iter = FALSE,
  gamma0 = NULL,
  phiA0 = NULL,
  phiB0 = NULL,
  n_post_samp = 1000,
  model_dynamics = survival_model("exponential")
)
```

## Arguments

- data:

  List, typically the output of [`prepare_data()`](prepare_data.md),
  with: `encodedA` (smaller source, PIVs encoded to natural numbers, `0`
  = missing), `encodedB` (larger source, encoded), `n_values`,
  `same_mistakes` (logical: same mistake parameter shared by A and B?),
  `PIVs_config` (named list, one entry per PIV, each a list with:
  `dynamics` (`"stable"`, `"flexible"`, or `"structured"`),
  `bound_mistakes` (length-2 numeric/NA, upper bound on the mistake
  probability in file 1 / file 2), `fix_mistakes` (length-2 numeric/NA,
  mistake probability fixed to this value in file 1 / file 2), and, only
  for `dynamics = "structured"`, `cond_hazard_cov` (a list with `cov1`
  and `cov2`, the names of covariates in file 1 / file 2 used to model
  the hazard of change).

- StEM_iter:

  Integer, total number of StEM iterations (including burn-in).

- StEM_burnin:

  Integer, number of StEM iterations discarded as burn-in.

- gibbs_iter:

  Integer, total number of Gibbs iterations per StEM step (including
  burn-in).

- gibbs_burnin:

  Integer, number of Gibbs iterations discarded as burn-in (`0` lets the
  algorithm auto-detect burn-in from the stabilisation of the linked
  count).

- music_on:

  Logical; if `TRUE`, opens a short tune in the browser when the
  algorithm finishes.

- new_directory:

  Path to an existing directory to save progress after each iteration,
  or `NULL` to disable.

- save_info_iter:

  Logical; save the environment at the end of each iteration (only used
  if `new_directory` is not `NULL`).

- gamma0:

  Optional starting values for `gamma`; at random if `NULL`.

- phiA0:

  Optional starting values for `phi`; at random if `NULL`.

- phiB0:

  Optional starting values for `phi`; at random if `NULL`.

- n_post_samp:

  Integer, number of posterior draws used to estimate the final linkage
  probabilities `Delta`. Default is set to 1000, we recommend not
  lowering it.

- model_dynamics:

  Object from [`survival_model()`](survival_model.md): the survival
  model for the change over time of the structured PIVs (default
  exponential).

## Value

A list with: `Delta` sparse-matrix summary (`i`, `j`, `x`) of posterior
linkage probabilities; a pair is a valid link candidate once `x > 0.5`
(one-to-one constraint), `gamma`, `eta`, `alpha`, `phi` the StEM chains
for each parameter.

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
                   
cond_hazard_params1 <- list( V1 = c(), V2 = c(), 
                             V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
gen_data <- simulate_data( PIVs_config, n_values, c(250, 300), 200, 
                           p_mistake, p_missing, cond_hazard_params1, TRUE )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> '2' is the larger source, saved as B; '1' saved as A.
fit1 <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
              gibbs_iter = 5, gibbs_burnin = 2, n_post_samp = 5,
              model_dynamics = survival_model("exponential") )
#> FlexRL
apply(fit1$alpha$V4[2:5,], 2, mean)
#>  intercept         Xe         Xf 
#> -2.1929366 -0.2846951 -0.2360663 
cond_hazard_params1$V4
#> [1] -0.3566749 -0.5108256 -0.6931472
head(fit1$Delta[fit1$Delta$x > 0.5, ])
#>    i  j   x
#> 1  1  1 1.0
#> 2  3  3 1.0
#> 4 15  5 1.0
#> 7 43  8 1.0
#> 8  9  9 1.0
#> 9 10 10 0.8

cond_hazard_params2 <- list( V1 = c(), V2 = c(), 
                             V3 = c(), V4 = log(c(0.3, 0.7, 0.6, 0.5)) )
gen_data <- simulate_data( PIVs_config, n_values, c(250, 300), 200, 
                           p_mistake, p_missing, cond_hazard_params2, TRUE,
                           model_dynamics = survival_model("weibull") )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> '2' is the larger source, saved as B; '1' saved as A.
fit2 <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
              gibbs_iter = 5, gibbs_burnin = 2, n_post_samp = 5,
              model_dynamics = survival_model("weibull"))
#> FlexRL
apply(fit2$alpha$V4[2:5,], 2, mean)
#>   log_shape   intercept          Xe          Xf 
#> -0.60357550 -0.88120552 -0.17256545  0.02085412 
cond_hazard_params2$V4
#> [1] -1.2039728 -0.3566749 -0.5108256 -0.6931472
head(fit2$Delta[fit2$Delta$x > 0.5, ])
#>    i j   x
#> 1  1 1 0.6
#> 2  2 2 0.6
#> 5  4 4 1.0
#> 6 50 5 1.0
#> 7  6 6 1.0
#> 8  7 7 1.0

cond_hazard_params3 <- list( V1 = c(), V2 = c(), 
                             V3 = c(), V4 = log(c(0.3, 0.7, 0.6, 0.5)) )
gen_data <- simulate_data( PIVs_config, n_values, c(250, 300), 200, 
                           p_mistake, p_missing, cond_hazard_params3, TRUE,
                           model_dynamics = survival_model("gompertz") )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> '2' is the larger source, saved as B; '1' saved as A.
fit3 <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
              gibbs_iter = 5, gibbs_burnin = 2, n_post_samp = 5,
              model_dynamics = survival_model("gompertz"))
#> FlexRL
apply(fit3$alpha$V4[2:5,], 2, mean)
#>            g    intercept           Xe           Xf 
#>  0.080275657 -1.994355256  0.007170501 -0.369704007 
cond_hazard_params3$V4
#> [1] -1.2039728 -0.3566749 -0.5108256 -0.6931472
head(fit3$Delta[fit3$Delta$x > 0.5, ])
#>     i j   x
#> 1   1 1 1.0
#> 2   3 3 0.8
#> 3   4 4 0.8
#> 4   6 6 1.0
#> 5   7 7 1.0
#> 6 129 8 1.0

cond_hazard_params4 <- list( V1 = c(), V2 = c(), 
                             V3 = c(), V4 = log(c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1)) )
gen_data <- simulate_data( PIVs_config, n_values, c(250, 300), 200, 
                           p_mistake, p_missing, cond_hazard_params4, TRUE,
                           model_dynamics = survival_model("piecewise", cuts = c(1,2,3)) )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> '2' is the larger source, saved as B; '1' saved as A.
fit4 <- StEM( data = prep_data, StEM_iter = 2, StEM_burnin = 1, 
              gibbs_iter = 2, gibbs_burnin = 1, n_post_samp = 2,
              model_dynamics = survival_model("piecewise", cuts = c(1,2,3)))
#> FlexRL
#> Warning: PIV 'V1', file A hit its mistake bound on 50% of StEM iterations.
#>   PIV 'V1', file B hit its mistake bound on 50% of StEM iterations.
apply(fit4$alpha$V4, 2, mean)
#>      log_h1      log_h2      log_h3      log_h4          Xe          Xf 
#>   0.1623534 -22.6890089 -18.7781821  -1.3672972   0.1358512   0.0770206 
cond_hazard_params4$V4
#> [1] -0.5108256 -0.6931472 -0.9162907 -1.2039728 -1.6094379 -2.3025851
head(fit4$Delta[fit4$Delta$x > 0.5, ])
#>     i  j x
#> 1  10  3 1
#> 2   4  4 1
#> 3 180  6 1
#> 4   8  8 1
#> 5   9  9 1
#> 6 214 10 1
```
