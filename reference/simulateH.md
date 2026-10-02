# Simulate the true PIV values underlying the registered records

Draw the latent true values of each PIV given the currently registered
(possibly mistaken or missing) values, the current linkage status, and
the current parameters.

## Usage

``` r
simulateH(data, links, survivalpSameH, sumRowD, sumColD, eta, phi)
```

## Arguments

- data:

  List with `encodedA`, `encodedB` (the two encoded data sources,
  missing values as `0`), `n_values`, `PIVs_config` and `same_mistakes`;
  see [`StEM()`](stEM.md).

- links:

  2-column matrix of (A, B) indices for the currently linked records.

- survivalpSameH:

  Matrix (n links x n PIVs); `1` for stable PIVs, and the survival
  probability that the true value is unchanged for unstable PIVs.

- sumRowD:

  Logical vector, one entry per record in A: does it form a link?

- sumColD:

  Logical vector, one entry per record in B: does it form a link?

- eta:

  List (one per PIV) of the distribution of true values.

- phi:

  List (one per PIV) of length-4 vectors: agreement prob. in A,
  agreement prob. in B, missing prob. in A, missing prob. in B.

## Value

List with `truepivsA` and `truepivsB`, matrices (same shape as the input
data) of simulated true PIV values.

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
gen_data <- simulate_data( PIVs_config, n_values, c(250, 300), 200, 
                           p_mistake, p_missing, cond_hazard_params, TRUE )
data_StEM <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs_stable <- sapply(data_StEM$PIVs_config, function(x)
                        x$dynamics != "structured")
FlexRL:::initDeltaMap()
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
```
