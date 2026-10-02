# Simulate the linkage matrix D

Given the current draw of true PIV values, samples the linkage matrix D
from its conditional distribution, and returns the updated
log-likelihood.

## Usage

``` r
simulateD(
  data,
  linksR,
  sumRowD,
  sumColD,
  truepivsA,
  truepivsB,
  gamma,
  eta,
  alpha,
  phi,
  model_dynamics = survival_model("exponential")
)
```

## Arguments

- data:

  List with `encodedA`, `encodedB` (the two encoded data sources,
  missing values as `0`), `n_values`, `PIVs_config` and `same_mistakes`;
  see [`StEM()`](stEM.md).

- linksR:

  2-column matrix (1-indexed) of the currently linked (A, B) indices.

- sumRowD:

  Logical vector, one entry per record in A: does it form a link?

- sumColD:

  Logical vector, one entry per record in B: does it form a link?

- truepivsA:

  Matrix of true PIV values, as returned by
  [`simulateH()`](simulateH.md).

- truepivsB:

  Matrix of true PIV values, as returned by
  [`simulateH()`](simulateH.md).

- gamma:

  Numeric, proportion of linked records as a fraction of the smaller
  file.

- eta:

  List (one per PIV) of the distribution of true values.

- alpha:

  List (one per PIV) of survival-model parameters (see
  [`survival_model()`](survival_model.md)).

- phi:

  List (one per PIV) of registration-error parameters.

- model_dynamics:

  Object from [`survival_model()`](survival_model.md) giving the
  survival function of the structured PIVs.

## Value

A list (`Dsample`, Rcpp `sampleD()` output) with: `links` updated set of
links, `sumRowD` updated sumRowD, `sumColD` updated sumColD, `loglik`
updated value of the complete log likelihood, `nlinkrec` updated number
of linked records

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
Dsample = simulateD(data=data_StEM, linksR=linksR, sumRowD=sumRowD,
                   sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB,
                   gamma=gamma, eta=eta, alpha=alpha, phi=phi)
linksCpp = Dsample$links
linksR = linksCpp + 1
```
