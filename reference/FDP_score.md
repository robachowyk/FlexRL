# Model specific, linkage score based, false discovery proportion at a given threshold

Model specific, linkage score based, false discovery proportion at a
given threshold

## Usage

``` r
FDP_score(LinkScore, threshold)
```

## Arguments

- LinkScore:

  Numeric vector, linkage scores of the candidate pairs.

- threshold:

  Numeric, score threshold above which a pair is declared linked.

## Value

A list with `FDP_score` Numeric, `1 - mean(score | score > threshold)`,
the estimated FDP at `threshold` and `n_linked` Integer, number of
linked records.

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
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> '2' is the larger source, saved as B; '1' saved as A.
fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 5,
             gibbs_iter = 10, gibbs_burnin = 5, n_post_samp = 10 )
#> FlexRL
#> Running StEM algorithm ■■■■■■■                           20% | iter 2/10 [1.1s]
#> Running StEM algorithm ■■■■■■■■■■                        30% | iter 3/10 [1.3s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | iter 10/10 [2.4s]
#> 
FDP_score(fit$Delta$x, 0.5)
#> $FDP_score
#> [1] 0.09032258
#> 
#> $n_linked
#> [1] 186
#> 
FDP_score(fit$Delta$x, 0.75)
#> $FDP_score
#> [1] 0.04556962
#> 
#> $n_linked
#> [1] 158
#> 
```
