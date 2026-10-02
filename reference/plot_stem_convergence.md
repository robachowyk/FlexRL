# Monte Carlo convergence (trace) plots for a fitted StEM model

Trace-plots the raw StEM chains of `gamma`, `eta`, `alpha` and `phi`
across iterations, to visually judge whether the chains have stabilised
and pick an adequate `StEM_burnin` for [`StEM()`](stEM.md). One plot per
parameter: `gamma` (proportion linked), `eta` (PIVs distribution),
`alpha` (hazard coefficients, for dynamic PIVs only), and `phi`
(agreement/missing rates).

## Usage

``` r
plot_StEM_convergence(fit)
```

## Arguments

- fit:

  List as returned by [`StEM()`](stEM.md), containing the raw chains
  `gamma`, `eta`, `alpha`, `phi` (StEM_iter rows each).

## Value

`NULL`, invisibly; called for its plotting side effect.

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
fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 3,
             gibbs_iter = 10, gibbs_burnin = 3, n_post_samp = 10 )
#> FlexRL
plot_StEM_convergence(fit)









```
