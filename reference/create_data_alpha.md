# Book-keeping data frame for parameters of PIVs dynamics

Internal helper used by [`StEM()`](stEM.md) to accumulate, across Gibbs
iterations, the covariates, true-value agreement indicator, and time
gaps needed to re-estimate the survival (hazard) parameters of a dynamic
structured PIV.

## Usage

``` r
create_data_alpha(n_coef_unstable, stable)
```

## Arguments

- n_coef_unstable:

  Integer, number of hazard coefficients for this PIV (1 for the
  baseline hazard, plus one per covariate from file A and file B).

- stable:

  Logical, whether this PIV is stable (`dynamics != "structured"`).

## Value

An empty data frame with `n_coef_unstable + 2` columns (the
covariates/intercept, plus `Hequal` and `times`) if `stable` is `FALSE`;
`NULL` if the PIV is stable (nothing to accumulate).

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
PIVs_stable <- sapply(PIVs_config, function(x) x$dynamics != "structured")
n_coef_unstable = c(0,0,0,3)
Valpha <- mapply(create_data_alpha, n_coef_unstable = n_coef_unstable,
                 stable = PIVs_stable, SIMPLIFY = FALSE)
```
