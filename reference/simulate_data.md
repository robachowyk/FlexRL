# Simulate two linked data sources for record linkage benchmarking

Creates two synthetic data sources of given sizes sharing a given number
of common entities ("links"), each described by a set of Partially
Identifying Variables (PIVs). For every PIV, choose the number of
possible values, the proportion of mistakes and of missing values, and
whether it is stable over time, flexible (may change but change is not
modelled), or structured (expected to change over time, with a survival
model for the hazard of change). For structured PIVs,
`enforce_estimability` forces half of the linked pairs to have a
near-zero time gap, which helps separate "mistake" from "change over
time" when fitting the model.

## Usage

``` r
simulate_data(
  PIVs_config,
  n_values,
  n_records,
  n_links,
  p_mistake,
  p_missing,
  cond_hazard_params,
  enforce_estimability,
  model_dynamics = survival_model("exponential")
)
```

## Arguments

- PIVs_config:

  Named list, one entry per PIV, each a list with: `dynamics`
  (`"stable"`, `"flexible"`, or `"structured"`), `bound_mistakes`
  (length-2 numeric/NA, upper bound on the mistake probability in file 1
  / file 2), `fix_mistakes` (length-2 numeric/NA, mistake probability
  fixed to this value in file 1 / file 2), and, only for
  `dynamics = "structured"`, `cond_hazard_cov` (a list with `cov1` and
  `cov2`, the names of covariates in file 1 / file 2 used to model the
  hazard of change).

- n_values:

  Integer vector, number of unique values per PIV (same order as
  `PIVs_config`).

- n_records:

  Integer vector of length 2, number of records to generate in file 1
  and file 2 (file 2 must be the larger of the two).

- n_links:

  Integer, number of records shared between the two files.

- p_mistake:

  Named list (one entry per PIV) of length-2 numeric vectors, proportion
  of mistakes to introduce in file 1 / file 2.

- p_missing:

  Named list (one entry per PIV) of length-2 numeric vectors, proportion
  of missing values to introduce in file 1 / file 2.

- cond_hazard_params:

  Named list (one entry per PIV) of numeric vectors for the survival
  model generating the changes (see
  [`survival_model()`](survival_model.md)); only used for
  `dynamics = "structured"` PIVs.

- enforce_estimability:

  Logical; if `TRUE`, half of the linked pairs are given a near-zero
  time gap to help estimate the instability parameters.

- model_dynamics:

  Object from [`survival_model()`](survival_model.md) used to generate
  the changes of the structured PIVs (default exponential).

## Value

A list with: `data1`, `data2` the two simulated (encoded) data frames,
`n_values` number of unique values per PIV, `time_difference` time gap
between linked records (`NA` if no PIV is structured), `proba_same_H`
matrix (n links x n PIVs) of probabilities that the true values
coincide, `true_pairs` data frame with the true `1`/`2` indices of the
linked records

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
str(gen_data, max.level = 1)
#> List of 6
#>  $ data1          :'data.frame': 250 obs. of  11 variables:
#>  $ data2          :'data.frame': 300 obs. of  9 variables:
#>  $ n_values       : int [1:4] 5 6 7 12
#>  $ time_difference: num [1:200] 0.00673 0.00642 0.00671 0.00503 0.00642 ...
#>  $ proba_same_H   : num [1:200, 1:4] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ true_pairs     :'data.frame': 200 obs. of  2 variables:
```
