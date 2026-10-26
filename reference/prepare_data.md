# Prepare two data sources for [`StEM()`](https://robachowyk.github.io/FlexRL/reference/StEM.md)

Wraps the data-preparation steps needed before calling
[`StEM()`](https://robachowyk.github.io/FlexRL/reference/StEM.md): tags
each source with a `source` column, labels the larger source as `B`,
adjusting `cond_hazard_cov`/`bound_mistakes`/`fix_mistakes` accordingly
if A and B are swapped), drops records whose PIV values fall outside the
support shared by both files, warns if two PIVs are strongly associated
(Cramer's V \> 0.3), which may degrade record linkage performance,
encodes every PIV to natural numbers using levels pooled across both
sources, and encodes missing values to `0`.

## Usage

``` r
prepare_data(
  data1,
  data2,
  label1,
  label2,
  PIVs_config,
  same_mistakes = TRUE,
  uniq_id = NULL,
  restrict_support_intersection = TRUE
)
```

## Arguments

- data1:

  Data frame, the raw data source (whichever has more rows becomes `B`).

- data2:

  Data frame, the raw data source (whichever has more rows becomes `B`).

- label1:

  Character, label recorded in the `source` column for `data1`.

- label2:

  Character, label recorded in the `source` column for `data2`.

- PIVs_config:

  Named list describing each PIV — see
  [`simulate_data()`](https://robachowyk.github.io/FlexRL/reference/simulate_data.md).

- same_mistakes:

  Logical, will A and B share one mistake-probability parameter per PIV.

- uniq_id:

  Optional column name (present in both files) with the true entity
  identifier, used to build `true_pairs` for evaluation; `NULL` if
  unavailable.

- restrict_support_intersection:

  Logical; if `TRUE` (default), records with an out-of-common-support
  PIV value are dropped (and a warning issued); if `FALSE`, only the
  warning is issued.

## Value

A list ready to use as the `data` argument of
[`StEM()`](https://robachowyk.github.io/FlexRL/reference/StEM.md):
`encodedA`, `encodedB`, `n_values`, `PIVs_config`, `same_mistakes`, and
`true_pairs` (`NULL` if `uniq_id` was not supplied).

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
gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
                           p_mistake, p_missing, cond_hazard_params, TRUE )
prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
                           PIVs_config, TRUE, "entity_id", TRUE )
#> Warning: PIVs 'V4' and 'V3' are associated (Cramer's V = 0.31 in data1, 0.22 in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.
#> '2' is the larger source, saved as B; '1' saved as A.
str(prep_data, max.level = 1)
#> List of 6
#>  $ encodedA     :'data.frame':   150 obs. of  11 variables:
#>  $ encodedB     :'data.frame':   200 obs. of  9 variables:
#>  $ n_values     : Named int [1:4] 5 6 7 12
#>   ..- attr(*, "names")= chr [1:4] "V1" "V2" "V3" "V4"
#>  $ PIVs_config  :List of 4
#>  $ same_mistakes: logi TRUE
#>  $ true_pairs   :'data.frame':   100 obs. of  2 variables:
```
