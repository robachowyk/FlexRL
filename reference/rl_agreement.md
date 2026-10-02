# Agreement rate of pairs on shared variables

For a set of pairs, computes how often the two records agree on each
variable in `common_vars`.

## Usage

``` r
RL_agreement(data1, data2, common_vars, pairs, na.rm = TRUE, na_match = NULL)
```

## Arguments

- data1:

  Data frame containing `common_vars`.

- data2:

  Data frame containing `common_vars`.

- common_vars:

  Character vector, names of the variables to compare (must exist in
  both `data1` and `data2`).

- pairs:

  Data frame/matrix/list with 2 columns of indices (into `data1`,
  `data2`) for the pairs to evaluate.

- na.rm:

  Logical; if `TRUE` (default), pairs with a missing value on a variable
  are excluded from that variable's agreement rate.

- na_match:

  Logical, required if `na.rm = FALSE`: should a missing value be
  treated as agreeing (`TRUE`) or disagreeing (`FALSE`) with any value?

## Value

List with `agreements` (named numeric vector, one entry per variable in
`common_vars`) and, if available, `true_agreements`.

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
#> Running StEM algorithm ■■■■■■■■■■■■■■■■                  50% | iter 5/10 [1s]
#> Running StEM algorithm ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | iter 10/10 [1.8s]
#> 
linked_pairs <- fit$Delta[fit$Delta$x > 0.5, ]
RL_agreement( prep_data$encodedA, prep_data$encodedB,
              names(PIVs_config), linked_pairs )
#> $agreements
#>        V1        V2        V3        V4 
#> 0.9716312 0.9645390 0.9716312 0.7092199 
#> 
RL_agreement( prep_data$encodedA, prep_data$encodedB,
              names(PIVs_config), prep_data$true_pairs )
#> $agreements
#>    V1    V2    V3    V4 
#> 0.940 0.955 0.860 0.690 
#> 
```
