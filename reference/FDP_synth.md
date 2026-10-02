# Model agnostic, based on synthetic data, false discovery proportion at a given threshold

Model agnostic, based on synthetic data, false discovery proportion at a
given threshold

## Usage

``` r
FDP_synth(idxA, idxB, LinkScore, threshold, n_records_A, n_records_B, n_synth)
```

## Arguments

- idxA:

  Integer vector, indices in A of linked records (for a previously set
  threshold).

- idxB:

  Integer vector, indices in B of linked records (for a previously set
  threshold).

- LinkScore:

  Numeric vector, linkage scores of the candidate pairs. If `NULL`,
  consider all given pairs as linked.

- threshold:

  Numeric, score threshold above which a pair is declared linked.

- n_records_A:

  Integer, number of records in A.

- n_records_B:

  Integer, number of records in B.

- n_synth:

  Integer, number of synthetic records to generate per iteration.

## Value

A list with `FDP_synth` Numeric, proportion of synthetic falsely linked
records, `n_linked_real` Integer, number of real linked records and
`n_linked_all` Integer, total number of linked records (real and
synthetic).

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
#> Warning: PIV 'V4': values out of the common support in data2 [5].
#> '2' is the larger source, saved as B; '1' saved as A.
PIVs <- names(PIVs_config)
n_synth <-as.integer(0.10 * nrow(prep_data$encodedB))
new_data <- synthesise("arf", prep_data$encodedA[, c(PIVs,"local_id","source",
                      "date",prep_data$PIVs_config$V4$cond_hazard_cov$covA)],
                      prep_data$encodedB[, c(PIVs,"local_id","source","date",
                      prep_data$PIVs_config$V4$cond_hazard_cov$covB)], PIVs, 
                      n_synth, TRUE)
#> Iteration: 0, Accuracy: 47.52%
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> Warning: PIV 'V1': synthetic records out of the common support were removed.
#> Warning: PIV 'V2': synthetic records out of the common support were removed.
#> Warning: PIV 'V4': synthetic records out of the common support were removed.
new_data$dataA[PIVs][is.na(new_data$dataA[PIVs])] <- 0
new_data$dataB[PIVs][is.na(new_data$dataB[PIVs])] <- 0
# cannot model dynamics for synthetic data, set dates to 0
new_data$dataA$date[is.na(new_data$dataA$date)] <- 0
new_data$dataB$date[is.na(new_data$dataB$date)] <- 0
arguments <- list(data = prep_data, StEM_iter = 5, StEM_burnin = 1, 
                  gibbs_iter = 5, gibbs_burnin = 2, n_post_samp = 5)
fit_flexrl <- link_with_FlexRL(new_data$dataA, new_data$dataB, arguments)  
#> FlexRL
FDP_synth(fit_flexrl$idxA, fit_flexrl$idxB, fit_flexrl$LinkScore, 0.5,
          nrow(prep_data$encodedA), nrow(prep_data$encodedB), 
          n_synth)
#> $FDP_synth
#> [1] 0.443997
#> 
#> $n_linked_real
#> [1] 137
#> 
#> $n_linked_all
#> [1] 143
#> 
```
