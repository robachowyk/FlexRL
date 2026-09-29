# FlexRL

<!-- badges: start -->
[![](https://cranlogs.r-pkg.org/badges/grand-total/FlexRL)](https://cran.r-project.org/web/packages/FlexRL/index.html)
<!-- badges: end -->

FlexRL is an R package for Flexible Record Linkage: it probabilistically link records that refer to the same entities across two data sources without a unique identifier, using Partially Identifying Variables (PIVs) such as product code, brand, category, birth year, sex or postal code. It applies wherever two sources are expected to overlap: healthcare monitoring studies at two time points, registries of casualties in conflict zones collected by distinct organisations, customer or product files of two retailers, survey waves, ...

FlexRL implements the [Stochastic Expectation Maximisation (StEM) approach to record linkage](https://doi.org/10.1093/jrsssc/qlaf016) of Robach et al. (2025). The model accounts for registration errors (missing values and mistakes) and for dynamic PIVs that evolve over time (e.g. postal code may change between data collections), and enforces one-to-one assignment. It returns the set of linked records together with their posterior linkage scores.

Since record linkage is rarely the end of the analysis, the package also provides tools for inference on the linked data: [estimators of the false discovery proportion of a linkage](https://doi.org/10.1002/sim.70292) and diagnostics comparing the linked sample with the source data. These tools also apply to the linkage output of other record linkage packages.

The algorithm can take time to run on large data sets, but it has a low memory footprint and runs on a standard computer.

Please [open an issue](https://github.com/robachowyk/FlexRL/issues) to report any bug, to make a request, or to ask for help :-)

## Installation

You can install FlexRL from CRAN with:

```r
install.packages("FlexRL")
library(FlexRL)
```

Or you can install the development version of FlexRL from its [GitHub](https://github.com/robachowyk/FlexRL) with one of the following:

```r
pak::pak("robachowyk/FlexRL")
```

```r
remotes::install_github("robachowyk/FlexRL")
```

```r
devtools::install_github("robachowyk/FlexRL")
```

FlexRL relies on Rcpp; when imported from Github, it may require gfortran and gcc.

## How to use `FlexRL`

A minimal example shipped with the package vignettes:

```r
library(FlexRL)

df2016 <- read.csv("FlexRL/vignettes/exSHIW16.csv", row.names = 1)
df2020 <- read.csv("FlexRL/vignettes/exSHIW20.csv", row.names = 1)

# one entry per PIV: stable (does not change over time), flexible (may change, change not modelled) or structured (change modelled with a survival model)
PIVs_config <- list(
  ANASCI = list(dynamics = "stable",   bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  SESSO  = list(dynamics = "stable",   bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  STACIV = list(dynamics = "flexible", bound_mistakes = c(NA, NA),     fix_mistakes = c(NA, NA)),
  STUDIO = list(dynamics = "flexible", bound_mistakes = c(NA, NA),     fix_mistakes = c(NA, NA))
)
PIVs <- names(PIVs_config)

# encode the PIVs, restrict to the common support, put the smaller file in A
prep_data <- prepare_data(df2016, df2020, "2016", "2020", PIVs_config,
                         same_mistakes = TRUE, uniq_id = "ID")

# gauge the difficulty of the task with exact matching
naive_linkage(PIVs, prep_data$encodedA, prep_data$encodedB)

# fit the model
fit <- StEM(data = prep_data, StEM_iter = 20, StEM_burnin = 10,
            gibbs_iter = 20, gibbs_burnin = 10)

# linked pairs: rows of encodedA (i), rows of encodedB (j), linkage score (x); a score > 0.5 enforces one-to-one assignment
linked <- fit$Delta[fit$Delta$x > 0.5, ]

# diagnostics for inference on the linked data
diag <- RL_diagnostics(fit, prep_data$encodedA, prep_data$encodedB, PIVs,
                       list(ANASCI = TRUE, SESSO = FALSE, STACIV = FALSE, STUDIO = FALSE),
                       true_pairs = prep_data$true_pairs, FDP_estimation = TRUE, 
                       RL_method = "FlexRL", data = prep_data,
                       StEM_iter = 10, StEM_burnin = 5, 
                       gibbs_iter = 10, gibbs_burnin = 5,
                       maxIter4CV = 3, n_repeats = 5)
diag # print(diag)
print(diag, threshold = 0.75)
plot(diag, "scores")
plot(diag, "distributions", threshold = 0.75)
plot(diag, "convergence")
plot(diag, "FDP")
plot(diag, "discrepancy")
```

More documentation is available on [CRAN](https://cran.r-project.org/web/packages/FlexRL/index.html) and in the repository [FlexRL-experiments](https://github.com/robachowyk/FlexRL-experiments).

For support requests, contact _robachowyk@gmail.com_.
