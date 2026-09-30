# FlexRL

<!-- badges: start -->
[![](https://cranlogs.r-pkg.org/badges/grand-total/FlexRL)](https://cran.r-project.org/web/packages/FlexRL/index.html)
<!-- badges: end -->

Linking cohort studies, census data, surveys or administrative records opens research opportunities with existing data and broadens the capabilities of statistical offices. It enables the reuse of existing data and supports sharing data without unique identifiers, giving access to more variables across wider populations and longer time periods. Without unique identifiers, however, probabilistic record linkage makes errors that propagate into subsequent analyses.

FlexRL is an open-source R package for Flexible Record Linkage. It links records that refer to the same entities across two data sources without a unique identifier, using Partially Identifying Variables (PIVs) such as birth year, sex, postal code, or product code, brand, category. It applies wherever two sources are expected to overlap: health monitoring studies at two time points, survey waves, customer or product files of two retailers...

FlexRL implements the [Stochastic Expectation Maximisation (StEM) approach to record linkage](https://doi.org/10.1093/jrsssc/qlaf016) of Robach et al. (JRSS-C, 2025). The model accounts for registration errors (missing values and mistakes) and for dynamic PIVs that change over time (e.g. a postal code between two data collections), and it enforces one-to-one assignment. It returns the linked records together with their posterior linkage scores.

Since record linkage is rarely the end of the analysis, FlexRL also provides tools for inference on linked data. They make the quality of linked data transparent:
- [estimators of the false discovery proportion in record linkage](https://doi.org/10.1002/sim.70292) (Robach et al., Statistics in Medicine, 2025);
- [diagnostics for assessing the divergence between linked sample and source population](https://kayanerobach.github.io/blog/2025/causal-record-linkage/) (Robach et al., arXiv, 2026).

These tools apply to the output of any record linkage method, including other packages, so any study using linked data can report linkage quality alongside its results.

The StEM algorithm may take time to run on large data sets, but it has a low memory footprint and runs on a standard computer.

More details are in the [software article](...) (Robach et al., arXiv, 2026) and on [CRAN](https://cran.r-project.org/web/packages/FlexRL/index.html).

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

df2016 <- read.csv("FlexRL/vignettes/ex-reg-SHIW-16.csv", row.names = 1)
df2020 <- read.csv("FlexRL/vignettes/ex-reg-SHIW-20.csv", row.names = 1)

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

For support requests, contact _robachowyk@gmail.com_.
