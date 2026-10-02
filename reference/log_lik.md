# Log-likelihood of the linkage matrix

Log-likelihood of the linkage matrix

## Usage

``` r
log_lik(LLL, LLA, LLB, links, sumRowD, sumColD, gamma)
```

## Arguments

- LLL:

  (Sparse) matrix of log-likelihood contributions for linked records.

- LLA:

  Numeric vector of log-likelihood contributions for non-linked records
  from A.

- LLB:

  Numeric vector of log-likelihood contributions for non-linked records
  from B.

- links:

  2-column matrix of indices (A, B) for the currently linked records.

- sumRowD:

  Logical vector, one entry per record in A: does it form a link?

- sumColD:

  Logical vector, one entry per record in B: does it form a link?

- gamma:

  Numeric, proportion of linked records as a fraction of the smaller
  file.

## Value

Numeric, the log-likelihood of the linkage matrix.

## Examples

``` r
LLL <- Matrix::Matrix(0, nrow = 13, ncol = 15, sparse = TRUE)
LLA <- stats::runif(13, 0, 2)
LLB <- stats::runif(15, 0, 2)
links <- as.matrix(data.frame(idxA = c(5, 9, 11, 12, 13),
                              idxB = c(5, 9, 11, 13, 15)))
LLL[links] <- 0.67
sumRowD <- (seq_len(13) %in% links[, 1])
sumColD <- (seq_len(15) %in% links[, 2])
gamma <- 0.5
log_lik(LLL, LLA, LLB, links, sumRowD, sumColD, gamma)
#> [1] 0.5597091
```
