# Maximum Mean Discrepancy between two sets of variables

Joint (multivariate) measure of distributional discrepancy between `x`
and `y`, using a Gaussian (RBF) kernel with bandwidth set by the median
pairwise distance (biased estimator, Gretton et al. 2012).

## Usage

``` r
mmd(x, y)
```

## Arguments

- x, y:

  Numeric matrices with the same number of columns.

## Value

Numeric, the (biased) MMD estimate.

## Examples

``` r
x <- matrix(rnorm(50), ncol = 5)
y <- matrix(rnorm(30) + 0.5, ncol = 5)
mmd(x, y)
#> [1] 0.5729587
```
