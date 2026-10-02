# Intersection-over-union of two histogram supports

Intersection-over-union of two histogram supports

## Usage

``` r
support_iou(h1, h2)
```

## Arguments

- h1:

  `histogram` object (as returned by
  [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html)).

- h2:

  `histogram` object (as returned by
  [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html)).

## Value

Numeric, IoU of the ranges over which `h1` and `h2` have non-zero
counts.

## Examples

``` r
h1 <- graphics::hist(rnorm(200), plot = FALSE)
h2 <- graphics::hist(rnorm(200) + 1, plot = FALSE)
support_iou(h1, h2)
#> [1] 0.6923077
```
