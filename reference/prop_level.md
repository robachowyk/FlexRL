# Proportion of a data set with a given variable at a given level

Proportion of a data set with a given variable at a given level

## Usage

``` r
prop_level(df, var, level)
```

## Arguments

- df:

  Data frame.

- var:

  Character, column name in `df`.

- level:

  Value to match against `df[[var]]`.

## Value

Numeric, proportion of rows where `df[[var]] == level` (na.rm=TRUE).

## Examples

``` r
df <- data.frame(colour = sample(c("orange", "purple"), 100, replace = TRUE))
prop_level(df, "colour", "purple")
#> [1] 0.49
```
