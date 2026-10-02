# Standardised mean difference between a selected set and a baseline

When missing values are encoded as 0m the smd will report information on
the missingness in the linked sample vs. the source.

## Usage

``` r
smd(data_select, data_baseline, var, continuous = TRUE)
```

## Arguments

- data_select:

  Data frame to compare (e.g. the linked set vs. the original file).

- data_baseline:

  Data frame to compare (e.g. the linked set vs. the original file).

- var:

  Character, column name to compare.

- continuous:

  Logical; if `TRUE`, compares means of `var` directly; if `FALSE`,
  compares the proportion at each observed level of `var` (default
  `TRUE`).

## Value

Named list, one SMD per variable (`continuous = TRUE`) or per level
(`continuous = FALSE`).

## Examples

``` r
base <- data.frame(age = rnorm(200, 40, 10), sex = sample(c("M", "F"),
                    200, replace = TRUE))
select <- data.frame(age = rnorm(80, 43, 10), sex = sample(c("M", "F"), 80,
                    replace = TRUE, prob = c(0.6, 0.4)))
smd(select, base, "age")
#> $age
#> [1] 0.2454256
#> 
smd(select, base, "sex", continuous = FALSE)
#> $sex_F
#> [1] -0.164884
#> 
#> $sex_M
#> [1] 0.164884
#> 
```
