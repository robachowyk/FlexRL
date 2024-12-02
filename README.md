
# FlexRL

<!-- badges: start -->
<!-- badges: end -->

FlexRL is a package for Flexible Record Linkage. Its goal is to find the common set of records among 2 data sources. The area of applications of Record Linkage is broad, it can be used to link data from healthcare monitoring studies, to count casualties in conflict zones by combining several registries, ... Therefore FlexRL models registration errors (missing values and mistakes in the data) and handles dynamic **P**artially **I**dentifying **V**ariable**s** that evolve over time (the zipcode can change between the 2 data collections for instance) in order to produce a final set of linked records with posterior probabilities to be linked. The algorithm can take time to run on large data sets but has a low memory footprint and can easily run on standard computers.

This package implements the **St**ochastic **EM** approach to Record Linkage described in '[https://arxiv.org/abs/2407.06835](A flexible model for Record Linkage)'. The main article and the supplementary material are available on arxiv.

## Installation

You can install the development version of FlexRL from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("robachowyk/FlexRL")
```

## How to use

This is a basic example which shows you how to solve a common problem:

``` r
library(FlexRL)
## basic example code
```

Documentation...

Vignettes...

FlexRL-experiments for reproducing the article experiments...
 
For support requests, contact _k dot c dot robach at amsterdamumc dot nl_.
