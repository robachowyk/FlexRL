# Survival models for the dynamics of a PIV

A structured PIV may change between the two registrations. The
probability that the true value of a linked pair is unchanged after a
time gap `t` is modelled by a survival function `S(t | X, alpha)`, where
`X` holds an intercept and the covariates given in `cond_hazard_cov` and
`alpha` the parameters estimated by [`StEM()`](stEM.md). This function
builds the model object used by [`StEM()`](stEM.md),
[`simulateD()`](simulateD.md) and [`simulate_data()`](simulate_data.md).

## Usage

``` r
survival_model(
  type = c("exponential", "weibull", "gompertz", "piecewise", "custom"),
  cuts = NULL,
  S = NULL,
  n_par = NULL,
  init = NULL
)
```

## Arguments

- type:

  One of `"exponential"`, `"weibull"`, `"gompertz"`, `"piecewise"`,
  `"custom"`.

- cuts:

  Numeric vector of cut points for `type = "piecewise"` (e.g. `c(1, 3)`
  gives three intervals `[0,1)`, `[1,3)`, `[3, Inf)`).

- S:

  For `type = "custom"`: function `(X, alpha, times)` returning the
  survival probability of each linked pair; `X` is a matrix with an
  intercept column first, then the covariates.

- n_par:

  For `type = "custom"`: function of `ncol(X)` returning the length of
  `alpha`.

- init:

  For `type = "custom"`: function of `ncol(X)` returning the starting
  values of `alpha`.

## Value

An object of class `"survival_model"`: a list with `type`, `S`,
`negloglik`, `n_par`, `init` and `par_names`.

## Details

Available models (`h` is the hazard, `lambda = exp(X alpha_cov)` the
proportional-hazards term):

- `"exponential"`:

  `S(t) = exp(-lambda t)`; `alpha` = coefficients of `X` (default, the
  model of the methodology paper).

- `"weibull"`:

  `S(t) = exp(-(lambda t)^k)` with shape `k = exp(alpha[1])`; `alpha` =
  log-shape, then coefficients of `X`.

- `"gompertz"`:

  `h(t) = lambda exp(g t)`, so `S(t) = exp(-lambda (exp(g t) - 1) / g)`
  with `g = alpha[1]`; `alpha` = `g`, then coefficients of `X`.

- `"piecewise"`:

  piecewise-constant baseline hazard on the intervals defined by `cuts`,
  times `lambda`; `alpha` = one log-hazard per interval, then
  coefficients of the covariates (no intercept).

- `"custom"`:

  `S`, `n_par` and `init` supplied by the user.

For every model the negative log-likelihood used in the M-step is
`-sum(Hequal log S + (1 - Hequal) log(1 - S))`, where `Hequal` indicates
that the true values of the linked pair agree. Parameters are estimated
with [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) (numerical
gradient).

## Examples

``` r
X <- cbind(intercept = rep(1, 5))
times <- c(0.001, 0.2, 1.3, 1.5, 2)
Hequal <- c(TRUE, TRUE, TRUE, FALSE, FALSE)

expo <- survival_model("exponential")
expo$S(X, alpha = log(0.3), times)
#> [1] 0.9997000 0.9417645 0.6770569 0.6376282 0.5488116
stats::nlminb(expo$init(ncol(X)), expo$negloglik, X = X, times = times,
              Hequal = Hequal)$par
#> [1] -0.3719637

weib <- survival_model("weibull")
weib$S(X, alpha = c(0.5, -1), times)
#> [1] 0.9999978 0.9865530 0.7435143 0.6871328 0.5471930
stats::nlminb(weib$init(ncol(X)), weib$negloglik, X = X, times = times,
              Hequal = Hequal)$par
#> [1]  5.3792727 -0.3564221

gomp <- survival_model("gompertz")
gomp$S(X, alpha = c(0.5, -1), times)
#> [1] 0.9996321 0.9255377 0.5098609 0.4396208 0.2824536
stats::nlminb(gomp$init(ncol(X)), gomp$negloglik, X = X, times = times,
              Hequal = Hequal)$par
#> [1]  112.3277 -160.6559
              
pw <- survival_model("piecewise", cuts = 1)
pw$S(X, alpha = c(-2, -1), times)
#> [1] 0.9998647 0.9732960 0.7821575 0.7266757 0.6045840
stats::nlminb(pw$init(ncol(X)), pw$negloglik, X = X, times = times,
              Hequal = Hequal)$par
#> [1] -22.4570466   0.9000345

# a custom model: log-logistic with scale exp(X alpha[-1]) 
# and shape exp(alpha[1])
loglogistic <- survival_model("custom",
  S = function(X, alpha, times) {
    lambda <- exp(as.matrix(X) %*% alpha[-1])
    1 / ( 1 + (lambda * times)^exp(alpha[1]) )
  },
  n_par = function(n_cov) n_cov + 1,
  init = function(n_cov) c(stats::rnorm(1, 0, 0.1), stats::runif(n_cov, log(0.01), log(1))))
loglogistic$S(X, alpha = c(0, -2), times)
#>           [,1]
#> [1,] 0.9998647
#> [2,] 0.9736463
#> [3,] 0.8503865
#> [4,] 0.8312532
#> [5,] 0.7869860
stats::nlminb(loglogistic$init(ncol(X)), loglogistic$negloglik, X = X, times = times,
              Hequal = Hequal)$par
#> [1]  5.6941282 -0.3342964
```
