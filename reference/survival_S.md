# Evaluate a survival model's S, C++ fast path for the built-in families

Internal dispatcher used by the Gibbs hot path
([`StEM()`](https://robachowyk.github.io/FlexRL/reference/StEM.md)
through
[`simulateD()`](https://robachowyk.github.io/FlexRL/reference/simulateD.md)
and the per-sweep `survivalpSameH` refresh): built-in models evaluate
through the fused, OpenMP-parallel `survivalS_cpp()` kernel;
`type = "custom"` models keep calling the user's R function, unchanged.

## Usage

``` r
survival_S(model, X, alpha, times, threads = 4)
```

## Arguments

- model:

  Object built by
  [`survival_model()`](https://robachowyk.github.io/FlexRL/reference/survival_model.md).

- X:

  Covariate matrix, intercept column first (same contract as `model$S`).

- alpha:

  Parameter vector (layout per family, see
  [`survival_model()`](https://robachowyk.github.io/FlexRL/reference/survival_model.md)).

- times:

  Numeric vector of time gaps.

- threads:

  Integer, number of OpenMP threads for the C++ kernel (default 4;
  ignored without OpenMP; identical results for any value).

## Value

Numeric vector of survival probabilities, one per row of `X`.
