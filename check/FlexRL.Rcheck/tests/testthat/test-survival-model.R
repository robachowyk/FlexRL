# survival_model(): closed-form survival functions, parameter bookkeeping,
# and the custom-model interface.

test_that("exponential survival is exp(-lambda * t)", {
  m <- survival_model("exponential")
  expect_equal(m$S(cbind(1), log(0.5), 2), exp(-1))
  # X = intercept + covariate, one ROW per linked pair; hazard = exp(X alpha):
  # record 1 (x = 1): hazard exp(log .5 + log 2) = 1, S = exp(-3);
  # record 2 (x = 2): hazard exp(log .5 + 2 log 2) = 2, S = exp(-2 * 1).
  expect_equal(m$S(cbind(1, c(1, 2)), c(log(0.5), log(2)), c(3, 1)),
               c(exp(-3), exp(-2)))
  expect_equal(m$n_par(3), 3L)
  expect_length(m$init(3), 3)
  expect_identical(m$par_names(c("intercept", "Xe", "Xf")), c("intercept", "Xe", "Xf"))
})

test_that("weibull survival is exp(-(lambda * t)^k)", {
  m <- survival_model("weibull")
  alpha <- c(log(2), log(0.5))          # log shape, log lambda
  expect_equal(m$S(cbind(1), alpha, 3), exp(-(0.5 * 3)^2))
  expect_equal(m$n_par(2), 3L)
  expect_length(m$init(2), 3)
})

test_that("gompertz survival is exp(-lambda * (exp(g t) - 1) / g)", {
  m <- survival_model("gompertz")
  alpha <- c(0.5, log(0.3))             # g, log lambda
  expect_equal(m$S(cbind(1), alpha, 2), exp(-0.3 * (exp(0.5 * 2) - 1) / 0.5))
  # g -> 0 limit must be continuous: exp(-lambda * t)
  m2 <- survival_model("gompertz")
  expect_equal(m2$S(cbind(1), c(1e-12, log(0.3)), 2), exp(-0.3 * 2), tolerance = 1e-6)
  expect_equal(m$n_par(2), 3L)
})

test_that("piecewise survival integrates interval hazards", {
  m <- survival_model("piecewise", cuts = c(1, 2))
  alpha <- c(log(0.5), log(1.5), log(3))  # log hazards of the 3 intervals
  # t = 1.5: 1 time unit at h=0.5, 0.5 time units at h=1.5 -> H0 = 1.25
  expect_equal(m$S(cbind(1), alpha, 1.5), exp(-1.25))
  # t beyond the last cut: full cumulative hazard
  expect_equal(m$S(cbind(1), alpha, 4), exp(-(0.5 * 1 + 1.5 * 1 + 3 * 2)))
  expect_equal(m$n_par(1), 3L)
  expect_error(survival_model("piecewise"), "cuts")
  expect_error(survival_model("piecewise", cuts = c(2, 1)), "increasing")
})

test_that("custom models pass through and negloglik clamps degenerate S", {
  m <- survival_model("custom",
                      S = function(X, alpha, times) rep(0.42, length(times)),
                      n_par = function(n_cov) 1L,
                      init = function(n_cov) 0)
  expect_equal(m$S(cbind(1), 0, 5), 0.42)
  expect_identical(m$par_names(c("intercept")), "alpha1")
  expect_error(survival_model("custom", S = NULL, n_par = function(n) 1L, init = function(n) 0), "functions")

  # negloglik must stay finite when S hits 0 or exceeds 1 (clamping)
  m_bad <- survival_model("custom",
                          S = function(X, alpha, times) rep(0, length(times)),
                          n_par = function(n_cov) 1L, init = function(n_cov) 0)
  expect_true(is.finite(m_bad$negloglik(0, cbind(1), c(1, 2), c(TRUE, FALSE))))
  m_over <- survival_model("custom",
                           S = function(X, alpha, times) rep(1.5, length(times)),
                           n_par = function(n_cov) 1L, init = function(n_cov) 0)
  expect_true(is.finite(m_over$negloglik(0, cbind(1), c(1, 2), c(TRUE, FALSE))))
})

test_that("negloglik matches a hand-computed Bernoulli likelihood", {
  m <- survival_model("exponential")
  alpha <- log(0.5)
  times <- c(1, 2)
  Hequal <- c(TRUE, FALSE)
  S <- c(exp(-0.5), exp(-1))
  expected <- -(log(S[1]) + log(1 - S[2]))
  expect_equal(m$negloglik(alpha, cbind(1), times, Hequal), expected)
})
