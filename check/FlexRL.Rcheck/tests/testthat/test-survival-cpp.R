# The C++ survival kernel (survivalS_cpp) and its dispatcher survival_S().
#
# The built-in families must agree with the R closures in survival_model()
# to rounding level: R's %*% is dgemv (its own dot-product summation order)
# and R's sum() accumulates in extended precision, so ~1e-16 differences are
# expected and everything here asserts < 1e-12. Custom models must keep
# calling the user's R function, unchanged.

test_that("survivalS_cpp matches the R closures for every built-in family", {
  set.seed(1)
  n <- 500
  X <- cbind(1, rnorm(n), rnorm(n))
  times <- rexp(n)
  fams <- list(
    exponential = list(m = survival_model("exponential"),
                       alpha = c(log(0.4), 0.05, -0.1)),
    weibull     = list(m = survival_model("weibull"),
                       alpha = c(-0.3, log(0.4), 0.05, -0.1)),
    gompertz    = list(m = survival_model("gompertz"),
                       alpha = c(0.4, log(0.4), 0.05, -0.1)),
    piecewise   = list(m = survival_model("piecewise", cuts = c(0.7, 2)),
                       alpha = c(log(0.3), log(0.5), log(0.2), 0.05, -0.1))
  )
  for (fam in names(fams)) {
    m <- fams[[fam]]$m
    alpha <- fams[[fam]]$alpha
    s_cpp <- FlexRL:::survival_S(m, X, alpha, times)
    s_r <- m$S(X, alpha, times)
    expect_length(s_cpp, n)
    expect_lt(max(abs(s_cpp - s_r)), 1e-12,
              label = paste("family", fam))
  }
})

test_that("survivalS_cpp edge cases agree with the R closures", {
  # t = 0: S = 1 for every family (hazard has had no time to act)
  for (args in list(
    list(survival_model("exponential"), c(log(0.4), 0.1)),
    list(survival_model("weibull"), c(-0.3, log(0.4), 0.1)),
    list(survival_model("gompertz"), c(0.4, log(0.4), 0.1)),
    list(survival_model("piecewise", cuts = 1), c(log(0.3), log(0.5), 0.1))
  )) {
    s <- FlexRL:::survival_S(args[[1]], cbind(1, 0.5), args[[2]], 0)
    expect_equal(s, 1)
  }
  # gompertz near g = 0 falls back to the exponential form (|g| < 1e-8) and
  # still matches the R closure, which branches the same way
  m <- survival_model("gompertz")
  a <- c(5e-9, log(0.3), 0.1)
  expect_lt(max(abs(FlexRL:::survival_S(m, cbind(1, 0.5), a, 1.5) -
                    m$S(cbind(1, 0.5), a, 1.5))), 1e-12)
  # weibull with shape k = 1 (log-shape = 0) reduces to the exponential model
  w <- FlexRL:::survival_S(survival_model("weibull"), cbind(1, 0.5),
                           c(0, log(0.3), 0.1), 1.5)
  e <- FlexRL:::survival_S(survival_model("exponential"), cbind(1, 0.5),
                           c(log(0.3), 0.1), 1.5)
  expect_equal(w, e, tolerance = 1e-12)
  # piecewise: S is continuous at the cut points (both sides of 0.7 and 2)
  m <- survival_model("piecewise", cuts = c(0.7, 2))
  a <- c(log(0.3), log(0.5), log(0.2))
  eps <- 1e-10
  for (cut in c(0.7, 2)) {
    s_lo <- FlexRL:::survival_S(m, cbind(1), a, cut - eps)
    s_mid <- FlexRL:::survival_S(m, cbind(1), a, cut)
    s_hi <- FlexRL:::survival_S(m, cbind(1), a, cut + eps)
    expect_lt(abs(s_lo - s_mid), 1e-6)
    expect_lt(abs(s_mid - s_hi), 1e-6)
    # and matches the R closure exactly at the cut
    expect_lt(abs(s_mid - m$S(cbind(1), a, cut)), 1e-12)
  }
  # piecewise with covariates but no intercept-style beta (nCov = 1: lam = 1)
  m <- survival_model("piecewise", cuts = 1)
  expect_equal(FlexRL:::survival_S(m, cbind(1), c(log(0.3), log(0.5)), 1.5),
               exp(-(0.3 * 1 + 0.5 * 0.5)), tolerance = 1e-12)
})

test_that("survivalS_cpp is bitwise identical across thread counts", {
  set.seed(3)
  n <- 2000
  X <- cbind(1, rnorm(n), rnorm(n))
  times <- rexp(n)
  # weibull exercises dot + exp + pow; piecewise exercises the interval sum
  s1 <- FlexRL:::survivalS_cpp(X, c(-0.3, log(0.4), 0.05, -0.1), times,
                               "weibull", numeric(0), 1)
  s2 <- FlexRL:::survivalS_cpp(X, c(-0.3, log(0.4), 0.05, -0.1), times,
                               "weibull", numeric(0), 2)
  expect_identical(s1, s2)
  p1 <- FlexRL:::survivalS_cpp(X, c(log(0.3), log(0.5), log(0.2), 0.05, -0.1),
                               times, "piecewise", c(0.7, 2), 1)
  p2 <- FlexRL:::survivalS_cpp(X, c(log(0.3), log(0.5), log(0.2), 0.05, -0.1),
                               times, "piecewise", c(0.7, 2), 2)
  expect_identical(p1, p2)
})

test_that("custom models keep calling the user's R function", {
  calls <- 0L
  m <- survival_model("custom",
                      S = function(X, alpha, times) {
                        calls <<- calls + 1L
                        rep(0.42, length(times))
                      },
                      n_par = function(n_cov) n_cov,
                      init = function(n_cov) rep(0, n_cov))
  expect_equal(FlexRL:::survival_S(m, cbind(1, 2), 0, c(1, 2, 3)), c(0.42, 0.42, 0.42))
  expect_equal(calls, 1L)  # dispatched straight to model$S, once
  # and StEM keeps working end-to-end on the tiny problem, whose model is custom
  data <- tiny_data()
  args <- list(data = data, StEM_iter = 2, StEM_burnin = 1, gibbs_iter = 3,
               gibbs_burnin = 1, n_post_sample = 5)
  set.seed(11)
  fit <- suppressWarnings(do.call(StEM, args))
  expect_s3_class(fit$Delta, "data.frame")
})

test_that("StEM runs end-to-end with the C++ survival path (weibull, piecewise)", {
  data <- tiny_data()
  # gompertz is covered by the formula-agreement tests above; on the tiny
  # problem its M-step can drive g to a degenerate value (S collapses to
  # 0/1 with ~2 links), which is a property of the data, not the kernel.
  for (md in list(survival_model("weibull"),
                  survival_model("piecewise", cuts = 1))) {
    args <- list(data = data, StEM_iter = 2, StEM_burnin = 1, gibbs_iter = 3,
                 gibbs_burnin = 1, n_post_sample = 5, model_dynamics = md)
    set.seed(11)
    fit <- suppressWarnings(do.call(StEM, args))
    expect_s3_class(fit$Delta, "data.frame")
    expect_true(all(fit$Delta$x > 0 & fit$Delta$x <= 1))
  }
})
