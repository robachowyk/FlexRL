# Regression tests for the tiny-input StEM crashes found in review.
#
# Two independent defects surfaced when StEM was run on very small files
# (which the rest of the suite never does):
#
# 1. sampleD() wrote through its LogicalVector& arguments. When a caller
#    passes a length-1 logical (the result of a scalar comparison, e.g.
#    `sumRowD > 0` with a single A record), R hands over its interned
#    TRUE/FALSE constant, so the write corrupted R's global constant:
#    R signalled "internal FALSE value has been modified" and afterwards
#    `1 > 2` evaluated to TRUE session-wide. Fixed by cloning the
#    arguments inside sampleD (src/FlexRLfunctions.cpp).
#
# 2. Matrix::Matrix(0, n, n, sparse = TRUE) returns a diagonalMatrix for
#    square inputs. Matrix::summary() of a diagonalMatrix is a diagSummary
#    (not the i/j/x triplet summary), and as.data.frame() on it fails with
#    "'dimnames' applied to non-array". This crashed StEM whenever
#    nrow(encodedA) == nrow(encodedB) and Delta stayed empty or
#    diagonal-only. Fixed by building Delta (and LLL) with
#    Matrix::sparseMatrix(..., x = numeric(0)), which is always a general
#    dgCMatrix.

test_that("StEM rejects hand-edited encoded data (NA in a PIV column) with a clear error", {
  # prepare_data() encodes missing values as 0; feeding it NA directly used to
  # crash inside the C++ Gibbs kernels ("Index out of bounds: index=2147483647",
  # from G(i)-1 overflowing on NA_INTEGER) instead of failing loudly.
  set.seed(1)
  cfg <- list(V1 = list(dynamics = "stable", bound_mistakes = c(NA, NA),
                        fix_mistakes = c(NA, NA)))
  gen <- simulate_data(cfg, c(4), c(5, 8), 3, list(V1 = c(0.02, 0.02)),
                       list(V1 = c(0.01, 0.01)), list(V1 = c()), TRUE)
  p <- prepare_data(gen$data1, gen$data2, "fileA", "fileB", cfg, uniq_id = "entity_id")
  p_bad <- p
  p_bad$encodedA$V1[c(1, 3)] <- NA
  expect_error(StEM(p_bad, StEM_iter = 2, StEM_burnin = 1, gibbs_iter = 3,
                    gibbs_burnin = 1, n_post_sample = 2),
               "0 = missing")
})

test_that("sampleD does not write through scalar logical arguments", {
  x <- rep(0, 1)
  FlexRL:::initDeltaMap()
  # One candidate pair and a near-certain link: the C++ code sets
  # sumRowD(0) = TRUE, which pre-fix corrupted R's interned FALSE.
  d <- FlexRL:::sampleD(
    matrix(1L, 1, 2), LLA = -1, LLB = -1, LLL = 100, gamma = 0.99,
    loglik = -10, nlinkrec = 0L, sumRowD = x > 0, sumColD = c(FALSE, FALSE, FALSE)
  )
  # Functionality preserved: the accepted link is reported back.
  expect_true(isTRUE(d$sumRowD[1]))
  expect_true(isTRUE(d$sumColD[1]))
  # The caller's own object is untouched and, crucially, R's global
  # interned constants are uncorrupted.
  expect_identical(x > 0, FALSE)
  expect_false(1 > 2)
})

test_that("StEM completes on tiny inputs, including square shapes", {
  for (shape in list(c(1, 5), c(2, 5), c(2, 2), c(3, 3))) {
    set.seed(1)
    cfg <- list(V1 = list(dynamics = "stable", bound_mistakes = c(NA, NA),
                          fix_mistakes = c(NA, NA)))
    gen <- simulate_data(cfg, c(4), shape, min(shape[1], 2),
                         list(V1 = c(0.02, 0.02)), list(V1 = c(0.01, 0.01)),
                         list(V1 = c()), TRUE)
    p <- prepare_data(gen$data1, gen$data2, "fileA", "fileB", cfg,
                      uniq_id = "entity_id")
    fit <- StEM(p, StEM_iter = 3, StEM_burnin = 1, gibbs_iter = 4,
                gibbs_burnin = 2, n_post_sample = 5)
    expect_s3_class(fit$Delta, "data.frame")
    expect_named(fit$Delta, c("i", "j", "x"))
    expect_true(all(fit$Delta$x > 0 & fit$Delta$x <= 1))
    # Corruption canary: pre-fix, a StEM run with a single A record flipped
    # R's interned FALSE, making this comparison TRUE for the rest of the
    # session (and wrecking the condition system).
    expect_false(1 > 2)
  }
})
