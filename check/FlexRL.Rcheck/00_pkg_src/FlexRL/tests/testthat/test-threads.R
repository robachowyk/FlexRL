# OpenMP thread-count invariance of the Gibbs kernels (Phase 3b).
#
# pairLogLik()'s parallel loops give every pair its own output slot and no
# reduction crosses threads, so results must be BITWISE identical for any
# `threads` value. Running with 2 threads also proves the parallel regions
# execute when the package is built with OpenMP: a data race would corrupt
# the chain and break identity. When the package is built without OpenMP the
# test passes trivially (both runs sequential).
#
# `threads = 2` on purpose: CRAN asks that tests use at most 2 cores.

test_that("StEM output is bitwise identical for threads = 1 and threads = 2", {
  std <- std_generate(seed = 21)
  set.seed(63)
  fit1 <- suppressWarnings(do.call(StEM, c(list(data = std$prep, threads = 1), std_stem_args())))
  set.seed(63)
  fit2 <- suppressWarnings(do.call(StEM, c(list(data = std$prep, threads = 2), std_stem_args())))
  expect_identical(fit1$Delta, fit2$Delta)
  expect_identical(fit1$gamma, fit2$gamma)
  expect_identical(fit1$eta, fit2$eta)
  expect_identical(fit1$alpha, fit2$alpha)
  expect_identical(fit1$phi, fit2$phi)
})
