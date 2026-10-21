# pasteIntoPattern: C++ pattern builder (guards FINDINGS #H).

test_that("pasteIntoPattern matches do.call(paste) on ordinary matrices", {
  A <- matrix(c(1L, 2L, 3L, 4L, 5L, 6L), nrow = 3)
  expect_equal(as.character(FlexRL:::pasteIntoPattern(A)),
               do.call(paste, c(as.data.frame(A), list(sep = "_"))))
})

test_that("pasteIntoPattern handles 0-column matrices without OOB access (#H)", {
  # A dataset with no stable PIVs produces a 0-column matrix here; the
  # pre-fix code read A(i, 0) out of bounds, yielding garbage patterns that
  # collapsed blocking to a spurious diagonal-only candidate set.
  zero_col <- matrix(integer(0), nrow = 4, ncol = 0)
  expect_equal(as.character(FlexRL:::pasteIntoPattern(zero_col)),
               rep("", 4))
  zero_col_0row <- matrix(integer(0), nrow = 0, ncol = 0)
  expect_length(FlexRL:::pasteIntoPattern(zero_col_0row), 0L)
})

test_that("a Gibbs run with no stable PIVs visits several link configs", {
  # With all patterns equal (""), every (i, j) is a candidate; the chain must
  # not collapse to the diagonal-only pathology seen pre-fix.
  set.seed(2024)
  source("helper-data.R", local = TRUE)
  pars <- tiny_pars()
  data <- tiny_data(pars); mod <- tiny_model()
  eta <- list(pars$eta); phi <- list(pars$phi); alpha <- list(0); gamma <- pars$gamma
  linksR <- matrix(integer(0), 0, 2)
  sumRowD <- rep(0L, 2L); sumColD <- rep(0L, 3L)
  survivalpSameH <- matrix(1, 0, 1)
  configs <- new.env(parent = emptyenv())
  linksCpp <- matrix(integer(0), 0, 2)
  for (s in 1:500) {
    h <- simulateH(data, linksCpp, survivalpSameH, sumRowD, sumColD, eta, phi)
    d <- simulateD(data, linksR, h$truepivsA, h$truepivsB,
                   gamma, eta, alpha, phi, mod)
    linksCpp <- d$links
    linksR <- linksCpp + 1L; sumRowD <- d$sumRowD; sumColD <- d$sumColD
    survivalpSameH <- matrix(1, nrow(linksR), 1)
    if (nrow(linksR) > 0) survivalpSameH[, 1] <- pars$S
    key <- if (nrow(linksR) == 0) "0" else
      paste(apply(linksR[order(linksR[, 1]), , drop = FALSE], 1,
                  function(r) paste(r, collapse = ">")), collapse = ";")
    configs[[key]] <- if (is.null(configs[[key]])) 1L else configs[[key]] + 1L
  }
  expect_gt(length(ls(configs)), 3L)
})
