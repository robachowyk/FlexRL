# Gibbs conditional-distribution tests.
#
# These pin the C++ Gibbs kernels to the exact conditionals implied by the
# paper's model (arXiv:2407.06835, Eq. 4/12):
#   P(H^A = a, H^B = b | G, t, Delta = 1)
#     proportional to eta[a] * [S * 1{a=b} + (1-S)/(n-1) * 1{a!=b}]
#                    * P(G^A | a) * P(G^B | b)
# NOTE: eta is indexed by H^A (file A), per the likelihood in simulateD()
# (which uses eta[h^A]) and paper Eq. 4. Regression test for the eta-indexing
# bug found in sampleL/sampleH (FINDINGS #A).
#
# sampleNL (non-linked records) is a control: P(h | g) proportional to
# eta[h] * P(g | h) — this one was correct before the fix.

sampleL_oracle <- function(GA, GB, S, nval, agreeA, agreeB, missA, missB, eta) {
  # Grid rows follow sampleH: ExpandGrid(values, values), a-major:
  # row r (1-based) has a = (r-1) %/% nval + 1, b = (r-1) %% nval + 1.
  w <- rep(NA_real_, nval * nval)
  for (a in 1:nval) for (b in 1:nval) {
    pjoint <- if (a == b) S else (1 - S) / (nval - 1)
    pga <- if (GA == 0) missA
           else if (GA == a) (1 - missA) * agreeA
           else (1 - missA) * (1 - agreeA) / (nval - 1)
    pgb <- if (GB == 0) missB
           else if (GB == b) (1 - missB) * agreeB
           else (1 - missB) * (1 - agreeB) / (nval - 1)
    w[(a - 1) * nval + b] <- eta[a] * pjoint * pga * pgb  # eta on H^A
  }
  w / sum(w)
}

test_that("sampleL: joint (H^A, H^B) draw follows the model conditional", {
  p <- tiny_pars()
  nval <- p$n_values
  choice_set <- FlexRL:::ExpandGrid(FlexRL:::generateSequence(nval),
                                    FlexRL:::generateSequence(nval))
  choice_equal <- choice_set[, 1] == choice_set[, 2]
  # sampleH passes eta indexed by the grid's A column (post-Findings-#A;
  # `eta_choice = eta_k[choice_set(_,0) - 1]` in the C++): the marginal of
  # H^A is eta (paper Eq. 4), expanded per grid row so it multiplies
  # elementwise against probH inside sampleL. We replicate that convention:
  eta_choice <- p$eta[choice_set[, 1]]
  phiA <- c(p$phi["agreeA"], p$phi["missA"])
  phiB <- c(p$phi["agreeB"], p$phi["missB"])
  n_draws <- 40000

  ## Case A: registered value missing in A, observed (=1) in B.
  set.seed(101)
  draws <- replicate(n_draws,
                     FlexRL:::sampleL(0L, 1L, p$S, choice_set, choice_equal,
                                      nval, phiA, phiB, eta_choice))
  a <- (draws - 1L) %/% nval + 1L
  b <- (draws - 1L) %% nval + 1L
  oracle <- sampleL_oracle(0L, 1L, p$S, nval, p$phi["agreeA"], p$phi["agreeB"],
                           p$phi["missA"], p$phi["missB"], p$eta)
  expect_lt(abs(mean(a == 1) - sum(oracle[choice_set[, 1] == 1])), 0.02)
  expect_lt(abs(mean(b == 1) - sum(oracle[choice_set[, 2] == 1])), 0.02)
  # joint cell check
  expect_lt(abs(mean(a == 1 & b == 1) - oracle[1]), 0.02)

  ## Case B: registered value observed (=1) in A, missing in B.
  set.seed(102)
  draws2 <- replicate(n_draws,
                      FlexRL:::sampleL(1L, 0L, p$S, choice_set, choice_equal,
                                       nval, phiA, phiB, eta_choice))
  a2 <- (draws2 - 1L) %/% nval + 1L
  b2 <- (draws2 - 1L) %% nval + 1L
  oracle2 <- sampleL_oracle(1L, 0L, p$S, nval, p$phi["agreeA"], p$phi["agreeB"],
                            p$phi["missA"], p$phi["missB"], p$eta)
  expect_lt(abs(mean(b2 == 2) - sum(oracle2[choice_set[, 2] == 2])), 0.02)
  expect_lt(abs(mean(a2 == 1) - sum(oracle2[choice_set[, 1] == 1])), 0.02)
})

test_that("sampleNL: true value of a non-linked record follows eta * P(G|H)", {
  p <- tiny_pars()
  nval <- p$n_values
  phiA <- c(p$phi["agreeA"], p$phi["missA"])
  n_draws <- 40000
  set.seed(103)
  # G = 0 (missing): H must follow eta exactly
  draws0 <- FlexRL:::sampleNL(rep(0L, n_draws), p$eta, phiA, runif(n_draws))
  expect_lt(abs(mean(draws0 == 1) - p$eta[1]), 0.02)
  # G = 1: P(h=1) proportional to eta[1] * P(G=1|h=1)
  p1 <- p$eta[1] * (1 - p$phi["missA"]) * p$phi["agreeA"]
  p2 <- p$eta[2] * (1 - p$phi["missA"]) * (1 - p$phi["agreeA"]) / (nval - 1)
  draws1 <- FlexRL:::sampleNL(rep(1L, n_draws), p$eta, phiA, runif(n_draws))
  expect_lt(abs(mean(draws1 == 1) - p1 / (p1 + p2)), 0.02)
})
