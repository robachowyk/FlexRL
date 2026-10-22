# Likelihood helpers and FDP estimators vs the papers' equations.

test_that("log_possible_config counts arrangements n_B!/(n_B-k)!", {
  expect_equal(FlexRL:::log_possible_config(5, 0), 0)
  expect_equal(FlexRL:::log_possible_config(5, 2), log(5 * 4))
  expect_equal(FlexRL:::log_possible_config(5, 5), log(factorial(5)))
  expect_equal(FlexRL:::log_possible_config(5, 2),
               lgamma(6) - lgamma(4))
})

test_that("log_lik matches a hand-computed complete-data log-likelihood", {
  # 2 A-records, 3 B-records, one link (A1-B1); LLL is an nA x nB matrix
  LLL <- matrix(-0.3, 2, 3)
  LLA <- c(-0.2, -0.4)                 # per A record (used when not linked)
  LLB <- c(-0.1, -0.5, -0.6)           # per B record (used when not linked)
  links <- matrix(c(1, 1), nrow = 1)
  sumRowD <- c(1, 0)
  sumColD <- c(1, 0, 0)
  gamma <- 0.5
  expected <- (log(0.5) * 1 + log(0.5) * 1) - FlexRL:::log_possible_config(3, 1) +
    (-0.4) + (-0.5 - 0.6) + (-0.3)
  expect_equal(log_lik(LLL, LLA, LLB, links, sumRowD, sumColD, gamma), expected)
})

test_that("log_lik warns when B is smaller than the number of links", {
  expect_warning(
    log_lik(matrix(-0.3, 2, 2), c(-0.2, -0.2), c(-0.1, -0.1),
            matrix(c(1, 1, 2, 1, 1, 2), nrow = 3, byrow = TRUE),
            c(1, 1), c(1, 1), 0.5),
    "fewer records"
  )
})

test_that("FDP_score equals 1 - mean(score | score > threshold) (paper Eq. 1)", {
  scores <- c(0.9, 0.8, 0.3)
  expect_equal(FDP_score(scores, 0.5)$FDP_score, 1 - mean(c(0.9, 0.8)))
  expect_equal(FDP_score(scores, 0.5)$n_linked, 2L)
  expect_equal(FDP_score(scores, 0.85)$FDP_score, 0.1)
  expect_equal(FDP_score(scores, 0.85)$n_linked, 1L)
  # no links above threshold -> NaN (documented behaviour)
  expect_true(is.nan(FDP_score(scores, 0.95)$FDP_score))
  expect_equal(FDP_score(scores, 0.95)$n_linked, 0L)
})

test_that("FDP_synth implements paper Eq. 2: FPsynth * (N_B/N_synth) / N_real_linked", {
  # 10 real A, 100 real B, 10 synthetic appended to B (indices 101..110)
  idxA <- c(1, 2, 3)
  idxB <- c(5, 105, 7)                 # one synthetic false positive
  res <- FDP_synth(idxA, idxB, c(0.9, 0.8, 0.7), 0.5, 10, 100, 10)
  expect_equal(res$FDP_synth, (1 * (100 / 10)) / (3 - 1))
  expect_equal(res$n_linked_real, 2L)
  expect_equal(res$n_linked_all, 3L)
})

test_that("FDP_synth handles LinkScore = NULL (score-less methods)", {
  idxA <- c(1, 2); idxB <- c(101, 102)
  res <- FDP_synth(idxA, idxB, NULL, NULL, 10, 100, 10)
  expect_equal(res$n_linked_all, 2L)
  expect_equal(res$n_linked_real, 0L)
  expect_equal(res$FDP_synth, (2 * 10) / 0)   # Inf: no real links at all
})

test_that("FDP_synth classifies real vs synthetic via idxB <= n_records_B", {
  idxA <- c(1, 1); idxB <- c(100, 101)
  res <- FDP_synth(idxA, idxB, NULL, NULL, 10, 100, 10)
  expect_equal(res$n_linked_real, 1L)  # idxB = 100 is the last real record
})
