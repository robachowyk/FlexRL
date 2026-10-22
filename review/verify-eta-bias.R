#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# FINDINGS #A in isolation: the eta-indexing bias of the linked-pair H-sampler.
#
# The upstream master build (@ fe9551c) samples the joint (H^A, H^B) of LINKED
# pairs with eta indexed by H^B (simulateH's `eta_choice =
# eta_k[choice_set(_,1) - 1]`), while the model (paper Eq. 4, and simulateD's
# own likelihood) indexes eta by H^A. This biases the H-sampler whenever the
# linked-pair support is not diagonal-only, i.e. for structured PIVs with
# survival probability < 1 (for stable/flexible PIVs S = 1, only the diagonal
# (h, h) carries weight, and there eta[h^B] == eta[h^A] — the bug is a no-op;
# that is why the master-vs-baseline seed clouds of the real-data problems are
# bitwise identical, see PARITY-REPORT.md).
#
# This script calls FlexRL:::sampleL DIRECTLY (no candidate-pair pipeline
# involved — that is a separate master defect, the missing 0-column guard of
# pasteIntoPattern), for the two canonical tiny-problem cases, and compares the
# empirical joint distribution of (H^A, H^B) against the exact conditional
#   P(a, b | G, Delta=1) ∝ eta[a] * [S * 1{a=b} + (1-S)/(n-1) * 1{a!=b}]
#                        * P(G^A | a) * P(G^B | b)
# (same oracle as tests/testthat/test-gibbs-conditionals.R).
#
# Arms (per case):
#   baseline builds (master @ fe9551c, glm-notes @ 4b09088): sampleL takes a
#     pre-expanded eta_choice grid, so BOTH indexings are evaluated on the SAME
#     kernel:
#       eta[h^B] = master's shipped simulateH behaviour (the bug),
#       eta[h^A] = the corrected behaviour (what glm-notes ships).
#   optimised build (glm-parallel): sampleL takes eta directly and indexes by
#     H^A internally (the fix, factorised kernel).
#
# Expected: eta[h^A] and the optimised kernel match the oracle to Monte-Carlo
# noise; eta[h^B] deviates by the bias. On the baseline kernel the two arms
# consume the same RNG stream, so their difference is purely the weighting.
#
# Usage: R_LIBS=<lib> Rscript review/verify-eta-bias.R [n_draws] [outdir]
# Default outdir: ../FlexRL-bench/parity next to the repo.
# ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
n_draws <- if (length(args) >= 1) as.integer(args[1]) else 200000L
a <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", a, value = TRUE)
repo <- if (length(f)) normalizePath(file.path(dirname(sub("^--file=", "", f[1])), "..")) else normalizePath(".")
outdir <- if (length(args) >= 2) args[2] else file.path(repo, "..", "FlexRL-bench", "parity")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# Optional label for the output file (which build R_LIBS selected); the
# baseline-kernel arms are build-independent (sampleL is identical in master
# and glm-notes), so the label only documents provenance.
label <- if (length(args) >= 3) args[3] else if ("unif" %in% names(formals(FlexRL:::sampleNL))) "new" else "baseline"

library(FlexRL)
is_new <- "unif" %in% names(formals(FlexRL:::sampleNL))

## ------------------------------------------------------------- tiny problem
## Same parameters as tests/testthat/helper-data.R tiny_pars():
## one structured PIV, 2 values, skewed eta, custom constant-S model.
nval <- 2L
eta <- c(0.7, 0.3)
S <- 0.6
phiA <- c(agreeA = 0.9, missA = 0.10)   # [agree, missing] for file A records
phiB <- c(agreeB = 0.85, missB = 0.05)  # [agree, missing] for file B records

#' Exact conditional joint of (H^A, H^B) for one linked pair (a-major grid,
#' same layout as sampleH's ExpandGrid: row r has a = (r-1) %/% nval + 1).
sampleL_oracle <- function(GA, GB) {
  w <- rep(NA_real_, nval * nval)
  for (aa in 1:nval) for (bb in 1:nval) {
    pjoint <- if (aa == bb) S else (1 - S) / (nval - 1)
    pga <- if (GA == 0) phiA["missA"]
           else if (GA == aa) (1 - phiA["missA"]) * phiA["agreeA"]
           else (1 - phiA["missA"]) * (1 - phiA["agreeA"]) / (nval - 1)
    pgb <- if (GB == 0) phiB["missB"]
           else if (GB == bb) (1 - phiB["missB"]) * phiB["agreeB"]
           else (1 - phiB["missB"]) * (1 - phiB["agreeB"]) / (nval - 1)
    w[(aa - 1) * nval + bb] <- eta[aa] * pjoint * pga * pgb  # eta on H^A
  }
  w / sum(w)
}

## -------------------------------------------------------------------- arms
## Case A: registered value missing (0) in A, observed (=1) in B.
## Case B: registered value observed (=1) in A, missing (0) in B.
cases <- list(
  list(name = "GA=0, GB=1", GA = 0L, GB = 1L),
  list(name = "GA=1, GB=0", GA = 1L, GB = 0L)
)

run_baseline_arm <- function(GA, GB, eta_index) {
  # eta_index: 1 = index by H^A (corrected), 2 = index by H^B (master's
  # shipped simulateH behaviour). Same kernel, same RNG stream for both.
  values <- FlexRL:::generateSequence(nval)
  choice_set <- FlexRL:::ExpandGrid(values, values)
  choice_equal <- choice_set[, 2] == choice_set[, 1]
  eta_choice <- eta[choice_set[, eta_index]]  # R is 1-based: eta[value]
  out <- FlexRL:::sampleL(rep(GA, n_draws), rep(GB, n_draws), rep(S, n_draws),
                          choice_set, choice_equal, nval, phiA, phiB, eta_choice)
  cbind(a = choice_set[out, 1], b = choice_set[out, 2])
}

run_new_arm <- function(GA, GB) {
  draws <- FlexRL:::sampleL(rep(GA, n_draws), rep(GB, n_draws), rep(S, n_draws),
                            nval, phiA, phiB, eta, runif(2 * n_draws))
  draws[, 1:2]
}

emp <- function(ab) {
  # empirical joint over the 4 cells, a-major (matches the oracle layout)
  n <- nrow(ab)
  v <- tabulate((ab[, 1] - 1L) * nval + ab[, 2], nbins = nval * nval) / n
  names(v) <- sprintf("a=%d,b=%d", rep(1:nval, each = nval), rep(1:nval, times = nval))
  v
}

res <- list(variant = if (is_new) "new (glm-parallel)" else "baseline",
            label = label,
            r_version = R.version.string, n_draws = n_draws)
for (cs in cases) {
  oracle <- sampleL_oracle(cs$GA, cs$GB)
  names(oracle) <- sprintf("a=%d,b=%d", rep(1:nval, each = nval), rep(1:nval, times = nval))
  arms <- list()
  if (is_new) {
    set.seed(2024)
    arms[["new kernel (eta[h^A] internal)"]] <- emp(run_new_arm(cs$GA, cs$GB))
  } else {
    set.seed(2024)
    arms[["baseline kernel, eta[h^B] (master's shipped simulateH)"]] <-
      emp(run_baseline_arm(cs$GA, cs$GB, 2L))
    set.seed(2024)
    arms[["baseline kernel, eta[h^A] (corrected, glm-notes' simulateH)"]] <-
      emp(run_baseline_arm(cs$GA, cs$GB, 1L))
  }
  res[[cs$name]] <- list(oracle = oracle, arms = arms)
}

## ------------------------------------------------------------------ report
say <- function(...) cat(sprintf(...), "\n", sep = "")
say("eta-indexing bias of the linked-pair H-sampler (FINDINGS #A), in isolation")
say("build: %s | %s | %d draws per arm", res$variant, res$r_version, n_draws)
for (cs in cases) {
  say("")
  say("== case %s (S = %.2f, eta = c(%.2f, %.2f))", cs$name, S, eta[1], eta[2])
  say("   exact oracle joint: %s",
      paste(sprintf("%s=%.4f", names(res[[cs$name]]$oracle), res[[cs$name]]$oracle), collapse = "  "))
  for (nm in names(res[[cs$name]]$arms)) {
    e <- res[[cs$name]]$arms[[nm]]
    tv <- 0.5 * sum(abs(e - res[[cs$name]]$oracle))
    say("   %-55s TV(emp, oracle) = %.4f  (max cell %.4f)", nm, tv, max(abs(e - res[[cs$name]]$oracle)))
    say("     %s", paste(sprintf("%s=%.4f", names(e), e), collapse = "  "))
  }
}
saveRDS(res, file.path(outdir, sprintf("eta-bias-%s.rds", label)))
say("")
say("-> written %s", file.path(outdir, sprintf("eta-bias-%s.rds", label)))
