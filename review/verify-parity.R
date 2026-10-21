#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# Old-vs-new parity harness for the performance work on `glm-parallel`.
#
# The two implementations are both called FlexRL, so they are exercised in
# SEPARATE R processes (R_LIBS pointing at the respective library) that dump
# their results to RDS files; a third process compares the dumps WITHOUT
# loading FlexRL at all.
#
# Usage:
#   R_LIBS=<baseline lib> Rscript review/verify-parity.R old  [outdir] [label]
#   R_LIBS=<current lib>  Rscript review/verify-parity.R new  [outdir] [label]
#   Rscript review/verify-parity.R compare [outdir] [oldlabel] [newlabel]
#                                                                     (no FlexRL load)
#
# The variant is detected at runtime ("unif" %in% names(formals(sampleNL))) and
# must match the requested mode. Default outdir: ../FlexRL-bench/parity next to
# the repo.
#
# The optional label names the build behind a dump: run modes write
# dump-<label>.rds (default label = the mode name, i.e. today's behaviour);
# compare mode reads dump-<oldlabel>.rds vs dump-<newlabel>.rds (default
# old/new). This lets e.g. the upstream master build stand in for the baseline
# (`old` mode, label "master") without touching the default invocations.
#
# Sections (see review/PERF-PLAN.md §6):
#   A. harness sanity      — simulate_data()/prepare_data() outputs, pure R-RNG
#                            code unchanged: must be EXACTLY equal old vs new.
#   B. deterministic grid  — candidate-pair SETS (old string pipeline vs
#                            buildCandidatePairs), LLA/LLB/LLL likelihood
#                            values on fixed states, log_lik() vs
#                            log_lik_pairs(): exact to FP (<= 1e-10).
#   C. stationary dist.    — tiny-problem Gibbs frequencies (30k sweeps, thin
#                            3, seed 2024) vs the exact posterior; TV distances
#                            old-vs-exact, new-vs-exact, old-vs-new.
#   D. StEM seed clouds    — >= 8 seeds per problem (synthetic std + three
#                            real-data problems): per-parameter z-scores of
#                            each cloud against the other, Delta Jaccard
#                            old-vs-old (reference) vs old-vs-new.
#
# The same fixed states and problems (shared code in
# tests/testthat/helper-golden.R + review/parity-problems.R) feed
# review/generate-goldens.R, which stores the baseline outputs as permanent
# regression fixtures validated by tests/testthat/test-golden-parity.R.
# ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: verify-parity.R old|new|compare [outdir] [label | newlabel oldlabel]")
mode <- args[1]
stopifnot(mode %in% c("old", "new", "compare"))
a <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", a, value = TRUE)
repo <- if (length(f)) normalizePath(file.path(dirname(sub("^--file=", "", f[1])), "..")) else normalizePath(".")
outdir <- if (length(args) >= 2) args[2] else file.path(repo, "..", "FlexRL-bench", "parity")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# Output label: names the build behind a dump (run modes write dump-<label>.rds,
# default = the mode name; compare reads dump-<oldlabel>.rds vs
# dump-<newlabel>.rds, default old/new). See the usage header above.
label <- if (mode != "compare" && length(args) >= 3) args[3] else mode
oldlabel <- if (mode == "compare" && length(args) >= 3) args[3] else "old"
newlabel <- if (mode == "compare" && length(args) >= 4) args[4] else "new"

# helper-data.R / helper-golden.R are pure definitions: sourceable without
# FlexRL (the compare step uses their state/cloud helpers and the exact
# posterior, but never loads the package).
source(file.path(repo, "tests", "testthat", "helper-data.R"))
source(file.path(repo, "tests", "testthat", "helper-golden.R"))

if (mode != "compare") {
  library(FlexRL)
  golden_assert_variant(mode)
  source(file.path(repo, "review", "parity-problems.R"))
  options(cli.progress_show_after = 3600) # keep the log readable
}

## ---------------------------------------------------------------- helpers

parity_tv <- golden_tv

parity_align_lll <- function(old_df, new_df) {
  # align by (iA, iB) pair; report max |diff| and whether the pair sets match
  ko <- paste(old_df$iA, old_df$iB, sep = "_")
  kn <- paste(new_df$iA, new_df$iB, sep = "_")
  mo <- match(ko, kn)
  list(
    same_pairs = setequal(ko, kn),
    n_pairs = length(ko),
    max_diff = if (!all(!is.na(mo))) NA_real_ else max(abs(old_df$LLL - new_df$LLL[mo]))
  )
}

## ------------------------------------------------------------------- main

if (mode != "compare") {
  dump <- list(mode = mode, label = label,
               variant = if (golden_is_new()) "new" else "old",
               r_version = R.version.string,
               seeds = GOLDEN_SEEDS)

  ## ---- A. harness sanity --------------------------------------------------
  set.seed(20241021)
  std <- std_generate(seed = 11) # internal set.seed(11): stream-independent
  set.seed(20261)                # externally-seeded second problem
  cfg2 <- std_config()
  gen2 <- simulate_data(cfg2, c(3, 4, 5, 6), c(60, 80), 40,
                        p_mistake = list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
                                         V3 = c(0.05, 0.05), V4 = c(0.02, 0.02)),
                        p_missing = list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
                                         V3 = c(0.005, 0.005), V4 = c(0.005, 0.005)),
                        cond_hazard_params = list(V1 = c(), V2 = c(), V3 = c(),
                                                  V4 = log(c(0.7, 0.6, 0.5))),
                        enforce_estimability = TRUE)
  prep2 <- prepare_data(gen2$data1, gen2$data2, "1", "2", cfg2, TRUE, "entity_id")
  dump$sanity <- list(
    std_A = std$prep$encodedA, std_B = std$prep$encodedB,
    std_n_values = std$prep$n_values, std_true_pairs = std$prep$true_pairs,
    gen2_A = prep2$encodedA, gen2_B = prep2$encodedB,
    gen2_n_values = prep2$n_values, gen2_true_pairs = prep2$true_pairs
  )

  ## ---- B. deterministic kernel parity on fixed states ---------------------
  problems <- parity_problems()
  states <- list()
  for (id in GOLDEN_STATE_IDS) {
    pid <- switch(id, std_links = "std", std_nolinks = "std", nltcs_links = "nltcs",
                  shiwreg_links = "shiwreg", shiwfull_links = "shiwfull",
                  tiny_links = NULL)
    prep <- if (is.null(pid)) {
      # tiny problem: hand-built data from helper-data.R
      td <- tiny_data()
      list(encodedA = td$encodedA, encodedB = td$encodedB, n_values = td$n_values,
           PIVs_config = td$PIVs_config, same_mistakes = td$same_mistakes)
    } else {
      problems[[pid]]$prep
    }
    state <- golden_build_state(id, prep)
    art <- golden_capture_artifacts(state)
    states[[id]] <- list(state = state, artifacts = art)
  }
  dump$states <- states

  ## ---- C. tiny stationary distribution ------------------------------------
  set.seed(2024) # same seed as tests/testthat/test-exact-posterior.R
  dump$tiny_freq <- tiny_gibbs_frequencies(n_sweeps = 30000, thin = 3L)

  ## ---- D. StEM seed clouds -------------------------------------------------
  dump$clouds <- lapply(problems, golden_run_cloud)

  saveRDS(dump, file.path(outdir, sprintf("dump-%s.rds", label)))
  cat(sprintf("verify-parity: mode=%s label=%s variant=%s -> %s\n", mode, label, dump$variant,
              file.path(outdir, sprintf("dump-%s.rds", label))))

} else {
  ## ---------------------------------------------------------------- compare
  old <- readRDS(file.path(outdir, sprintf("dump-%s.rds", oldlabel)))
  new <- readRDS(file.path(outdir, sprintf("dump-%s.rds", newlabel)))
  report <- character(0)
  say <- function(...) {
    line <- sprintf(...)
    report <<- c(report, line)
    cat(line, "\n")
  }

  say("FlexRL parity report — A=%s [variant %s] vs B=%s [variant %s]",
      oldlabel, old$variant, newlabel, new$variant)
  say("")

  ## A. sanity
  say("== A. harness sanity (simulate_data + prepare_data; must be EXACTLY equal)")
  sanity_ok <- TRUE
  for (nm in names(old$sanity)) {
    eq <- identical(old$sanity[[nm]], new$sanity[[nm]])
    sanity_ok <- sanity_ok && eq
    say("   %-16s identical: %s", nm, eq)
  }
  say("   -> harness sanity: %s", if (sanity_ok) "PASS" else "FAIL")
  say("")

  ## B. deterministic kernels
  say("== B. deterministic kernel parity (fixed states; tolerance <= 1e-10)")
  b_ok <- TRUE
  for (id in names(old$states)) {
    o <- old$states[[id]]
    n <- new$states[[id]]
    st_o <- o$state; st_n <- n$state
    state_eq <- identical(st_o$truepivsA, st_n$truepivsA) &&
      identical(st_o$truepivsB, st_n$truepivsB) &&
      identical(st_o$linksR, st_n$linksR) &&
      identical(st_o$sumRowD, st_n$sumRowD) &&
      identical(st_o$eta, st_n$eta) && identical(st_o$phi, st_n$phi) &&
      identical(st_o$data, st_n$data)
    sel_o <- golden_sort_pairs(o$artifacts$select)
    sel_n <- golden_sort_pairs(n$artifacts$select)
    pairs_eq <- identical(sel_o, sel_n)
    roworder_eq <- identical(o$artifacts$select, n$artifacts$select)
    ali <- parity_align_lll(o$artifacts$LLL_pairs, n$artifacts$LLL_pairs)
    dLLA <- max(abs(o$artifacts$LLA - n$artifacts$LLA))
    dLLB <- max(abs(o$artifacts$LLB - n$artifacts$LLB))
    lls <- c(old_internal = o$artifacts$ll0_internal,
             old_sparse = o$artifacts$ll0_sparse,
             new_sparse = n$artifacts$ll0_sparse,
             new_pairs = if (is.null(n$artifacts$ll0_pairs)) NA_real_ else n$artifacts$ll0_pairs)
    dll0 <- max(abs(lls - mean(lls, na.rm = TRUE)), na.rm = TRUE)
    ok <- state_eq && pairs_eq && ali$same_pairs && ali$max_diff <= 1e-10 &&
      dLLA <= 1e-10 && dLLB <= 1e-10 && dll0 <= 1e-10
    b_ok <- b_ok && ok
    say("   [%s] %s", if (ok) "PASS" else "FAIL", id)
    say("      state rebuilt identically: %s | candidate-pair SET equal: %s (row order equal: %s, expected FALSE)",
        state_eq, pairs_eq, roworder_eq)
    say("      n pairs: %d | max|dLLL| = %.3g | max|dLLA| = %.3g | max|dLLB| = %.3g | spread of log_lik values = %.3g",
        ali$n_pairs, ali$max_diff, dLLA, dLLB, dll0)
    say("      log_lik values: old-internal=%.10g old-sparse=%.10g new-sparse=%.10g new-pairs=%.10g",
        lls[1], lls[2], lls[3], lls[4])
  }
  say("   -> deterministic kernels: %s", if (b_ok) "PASS" else "FAIL")
  say("")

  ## C. tiny stationary distribution
  say("== C. tiny-problem stationary distribution (30k sweeps, thin 3, seed 2024)")
  # exact posterior: pure R code from helper-data.R, computed WITHOUT FlexRL
  exact <- tiny_exact_posterior()
  cfgs <- tiny_configs()
  keys <- vapply(cfgs, function(cfg) {
    if (length(cfg[[1]]) == 0) return("0")
    paste(apply(cbind(cfg[[1]], cfg[[2]]), 1, paste, collapse = ">"), collapse = ";")
  }, character(1))
  freqo <- old$tiny_freq; freqn <- new$tiny_freq
  emp <- function(f) {
    v <- f[match(keys, names(f))]
    v[is.na(v)] <- 0
    v
  }
  eo <- emp(freqo); en <- emp(freqn)
  tv_o <- 0.5 * sum(abs(eo - exact))
  tv_n <- 0.5 * sum(abs(en - exact))
  tv_on <- parity_tv(freqo, freqn)
  say("   TV(%s, exact)   = %.4f  (max cell %.4f)", oldlabel, tv_o, max(abs(eo - exact)))
  say("   TV(%s, exact)   = %.4f  (max cell %.4f)", newlabel, tv_n, max(abs(en - exact)))
  say("   TV(%s, %s)     = %.4f  (MC noise level: both are 10k-draw estimates)", oldlabel, newlabel, tv_on)
  say("   test tolerances: max cell 0.025, TV 0.05 vs exact; same-implementation pairs differ by MC error only")
  c_ok <- tv_o < 0.05 && tv_n < 0.05 && max(abs(eo - exact)) < 0.025 &&
    max(abs(en - exact)) < 0.025 && tv_on < 0.05
  say("   -> stationary distribution: %s", if (c_ok) "PASS" else "FLAG (inspect)")
  say("")

  ## D. StEM seed clouds
  say("== D. StEM seed-cloud parity (%d seeds: %s)", length(GOLDEN_SEEDS),
      paste(GOLDEN_SEEDS, collapse = ", "))
  d_ok <- TRUE
  for (p in names(old$clouds)) {
    co <- old$clouds[[p]]
    cn <- new$clouds[[p]]
    zs <- golden_zstats(golden_cloud_params(co), golden_cloud_params(cn))
    zmax_no <- max(zs$new_vs_old)
    zmax_on <- max(zs$old_vs_new)
    worst_no <- names(which.max(zs$new_vs_old))
    worst_on <- names(which.max(zs$old_vs_new))
    jo <- c()
    for (a in seq_along(co)) for (b in seq_len(a - 1L)) {
      jo <- c(jo, golden_jaccard(co[[a]]$links, co[[b]]$links))
    }
    jn <- c()
    for (a in co) for (b in cn) jn <- c(jn, golden_jaccard(a$links, b$links))
    nlo <- vapply(co, function(x) x$n_links, numeric(1))
    nln <- vapply(cn, function(x) x$n_links, numeric(1))
    ok <- zmax_no <= 4 && zmax_on <= 4
    d_ok <- d_ok && ok
    say("   [%s] %s", if (ok) "PASS" else "FLAG", p)
    say("      n links (Delta x>0.5): %s %d-%d, %s %d-%d",
        oldlabel, min(nlo), max(nlo), newlabel, min(nln), max(nln))
    say("      max|z| %s-vs-%s-cloud = %.2f (%s) | %s-vs-%s-cloud = %.2f (%s)",
        newlabel, oldlabel, zmax_no, worst_no, oldlabel, newlabel, zmax_on, worst_on)
    say("      Delta Jaccard %s-vs-%s: mean %.3f sd %.3f min %.3f | %s-vs-%s: mean %.3f sd %.3f min %.3f",
        oldlabel, oldlabel, mean(jo), stats::sd(jo), min(jo),
        oldlabel, newlabel, mean(jn), stats::sd(jn), min(jn))
  }
  say("   -> seed clouds: %s", if (d_ok) "PASS" else "FLAG (inspect)")
  say("")

  overall <- sanity_ok && b_ok && c_ok && d_ok
  say("OVERALL: %s", if (overall) "PASS" else "ATTENTION NEEDED")
  logfile <- if (oldlabel == "old" && newlabel == "new") "parity-compare.log" else
    sprintf("parity-compare-%s-vs-%s.log", oldlabel, newlabel)
  writeLines(report, file.path(outdir, logfile))
}
