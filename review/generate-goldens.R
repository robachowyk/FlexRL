#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# Generate the golden-parity fixtures from the BASELINE implementation.
#
# NOT run by testthat. The committed fixtures under tests/testthat/golden/ were
# generated ONCE from the pre-optimisation baseline (branch `glm-notes`
# @ 4b09088, the commit `glm-parallel` was branched from, i.e. including all
# review-phase statistical fixes) built in a separate library:
#
#   R_LIBS=/mnt/halveGig/kaya/Broncode/FlexRL-bench/rlib-old \
#     Rscript review/generate-goldens.R
#
# Regenerating requires the baseline build (worktree at
# /mnt/halveGig/kaya/Broncode/FlexRL-baseline, installed with
#   R CMD INSTALL --library=/mnt/halveGig/kaya/Broncode/FlexRL-bench/rlib-old \
#     /mnt/halveGig/kaya/Broncode/FlexRL-baseline
# ). The script REFUSES to run against the optimised variant (detected via
# "unif" %in% names(formals(FlexRL:::sampleNL))).
#
# Fixtures written (full precision, gzip-compressed saveRDS, nothing rounded):
#   state_<id>.rds        the fixed state (reduced prepared data, latent true
#                         values, links, parameters, model spec) — plain R
#                         objects, fully reconstructible by the new package
#   pairs_<id>.rds        candidate-pair SET from the baseline string pipeline
#                         (rows sorted canonically; row order is NOT part of
#                         the contract)
#   likelihood_<id>.rds   LLA, LLB, per-pair LLL aligned by (iA, iB), and the
#                         baseline log_lik values (internal + sparse route)
#   stem_cloud_<pid>.rds  per-seed StEM summaries (final gamma/eta/phi/alpha,
#                         Delta x > 0.5 link set) for the 12 seeds
#   tiny_frequencies.rds  tiny-problem Gibbs frequencies (30k sweeps, thin 3,
#                         seed 2024 — same seed as test-exact-posterior.R)
#   meta.rds              provenance (variant, baseline ref, R version, seeds)
#
# State/cloud/tiny definitions are SHARED with the validation side
# (tests/testthat/helper-golden.R sourced by test-golden-parity.R), so the
# tests rebuild the exact same states and rerun the exact same seeds.
# ----------------------------------------------------------------------------

suppressMessages(library(FlexRL))

a <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", a, value = TRUE)
repo <- if (length(f)) normalizePath(file.path(dirname(sub("^--file=", "", f[1])), "..")) else normalizePath(".")
source(file.path(repo, "tests", "testthat", "helper-data.R"))
source(file.path(repo, "tests", "testthat", "helper-golden.R"))
source(file.path(repo, "review", "parity-problems.R"))
golden_assert_variant("old")

outdir <- file.path(repo, "tests", "testthat", "golden")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
options(cli.progress_show_after = 3600)

problems <- parity_problems()

state_prep <- function(id, problems) {
  pid <- switch(id,
                std_links = "std", std_nolinks = "std", nltcs_links = "nltcs",
                shiwreg_links = "shiwreg", shiwfull_links = "shiwfull",
                tiny_links = NULL)
  if (is.null(pid)) {
    td <- tiny_data()
    list(encodedA = td$encodedA, encodedB = td$encodedB, n_values = td$n_values,
         PIVs_config = td$PIVs_config, same_mistakes = td$same_mistakes)
  } else {
    problems[[pid]]$prep
  }
}

## ---- fixed states: pairs + likelihood artifacts ---------------------------
for (id in GOLDEN_STATE_IDS) {
  state <- golden_build_state(id, state_prep(id, problems))
  art <- golden_capture_artifacts(state)
  saveRDS(state, file.path(outdir, sprintf("state_%s.rds", id)))
  saveRDS(
    list(select_sorted = golden_sort_pairs(art$select), n_pairs = nrow(art$select)),
    file.path(outdir, sprintf("pairs_%s.rds", id))
  )
  saveRDS(
    list(
      LLA = art$LLA, LLB = art$LLB, gamma = art$gamma,
      LLL_pairs = art$LLL_pairs,
      ll0_sparse = art$ll0_sparse,
      ll0_internal = art$ll0_internal   # simulateD()'s own value (baseline only)
    ),
    file.path(outdir, sprintf("likelihood_%s.rds", id))
  )
  cat("state:", id, sprintf("(%d pairs)\n", nrow(art$select)))
}

## ---- StEM seed clouds ------------------------------------------------------
for (pid in names(problems)) {
  pr <- problems[[pid]]
  cloud <- golden_run_cloud(pr)
  saveRDS(
    list(
      problem = list(model_spec = pr$model_spec, StEM_iter = pr$StEM_iter,
                     StEM_burnin = pr$StEM_burnin, gibbs_iter = pr$gibbs_iter,
                     gibbs_burnin = pr$gibbs_burnin, n_post_sample = pr$n_post_sample,
                     seeds = GOLDEN_SEEDS),
      cloud = cloud
    ),
    file.path(outdir, sprintf("stem_cloud_%s.rds", pid))
  )
  cat("cloud:", pid, sprintf("(%d seeds)\n", length(cloud)))
}

## ---- tiny-problem stationary frequencies -----------------------------------
set.seed(2024) # same seed as tests/testthat/test-exact-posterior.R
freq <- tiny_gibbs_frequencies(n_sweeps = 30000, thin = 3L)
saveRDS(freq, file.path(outdir, "tiny_frequencies.rds"))
cat("tiny frequencies:", length(freq), "observed configurations\n")

## ---- provenance -------------------------------------------------------------
saveRDS(
  list(
    generator = "review/generate-goldens.R",
    baseline_ref = "glm-notes @ 4b09088",
    baseline_description = "pre-optimisation implementation incl. review-phase statistical fixes",
    variant = "old",
    r_version = R.version.string,
    date = as.character(Sys.Date()),
    seeds = GOLDEN_SEEDS,
    state_seeds = GOLDEN_STATE_SEED
  ),
  file.path(outdir, "meta.rds")
)
cat("goldens written to", outdir, "\n")
