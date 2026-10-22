# Problem builders for the parity harness (review/verify-parity.R) and the
# golden-fixture generation (review/generate-goldens.R). NOT sourced by the
# testthat suite: the tests consume the stored goldens and never read the
# vignette CSVs.
#
# The problems mirror the vignettes' own FlexRL workflows:
#   * std      : the synthetic 4-PIV problem of tests/testthat/helper-data.R
#                (stable, stable, flexible, structured with covariates).
#   * shiwreg  : vignettes/ex-reg-SHIW-16.csv / ex-reg-SHIW-20.csv with the
#                regional-SHIW PIVs_config of the block vignette (4 stable
#                PIVs: ANASCI, SESSO, STUDIO, STACIV).
#   * nltcs    : vignettes/ex-reg-NLTCS-82.csv / ex-reg-NLTCS-94.csv with the
#                block vignette's config (SEX, DOB_MONTH, DOB_YEAR stable;
#                STATE flexible) — exercises the survival/alpha machinery on
#                real data.
#   * shiwfull : vignettes/SHIW2016.csv / SHIW2020.csv (deduplicated by ID as
#                the vignette does), then subsampled to the FIRST N rows
#                (deterministic; N = 2500, chosen so a full seed-cloud also
#                completes under the slow baseline build) with the old block
#                vignette's full-SHIW config (ANASCI, SESSO, STUDIO, STACIV,
#                IREG, all stable).
#
# StEM settings are scaled down (a full seed-cloud must finish in minutes even
# under the baseline build): 6 StEM iterations (3 burn-in), 8 Gibbs iterations
# (4 burn-in), 20 posterior samples for the real-data problems; the synthetic
# problem uses the test suite's own std_stem_args().

parity_repo_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) > 0) {
    normalizePath(file.path(dirname(sub("^--file=", "", f[1])), ".."))
  } else {
    normalizePath(".")
  }
}

parity_vignettes <- function() file.path(parity_repo_root(), "vignettes")

#' Reduce a prepare_data() output to the columns StEM()/simulateD() actually
#' touch (PIVs, `date` when any PIV is unstable, and the conditional-hazard
#' covariates — cov1 columns in file A, cov2 columns in file B). Keeps the
#' golden fixtures small and self-describing; has no effect on any computation.
parity_reduce_prep <- function(prep) {
  PIVs <- names(prep$PIVs_config)
  unstable <- !vapply(prep$PIVs_config, function(x) x$dynamics != "structured", logical(1))
  needA <- PIVs
  needB <- PIVs
  if (any(unstable)) {
    needA <- c(needA, "date")
    needB <- c(needB, "date")
  }
  for (idx in which(unstable)) {
    cc <- prep$PIVs_config[[idx]]$cond_hazard_cov
    needA <- c(needA, cc$covA)
    needB <- c(needB, cc$covB)
  }
  needA <- unique(needA)
  needB <- unique(needB)
  stopifnot(all(needA %in% colnames(prep$encodedA)),
            all(needB %in% colnames(prep$encodedB)))
  list(
    encodedA = prep$encodedA[, needA, drop = FALSE],
    encodedB = prep$encodedB[, needB, drop = FALSE],
    n_values = prep$n_values,
    PIVs_config = prep$PIVs_config,
    same_mistakes = prep$same_mistakes
  )
}

parity_problems <- function(N_shiwfull = 2500) {
  v <- parity_vignettes()

  ## --- synthetic std problem (helper-data.R; fixed seed inside) ------------
  std <- std_generate(seed = 11)

  ## --- regional SHIW (ex-reg-SHIW-16/20): block vignette config ------------
  df16 <- read.csv(file.path(v, "ex-reg-SHIW-16.csv"), row.names = 1)
  df20 <- read.csv(file.path(v, "ex-reg-SHIW-20.csv"), row.names = 1)
  cfg_shiw <- list(
    ANASCI = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
    SESSO  = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
    STUDIO = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
    STACIV = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA))
  )
  prep_shiwreg <- prepare_data(df16, df20, "2016", "2020", cfg_shiw,
                               same_mistakes = TRUE, uniq_id = "ID")

  ## --- regional NLTCS (ex-reg-NLTCS-82/94): block vignette config ----------
  df82 <- read.csv(file.path(v, "ex-reg-NLTCS-82.csv"), row.names = 1)
  df94 <- read.csv(file.path(v, "ex-reg-NLTCS-94.csv"), row.names = 1)
  cfg_nltcs <- list(
    SEX       = list(dynamics = "stable",   bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
    DOB_MONTH = list(dynamics = "stable",   bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
    DOB_YEAR  = list(dynamics = "stable",   bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
    STATE     = list(dynamics = "flexible", bound_mistakes = c(NA, NA),     fix_mistakes = c(NA, NA))
  )
  prep_nltcs <- prepare_data(df82, df94, "1982", "1994", cfg_nltcs,
                             same_mistakes = TRUE, uniq_id = "ID")

  ## --- full SHIW 2016/2020 (subsampled): old block vignette config ---------
  df16f <- read.csv(file.path(v, "SHIW2016.csv"), row.names = 1)
  df20f <- read.csv(file.path(v, "SHIW2020.csv"), row.names = 1)
  df16f <- df16f[!df16f$ID %in% df16f$ID[duplicated(df16f$ID)], ]
  df20f <- df20f[!df20f$ID %in% df20f$ID[duplicated(df20f$ID)], ]
  # deterministic subsample: first N rows of each deduplicated file
  df16f <- df16f[seq_len(min(N_shiwfull, nrow(df16f))), ]
  df20f <- df20f[seq_len(min(N_shiwfull, nrow(df20f))), ]
  cfg_shiwfull <- c(cfg_shiw, list(
    IREG = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA))
  ))
  prep_shiwfull <- prepare_data(df16f, df20f, "2016", "2020", cfg_shiwfull,
                                same_mistakes = TRUE, uniq_id = "ID")

  list(
    std = list(
      prep = parity_reduce_prep(std$prep), model_spec = "exponential",
      StEM_iter = 10, StEM_burnin = 5, gibbs_iter = 8, gibbs_burnin = 4,
      n_post_sample = 50
    ),
    shiwreg = list(
      prep = parity_reduce_prep(prep_shiwreg), model_spec = "exponential",
      StEM_iter = 6, StEM_burnin = 3, gibbs_iter = 8, gibbs_burnin = 4,
      n_post_sample = 20
    ),
    nltcs = list(
      prep = parity_reduce_prep(prep_nltcs), model_spec = "exponential",
      StEM_iter = 6, StEM_burnin = 3, gibbs_iter = 8, gibbs_burnin = 4,
      n_post_sample = 20
    ),
    shiwfull = list(
      prep = parity_reduce_prep(prep_shiwfull), model_spec = "exponential",
      StEM_iter = 6, StEM_burnin = 3, gibbs_iter = 8, gibbs_burnin = 4,
      n_post_sample = 20
    )
  )
}
