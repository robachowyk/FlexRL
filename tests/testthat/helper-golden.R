# Shared helpers for the golden-parity fixtures: regression tests that pin the
# current implementation to reference outputs produced by the PRE-optimisation
# baseline (branch `glm-notes` @ 4b09088, built in FlexRL-bench/rlib-old).
#
# This file is deliberately version-agnostic: it is sourced both by
# `review/verify-parity.R` / `review/generate-goldens.R` (which run under
# whichever FlexRL build `R_LIBS` points at) and by the testthat suite (which
# runs under the current implementation). Where the two implementations differ
# (candidate-pair construction, likelihood bookkeeping), the code branches on
# the variant detected at runtime.
#
# Statistical ground rules (review/PERF-PLAN.md §5-§6): the optimisation changed
# HOW uniforms are consumed (pre-drawn blocks, inverse-CDF, factorised
# sampleL), so two implementations given the same seed produce different
# chains BY DESIGN. What must hold is:
#   * deterministic kernels (likelihood values, candidate-pair SETS) agree
#     exactly (up to floating-point, <= 1e-10);
#   * stochastic outputs agree in distribution: the seed-cloud statistics of
#     one implementation must look like just another seed of the other.

## ---------------------------------------------------------------- variant

#' Is the loaded FlexRL the performance-optimised (glm-parallel) variant?
#'
#' The optimised kernels consume pre-drawn uniforms (the `unif` argument of
#' `sampleNL()`); the baseline kernels draw their own uniforms through
#' `Rcpp::sample()`.
golden_is_new <- function() {
  "unif" %in% names(formals(FlexRL:::sampleNL))
}

#' Stop unless the loaded variant matches `mode` ("old" / "new")
golden_assert_variant <- function(mode) {
  new <- golden_is_new()
  if (mode == "new" && !new) {
    stop("Expected the optimised (glm-parallel) FlexRL but found the baseline variant.", call. = FALSE)
  }
  if (mode == "old" && new) {
    stop("Expected the baseline (glm-notes) FlexRL but found the optimised variant.", call. = FALSE)
  }
  invisible(NULL)
}

## ------------------------------------------------------------------ seeds

## Seeds for the StEM seed-cloud runs. The SAME set is used for every problem
## and both implementations: the cloud (not a single run) is the unit of
## comparison, because the RNG streams differ between implementations by design.
GOLDEN_SEEDS <- c(11, 101, 202, 303, 404, 505, 606, 707, 808, 909, 1010, 1111)

## Seeds for the deterministic fixed states of the kernel-level checks. The
## state builders use only base-R randomness (set.seed + sample.int), so both
## implementations rebuild bit-identical states from these seeds.
GOLDEN_STATE_SEED <- list(
  std_links = 501, std_nolinks = 502, nltcs_links = 503,
  shiwreg_links = 504, shiwfull_links = 505, tiny_links = 506
)
GOLDEN_STATE_NLINKS <- list(
  std_links = 40L, std_nolinks = 0L, nltcs_links = 100L,
  shiwreg_links = 8L, shiwfull_links = 150L, tiny_links = NULL
)
GOLDEN_STATE_GAMMA <- list(
  std_links = 0.4, std_nolinks = 0.4, nltcs_links = 0.5,
  shiwreg_links = 0.6, shiwfull_links = 0.5, tiny_links = 0.6
)
GOLDEN_STATE_IDS <- c("std_links", "std_nolinks", "nltcs_links", "shiwreg_links",
                      "shiwfull_links", "tiny_links")

## ------------------------------------------- deterministic parameter presets
## (no RNG involved: identical under both implementations by construction)

golden_eta <- function(nv) {
  v <- seq_len(nv)^(-0.3)
  v / sum(v)
}

# (agreement A, agreement B, missing A, missing B)
golden_phi <- function() c(0.90, 0.86, 0.08, 0.10)

golden_alpha <- function(prep) {
  PIVs <- names(prep$PIVs_config)
  stable <- vapply(prep$PIVs_config, function(x) x$dynamics != "structured", logical(1))
  lapply(seq_along(PIVs), function(idx) {
    if (stable[idx]) return(-Inf)
    cc <- prep$PIVs_config[[idx]]$cond_hazard_cov
    n <- ncol(prep$encodedA[, cc$covA, drop = FALSE]) +
      ncol(prep$encodedB[, cc$covB, drop = FALSE]) + 1L
    rep(-0.4, n)
  })
}

#' Rebuild a survival-model object from its small serialisable spec
golden_model <- function(spec = "exponential") {
  switch(spec,
    exponential = FlexRL::survival_model("exponential"),
    custom06 = FlexRL::survival_model(
      "custom",
      S     = function(X, alpha, times) rep(0.6, length(times)),
      n_par = function(n_cov) 1L,
      init  = function(n_cov) 0
    ),
    stop("unknown model spec: ", spec, call. = FALSE)
  )
}

## --------------------------------------------------------- candidate pairs

#' The baseline's candidate-pair pipeline, verbatim from the pre-optimisation
#' `simulateD()`: pasteIntoPattern -> unique -> factor -> indexPatterns ->
#' pairPatterns. Still present (and used by naive_linkage()) in both versions.
golden_candidate_pairs_old <- function(truepivsA, truepivsB, PIVs_stable) {
  UA <- FlexRL:::pasteIntoPattern(as.matrix(truepivsA[, PIVs_stable]))
  UB <- FlexRL:::pasteIntoPattern(as.matrix(truepivsB[, PIVs_stable]))
  valuesU <- unique(c(UA, UB))
  UA <- as.numeric(factor(UA, levels = valuesU))
  UB <- as.numeric(factor(UB, levels = valuesU))
  tmpA <- FlexRL:::indexPatterns(UA, length(valuesU))
  tmpB <- FlexRL:::indexPatterns(UB, length(valuesU))
  FlexRL:::pairPatterns(tmpA, tmpB, length(tmpA))
}

#' Candidate pairs under the loaded implementation. NOTE: the ROW ORDER of the
#' two implementations differs by design (baseline: pattern-bucket order of
#' first appearance; optimised: lexicographic by pattern tuple). Only the SET
#' of pairs is comparable.
golden_candidate_pairs <- function(truepivsA, truepivsB, PIVs_stable) {
  if (golden_is_new()) {
    FlexRL:::buildCandidatePairs(
      as.matrix(truepivsA[, PIVs_stable, drop = FALSE]),
      as.matrix(truepivsB[, PIVs_stable, drop = FALSE])
    )
  } else {
    golden_candidate_pairs_old(truepivsA, truepivsB, PIVs_stable)
  }
}

golden_sort_pairs <- function(select) {
  select[order(select[, 1], select[, 2]), , drop = FALSE]
}

golden_pair_key <- function(select) paste(select[, 1], select[, 2], sep = "_")

## ------------------------------------------------------- tiny frequencies

#' Total variation distance between two named frequency vectors (unobserved
#' categories count as 0)
golden_tv <- function(a, b) {
  nms <- union(names(a), names(b))
  a <- a[nms]
  b <- b[nms]
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  0.5 * sum(abs(a - b))
}

#' Align an empirical tiny-problem frequency vector (keys as produced by
#' tiny_gibbs_frequencies) with the 13 enumerated configurations; missing
#' keys count as 0.
golden_tiny_empirical <- function(freq) {
  cfgs <- tiny_configs()
  keys <- vapply(cfgs, function(cfg) {
    if (length(cfg[[1]]) == 0) return("0")
    paste(apply(cbind(cfg[[1]], cfg[[2]]), 1, paste, collapse = ">"), collapse = ";")
  }, character(1))
  v <- freq[match(keys, names(freq))]
  v[is.na(v)] <- 0
  v
}

## ------------------------------------------------------------ fixed states

#' Build a deterministic fixed state (data + latent true values + links +
#' parameters) for the kernel-level likelihood checks.
#'
#' `prep` must be a reduced prepare_data() output (see parity_reduce_prep()).
#' All randomness is base-R (`set.seed` + `sample.int`) drawn in a fixed order,
#' so both implementations rebuild bit-identical states.
#'
#' The state's links are drawn from the candidate-pair set itself (sorted into
#' the canonical (i, j) order first, so the sampled links do not depend on the
#' implementations' differing candidate-pair row orders).
golden_build_state <- function(id, prep) {
  PIVs <- names(prep$PIVs_config)
  PIVs_stable <- vapply(prep$PIVs_config, function(x) x$dynamics != "structured", logical(1))
  nA <- nrow(prep$encodedA)
  nB <- nrow(prep$encodedB)
  nvs <- prep$n_values
  if (is.null(names(nvs))) {
    nvs <- stats::setNames(rep(nvs, length.out = length(PIVs)), PIVs)
  } else {
    nvs <- nvs[PIVs]
  }

  set.seed(GOLDEN_STATE_SEED[[id]])
  # The latent true values start from the REGISTERED values (0 = missing is
  # imputed by resampling from the pooled observed values), then ~10% of the
  # records are perturbed to a DIFFERENT value per PIV. This keeps the state
  # realistic (true matches share patterns, so candidate sets have the same
  # order of magnitude as the naive exact-match linkage) while still exercising
  # every branch of the likelihood (agreement, typo, missing, survival).
  # Uniform draws over 1..n_values would instead put almost every record in its
  # own pattern bucket for high-cardinality real PIVs.
  draw_tp <- function(n, df) {
    m <- as.matrix(df[, PIVs, drop = FALSE])
    storage.mode(m) <- "integer"
    for (p in seq_along(PIVs)) {
      nv <- as.integer(nvs[p])
      pool <- c(prep$encodedA[[PIVs[p]]], prep$encodedB[[PIVs[p]]])
      pool <- pool[pool > 0]
      z <- which(m[, p] == 0L)
      if (length(z) > 0L) {
        m[z, p] <- if (length(pool) > 0L) {
          pool[sample.int(length(pool), length(z), replace = TRUE)]
        } else {
          sample.int(nv, length(z), replace = TRUE)
        }
      }
      r <- sample.int(n, max(1L, as.integer(round(0.10 * n))))
      alt <- sample.int(nv - 1L, length(r), replace = TRUE)
      alt <- alt + (alt >= m[r, p])
      m[r, p] <- alt
    }
    colnames(m) <- PIVs
    m
  }
  if (id == "tiny_links") {
    # hand-picked tiny state (1 structured PIV, 2 values, custom survival model)
    truepivsA <- matrix(c(2L, 1L), ncol = 1L, dimnames = list(NULL, PIVs))
    truepivsB <- matrix(c(1L, 2L, 1L), ncol = 1L, dimnames = list(NULL, PIVs))
  } else {
    truepivsA <- draw_tp(nA, prep$encodedA)
    truepivsB <- draw_tp(nB, prep$encodedB)
  }

  select <- golden_sort_pairs(golden_candidate_pairs(truepivsA, truepivsB, PIVs_stable))
  if (id == "tiny_links") {
    linksR <- matrix(c(1L, 2L, 1L, 3L), ncol = 2L)
  } else {
    k <- GOLDEN_STATE_NLINKS[[id]]
    stopifnot(k >= 0L)
    # Greedy one-to-one pick of k candidate pairs: random order over the A
    # records that appear in the candidate set, random still-free partner for
    # each. Deterministic given the state seed.
    picked <- character(0)
    usedB <- integer(0)
    uA <- unique(select[, 1])
    for (i in uA[sample.int(length(uA))]) {
      if (length(picked) >= k) break
      js <- setdiff(select[select[, 1] == i, 2], usedB)
      if (length(js) > 0) {
        j <- js[sample.int(length(js), 1L)]
        usedB <- c(usedB, j)
        picked <- c(picked, paste(i, j, sep = "_"))
      }
    }
    stopifnot(length(picked) == k)
    linksR <- golden_sort_pairs(matrix(as.integer(unlist(strsplit(picked, "_"))),
                                       ncol = 2L, byrow = TRUE))
    storage.mode(linksR) <- "integer"
  }
  sumRowD <- seq_len(nA) %in% linksR[, 1]
  sumColD <- seq_len(nB) %in% linksR[, 2]

  list(
    id = id, data = prep,
    model_spec = if (id == "tiny_links") "custom06" else "exponential",
    linksR = linksR, sumRowD = sumRowD, sumColD = sumColD,
    truepivsA = truepivsA, truepivsB = truepivsB,
    gamma = GOLDEN_STATE_GAMMA[[id]],
    eta = lapply(nvs, golden_eta),
    phi = lapply(nvs, function(nv) golden_phi()),
    alpha = if (id == "tiny_links") list(0) else golden_alpha(prep)
  )
}

## --------------------------------------------- likelihood-artifact capture

#' Run the loaded implementation's `simulateD()` on a fixed state and capture
#' the deterministic likelihood artifacts.
#'
#' `simulateD()` computes LLA / LLB / LLL internally and hands them to the C++
#' `sampleD()`. To capture them without reimplementing anything, `sampleD` is
#' shadowed by an R shim in a child environment of the package namespace: the
#' shim records its arguments and returns a self-consistent "no links" result
#' so `simulateD()`'s own log-likelihood sanity check still passes. The
#' returned Dsample is NOT a valid state update and must not be used as one.
#'
#' Returns:
#'   select      candidate-pair matrix in the implementation's own row order
#'   LLA, LLB    per-record non-linked log-likelihood contributions
#'   LLL_pairs   data.frame(iA, iB, LLL): per-candidate-pair linked contribution
#'   gamma       the pLink vector passed to sampleD
#'   ll0_sparse  log_lik() on the state (sparse LLL route; both implementations)
#'   ll0_internal (baseline only) the value simulateD computed internally
#'   ll0_pairs   (optimised only) log_lik_pairs() on the state (vector route)
golden_capture_artifacts <- function(state) {
  FlexRL:::initDeltaMap()
  simD <- FlexRL:::simulateD
  captured <- new.env(parent = emptyenv())
  # Variant-specific signature: simulateD() passes every argument positionally,
  # so only a shim whose formals match the caller's signature binds them.
  #   baseline : sampleD(S, LLA, LLB, LLL, gamma, loglik, nlinkrec, sumRowD, sumColD)
  #   optimised: sampleD(S, LLA, LLB, LLL, gamma, nlinkrec, sumRowD, sumColD)
  shim <- if (golden_is_new()) {
    function(S, LLA, LLB, LLL, gamma, nlinkrec, sumRowD, sumColD) {
      captured$S <- S
      captured$LLA <- LLA
      captured$LLB <- LLB
      captured$LLL <- LLL
      captured$gamma <- gamma
      captured$ll0_internal <- NULL
      # Self-consistent loglik for the returned empty-links state: the ll_check
      # inside simulateD() then compares equal values and does not stop.
      ll <- FlexRL:::log_lik_pairs(LLL, LLA, LLB, integer(0), 0L, sumRowD, sumColD, gamma)
      list(links = matrix(0L, 0L, 2L), linksQ = integer(0),
           sumRowD = sumRowD, sumColD = sumColD, loglik = ll, nlinkrec = 0L)
    }
  } else {
    function(S, LLA, LLB, LLL, gamma, loglik, nlinkrec, sumRowD, sumColD) {
      captured$S <- S
      captured$LLA <- LLA
      captured$LLB <- LLB
      captured$LLL <- LLL
      captured$gamma <- gamma
      captured$ll0_internal <- loglik
      ll <- FlexRL:::log_lik(LLL, LLA, LLB, matrix(0L, 0L, 2L), sumRowD, sumColD, gamma)
      list(links = matrix(0L, 0L, 2L), linksQ = integer(0),
           sumRowD = sumRowD, sumColD = sumColD, loglik = ll, nlinkrec = 0L)
    }
  }
  myenv <- new.env(parent = environment(simD))
  assign("sampleD", shim, envir = myenv)
  environment(simD) <- myenv
  Dsample <- simD(
    data = state$data, linksR = state$linksR, sumRowD = state$sumRowD,
    sumColD = state$sumColD, truepivsA = state$truepivsA, truepivsB = state$truepivsB,
    gamma = state$gamma, eta = state$eta, alpha = state$alpha, phi = state$phi,
    model_dynamics = golden_model(state$model_spec)
  )
  rm(Dsample) # shim output is not a valid state update

  select <- captured$S
  LLA <- as.numeric(captured$LLA)
  LLB <- as.numeric(captured$LLB)
  LLL <- as.numeric(captured$LLL)
  gamma <- as.numeric(captured$gamma)
  nA <- nrow(state$data$encodedA)
  nB <- nrow(state$data$encodedB)
  LLL_sparse <- Matrix::sparseMatrix(i = select[, 1], j = select[, 2], x = LLL,
                                     dims = c(nA, nB))
  out <- list(
    id = state$id, select = select, LLA = LLA, LLB = LLB, gamma = gamma,
    LLL_pairs = data.frame(iA = as.integer(select[, 1]),
                           iB = as.integer(select[, 2]),
                           LLL = LLL),
    ll0_sparse = FlexRL::log_lik(LLL_sparse, LLA, LLB, state$linksR,
                                 state$sumRowD, state$sumColD, gamma),
    ll0_internal = captured$ll0_internal
  )
  if (golden_is_new()) {
    # linksQ of the fixed state: q-rows of `select` matching the state's links.
    # A state link that is not a candidate pair has no LLL contribution in
    # either route (sparse LLL is 0 there; linksQ cannot name it), so NA -> drop
    # keeps the two routes numerically identical.
    q <- match(golden_pair_key(state$linksR), golden_pair_key(select))
    out$ll0_pairs <- FlexRL:::log_lik_pairs(LLL, LLA, LLB, q[!is.na(q)],
                                            nrow(state$linksR),
                                            state$sumRowD, state$sumColD, gamma)
  }
  out
}

## ------------------------------------------------------------- StEM clouds

#' Run StEM() for every seed in GOLDEN_SEEDS on a prepared problem and return
#' the per-seed summaries (see golden_summarize_stem()).
golden_run_cloud <- function(problem) {
  model <- golden_model(problem$model_spec)
  out <- list()
  for (s in GOLDEN_SEEDS) {
    set.seed(s)
    FlexRL:::initDeltaMap()
    fit <- FlexRL::StEM(
      data = problem$prep,
      StEM_iter = problem$StEM_iter, StEM_burnin = problem$StEM_burnin,
      gibbs_iter = problem$gibbs_iter, gibbs_burnin = problem$gibbs_burnin,
      n_post_sample = problem$n_post_sample,
      model_dynamics = model
    )
    out[[as.character(s)]] <- golden_summarize_stem(fit, problem$prep)
  }
  out
}

#' Final-iteration summary of one StEM() fit: gamma, eta, phi, alpha (unstable
#' PIVs only) and the Delta link set at threshold x > 0.5.
golden_summarize_stem <- function(fit, prep) {
  it <- NROW(fit$gamma)
  PIVs <- names(prep$PIVs_config)
  unstable <- !vapply(prep$PIVs_config, function(x) x$dynamics != "structured", logical(1))
  keep <- fit$Delta$x > 0.5
  links <- if (any(keep)) {
    unname(cbind(as.integer(fit$Delta$i[keep]), as.integer(fit$Delta$j[keep])))
  } else {
    matrix(integer(0), 0L, 2L)
  }
  links <- links[order(links[, 1], links[, 2]), , drop = FALSE]
  list(
    gamma = as.numeric(fit$gamma[it, ]),
    eta = lapply(fit$eta, function(m) as.numeric(m[it, ])),
    phi = lapply(fit$phi, function(m) as.numeric(m[it, ])),
    alpha = stats::setNames(
      lapply(which(unstable), function(k) as.numeric(fit$alpha[[k]][it, ])),
      PIVs[unstable]
    ),
    n_links = sum(keep),
    links = links
  )
}

golden_jaccard <- function(a, b) {
  ka <- if (nrow(a) == 0L) character(0) else paste(a[, 1], a[, 2], sep = "_")
  kb <- if (nrow(b) == 0L) character(0) else paste(b[, 1], b[, 2], sep = "_")
  if (length(ka) == 0L && length(kb) == 0L) return(1)
  length(intersect(ka, kb)) / length(union(ka, kb))
}

#' Collect the scalar parameter components of a seed-cloud into named matrices
#' (rows = seeds) for z-score comparisons.
golden_cloud_params <- function(cloud) {
  comps <- list()
  comps$gamma <- do.call(rbind, lapply(cloud, function(s) matrix(s$gamma, 1L)))
  one <- cloud[[1]]
  for (k in seq_along(one$eta)) {
    comps[[paste0("eta[[", k, "]]")]] <-
      do.call(rbind, lapply(cloud, function(s) matrix(s$eta[[k]], 1L)))
  }
  for (k in seq_along(one$phi)) {
    comps[[paste0("phi[[", k, "]]")]] <-
      do.call(rbind, lapply(cloud, function(s) matrix(s$phi[[k]], 1L)))
  }
  if (length(one$alpha) > 0) {
    for (k in seq_along(one$alpha)) {
      comps[[paste0("alpha[[", names(one$alpha)[k], "]]")]] <-
        do.call(rbind, lapply(cloud, function(s) matrix(s$alpha[[k]], 1L)))
    }
  }
  comps
}

#' Max |z| of one cloud's seeds against the other cloud's per-component
#' mean/sd, and vice versa.
#'
#' The scale for each column is the LARGER of the two clouds' sds. The old
#' cloud alone is not a valid scale for near-degenerate columns: e.g. a
#' zero-frequency eta category can be exactly 0 in every seed of one cloud and
#' occasionally 1/(draws) in a seed of the other — pure MC noise on a rare
#' category, which a 0-sd reference would report as |z| ~ 1e5. With the larger
#' sd as scale such columns stay at MC level, while a genuine level shift
#' between two constant clouds (both sds ~ 0) still yields a huge z through the
#' 1e-7 absolute floor.
golden_zstats <- function(old_params, new_params) {
  stopifnot(identical(names(old_params), names(new_params)))
  zscore <- function(target, reference, scale) {
    # per-COLUMN standardisation: sweep() both steps (plain `/` on a matrix
    # would recycle the sd vector column-major, i.e. mismatched)
    d <- sweep(target, 2, colMeans(reference), "-")
    sweep(d, 2, pmax(scale, 1e-7), "/")
  }
  z_n_o <- vapply(names(old_params), function(nm) {
    o <- old_params[[nm]]
    n <- new_params[[nm]]
    s <- pmax(apply(o, 2, sd), apply(n, 2, sd))
    max(abs(zscore(n, o, s)))
  }, numeric(1))
  z_o_n <- vapply(names(old_params), function(nm) {
    o <- old_params[[nm]]
    n <- new_params[[nm]]
    s <- pmax(apply(o, 2, sd), apply(n, 2, sd))
    max(abs(zscore(o, n, s)))
  }, numeric(1))
  list(new_vs_old = z_n_o, old_vs_new = z_o_n)
}
