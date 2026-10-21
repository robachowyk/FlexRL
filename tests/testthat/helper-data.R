# Shared helpers for the FlexRL test suite.
#
# Two families of fixtures:
#  1. `std_*`  : the vignette-style 4-PIV problem (moderate size) for
#                pipeline/integration tests.
#  2. `tiny_*` : hand-built, fully known-truth micro problems for exact
#                mathematical checks (Gibbs conditionals, enumeration).
#
# All generators take a seed; every test must set.seed() before use.

## ---------------------------------------------------------------- tiny problem
## One structured PIV, 2 values, so the exact posterior over
## (H^A, H^B, Delta) is enumerable (13 linkage configurations x 32 true-value
## combinations). Used by test-exact-posterior.R and test-gibbs-conditionals.R.

tiny_pars <- function() {
  list(
    n_values = 2L,
    eta      = c(0.7, 0.3),                # skewed on purpose
    phi      = c(agreeA = 0.9, agreeB = 0.85, missA = 0.10, missB = 0.05),
    S        = 0.6,                        # P(true values still equal | linked)
    gamma    = 0.6,
    GA       = c(1L, 0L),                  # A record 2 has a missing value
    GB       = c(1L, 2L, 2L)
  )
}

tiny_data <- function(pars = tiny_pars()) {
  list(
    encodedA = data.frame(V1 = pars$GA, date = 0),
    encodedB = data.frame(V1 = pars$GB, date = 0),
    n_values = pars$n_values,
    PIVs_config = list(V1 = list(dynamics = "structured",
                                 bound_mistakes = c(NA, NA),
                                 fix_mistakes = c(NA, NA))),
    same_mistakes = TRUE
  )
}

tiny_model <- function() {
  survival_model("custom",
                 S     = function(X, alpha, times) rep(0.6, length(times)),
                 n_par = function(n_cov) 1L,
                 init  = function(n_cov) 0)
}

## Exact per-record registered-value likelihood P(G = g | H = h)
tiny_phi_logp <- function(g, h, agree, miss, n) {
  if (is.na(g) || g == 0) return(log(miss))
  if (g == h) return(log((1 - miss) * agree))
  log((1 - miss) * (1 - agree) / (n - 1))
}

## Exact joint log-likelihood of (Delta, H^A, H^B) given registered data,
## following the paper's complete-data likelihood:
##   L_phi(G^A|H^A) * L_phi(G^B|H^B) * [eta[h^A] * S-part]_linked
##   * eta[h^A]_nonlinked-A * eta[h^B]_nonlinked-B * L_gamma(Delta)
tiny_loglik <- function(linksA, linksB, hA, hB, pars = tiny_pars()) {
  n <- pars$n_values
  linked <- rep(FALSE, length(pars$GA))
  linked[linksA] <- TRUE
  linkedB <- rep(FALSE, length(pars$GB))
  linkedB[linksB] <- TRUE
  ll <- 0
  # eta on H^A for all A records (linked + non-linked)
  ll <- ll + sum(log(pars$eta)[hA])
  # eta on H^B for non-linked B records only
  if (any(!linkedB)) ll <- ll + sum(log(pars$eta)[hB[!linkedB]])
  # survival part for linked pairs: S if equal, (1-S)/(n-1) otherwise
  for (l in seq_along(linksA)) {
    a <- hA[linksA[l]]; b <- hB[linksB[l]]
    ll <- ll + log(if (a == b) pars$S else (1 - pars$S) / (n - 1))
  }
  # registered-value likelihoods
  ll <- ll + sum(mapply(function(g, h) tiny_phi_logp(g, h, pars$phi["agreeA"], pars$phi["missA"], n), pars$GA, hA))
  ll <- ll + sum(mapply(function(g, h) tiny_phi_logp(g, h, pars$phi["agreeB"], pars$phi["missB"], n), pars$GB, hB))
  # L_gamma: Bernoulli row sums + number of arrangements n_B!/(n_B - k)!
  k <- length(linksA)
  ll <- ll + k * log(pars$gamma) + (length(pars$GA) - k) * log(1 - pars$gamma)
  if (k > 0) {
    ll <- ll - (lgamma(length(pars$GB) + 1) - lgamma(length(pars$GB) - k + 1))
  }
  ll
}

## All 13 one-to-one linkage configurations on the tiny problem:
## list of (linksA, linksB) index vectors.
tiny_configs <- function() {
  cfgs <- list(list(integer(0), integer(0)))
  for (j in 1:3) cfgs[[length(cfgs) + 1]] <- list(1L, j)           # A1 -> Bj
  for (j in 1:3) cfgs[[length(cfgs) + 1]] <- list(2L, j)           # A2 -> Bj
  for (j in 1:3) for (k in setdiff(1:3, j))                         # A1->Bj, A2->Bk
    cfgs[[length(cfgs) + 1]] <- list(c(1L, 2L), c(j, k))
  cfgs
}

## Exact posterior over the 13 configurations (sums out H^A, H^B).
tiny_exact_posterior <- function(pars = tiny_pars()) {
  cfgs <- tiny_configs()
  p <- numeric(length(cfgs))
  for (c_i in seq_along(cfgs)) {
    acc <- 0
    for (hA in 1:2) for (hB in 1:2) for (hB2 in 1:2) for (hB3 in 1:2) {
      acc <- acc + exp(tiny_loglik(cfgs[[c_i]][[1]], cfgs[[c_i]][[2]],
                                   c(hA, hB2), c(hB, hB2, hB3), pars))
    }
    p[c_i] <- acc
  }
  p / sum(p)
}

## Run the package's own Gibbs sampler (simulateH / simulateD alternating,
## exactly as StEM's gibbs_step does) on the tiny problem and return the
## empirical frequencies of the 13 configurations.
tiny_gibbs_frequencies <- function(n_sweeps, thin = 1L, pars = tiny_pars()) {
  data <- tiny_data(pars)
  mod  <- tiny_model()
  eta   <- list(pars$eta)
  phi   <- list(pars$phi)
  alpha <- list(0)
  gamma <- pars$gamma
  linksCpp <- matrix(integer(0), 0, 2)
  linksR   <- matrix(integer(0), 0, 2)
  sumRowD  <- rep(0L, length(pars$GA))
  sumColD  <- rep(0L, length(pars$GB))
  survivalpSameH <- matrix(1, 0, 1)
  key_of <- function(linksR) {
    if (nrow(linksR) == 0) return("0")
    ord <- order(linksR[, 1])
    paste(apply(linksR[ord, , drop = FALSE], 1, paste, collapse = ">"),
          collapse = ";")
  }
  counts <- new.env(parent = emptyenv())
  for (s in seq_len(n_sweeps)) {
    h <- simulateH(data, linksCpp, survivalpSameH, sumRowD, sumColD, eta, phi)
    d <- simulateD(data, linksR, h$truepivsA, h$truepivsB,
                   gamma, eta, alpha, phi, mod)
    linksCpp <- d$links
    linksR   <- linksCpp + 1L
    sumRowD  <- d$sumRowD
    sumColD  <- d$sumColD
    survivalpSameH <- matrix(1, nrow(linksR), 1)
    if (nrow(linksR) > 0) {
      times <- abs(data$encodedB[linksR[, 2], "date"] - data$encodedA[linksR[, 1], "date"])
      survivalpSameH[, 1] <- mod$S(cbind(rep(1, nrow(linksR))), alpha[[1]], times)
    }
    if (s %% thin == 0L) {
      key <- key_of(linksR)
      counts[[key]] <- if (is.null(counts[[key]])) 1L else counts[[key]] + 1L
    }
  }
  # NB: eapply() returns entries in hash order; unlist() keeps those names.
  # Do NOT relabel with ls(counts) (sorted) - that permutes key -> value.
  freq <- unlist(eapply(counts, function(x) x))
  freq / sum(freq)
}

## ---------------------------------------------------------------- std problem
## Vignette-style 4-PIV problem (stable, stable, flexible, structured),
## moderate size, used for pipeline tests.

std_config <- function() {
  list(
    V1 = list(dynamics = "stable",    bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
    V2 = list(dynamics = "stable",    bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
    V3 = list(dynamics = "flexible",  bound_mistakes = c(NA, NA),     fix_mistakes = c(NA, NA)),
    V4 = list(dynamics = "structured", bound_mistakes = c(NA, NA),    fix_mistakes = c(0.03, 0.03),
              cond_hazard_cov = list(cov1 = c("Xe", "Xf"), cov2 = c()))
  )
}

std_generate <- function(n = c(150, 200), n_links = 100, seed = 11) {
  PIVs_config <- std_config()
  n_values <- c(5, 6, 7, 12)
  p_mistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
  p_missing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
  cond_hazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = log(c(0.7, 0.6, 0.5)))
  set.seed(seed)
  gen <- simulate_data(PIVs_config, n_values, n, n_links,
                       p_mistake, p_missing, cond_hazard_params, TRUE)
  prep <- prepare_data(gen$data1, gen$data2, "fileA", "fileB", PIVs_config,
                       uniq_id = "entity_id")
  list(gen = gen, prep = prep, PIVs = names(PIVs_config), PIVs_config = PIVs_config)
}

## Small StEM settings for fast pipeline tests.
std_stem_args <- function() {
  list(StEM_iter = 10, StEM_burnin = 5, gibbs_iter = 8, gibbs_burnin = 4,
       n_post_sample = 50)
}
