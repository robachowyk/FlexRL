# Structural invariants of the Gibbs sampler and the StEM output.

test_that("Gibbs sweeps keep the linkage one-to-one at all times", {
  std <- std_generate(seed = 21)
  data <- std$prep
  PIVs <- std$PIVs
  eta   <- lapply(data$n_values, function(x) rep(1 / x, x))
  phi   <- lapply(data$n_values, function(x) c(0.9, 0.9, 0.1, 0.1))
  alpha <- lapply(seq_along(PIVs), function(k) {
    if (data$PIVs_config[[k]]$dynamics == "structured") rep(log(0.05), 3) else c(-Inf)
  })
  gamma <- 0.5
  linksCpp <- matrix(integer(0), 0, 2)
  linksR   <- matrix(integer(0), 0, 2)
  sumRowD  <- rep(0L, nrow(data$encodedA))
  sumColD  <- rep(0L, nrow(data$encodedB))
  survivalpSameH <- matrix(1, 0, length(PIVs))
  mod <- survival_model("exponential")
  set.seed(31)
  for (s in 1:50) {
    h <- simulateH(data, linksCpp, survivalpSameH, sumRowD, sumColD, eta, phi)
    d <- simulateD(data, linksR, h$truepivsA, h$truepivsB,
                   gamma, eta, alpha, phi, mod)
    linksCpp <- d$links
    linksR   <- linksCpp + 1L
    sumRowD  <- d$sumRowD
    sumColD  <- d$sumColD
    survivalpSameH <- matrix(1, nrow(linksR), length(PIVs))
    if (nrow(linksR) > 0) {
      times <- abs(data$encodedB[linksR[, 2], "date"] - data$encodedA[linksR[, 1], "date"])
      dyna <- which(vapply(data$PIVs_config, function(x) x$dynamics == "structured", logical(1)))
      for (k in dyna) {
        cov_k <- cbind(rep(1, nrow(linksR)),
                       data$encodedA[linksR[, 1], data$PIVs_config[[k]]$cond_hazard_cov$covA, drop = FALSE],
                       data$encodedB[linksR[, 2], data$PIVs_config[[k]]$cond_hazard_cov$covB, drop = FALSE])
        survivalpSameH[, k] <- mod$S(cov_k, alpha[[k]], times)
      }
    }
    expect_true(all(sumRowD <= 1L), info = paste("sweep", s, "row sums"))
    expect_true(all(sumColD <= 1L), info = paste("sweep", s, "col sums"))
    expect_equal(sum(sumRowD), nrow(linksR))
    expect_equal(sum(sumColD), nrow(linksR))
    expect_true(!any(duplicated(linksR[, 1])) && !any(duplicated(linksR[, 2])))
  }
})

test_that("StEM Delta output is a valid posterior summary", {
  std <- std_generate(seed = 21)
  set.seed(41)
  fit1 <- suppressWarnings(do.call(StEM, c(list(data = std$prep), std_stem_args())))
  set.seed(41)
  fit2 <- suppressWarnings(do.call(StEM, c(list(data = std$prep), std_stem_args())))

  d <- fit1$Delta
  expect_true(all(d$x > 0 & d$x <= 1))
  # at most one pair per A-record / B-record above the 0.5 decision threshold
  hi <- d[d$x > 0.5, ]
  expect_false(any(duplicated(hi$i)))
  expect_false(any(duplicated(hi$j)))
  # reproducibility under a fixed seed
  expect_identical(fit1$Delta, fit2$Delta)

  # posterior mass per A-record must not exceed 1
  mass_per_row <- tapply(d$x, d$i, sum)
  expect_true(all(mass_per_row <= 1 + 1e-12))
})

test_that("StEM accepts gibbs_burnin = 0 (auto burn-in) on the tiny problem", {
  data <- tiny_data()
  args <- list(data = data, StEM_iter = 3, StEM_burnin = 1,
               gibbs_iter = 4, gibbs_burnin = 0, n_post_sample = 10)
  set.seed(51)
  fit <- suppressWarnings(do.call(StEM, args))
  expect_equal(length(fit$gamma), 3L)
  expect_true(all(fit$Delta$x > 0 & fit$Delta$x <= 1))
})
