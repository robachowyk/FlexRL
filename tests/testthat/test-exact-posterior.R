# Exact-posterior test: the Gibbs sampler's stationary distribution must match
# the exact posterior over (Delta, H^A, H^B) for a tiny fully enumerable
# problem. This is the test that catches ANY mis-specification of the Gibbs
# conditionals (it would have caught the eta-indexing bug, FINDINGS #A).

test_that("Gibbs stationary distribution matches the exact posterior", {
  set.seed(2024)
  pars  <- tiny_pars()
  exact <- tiny_exact_posterior(pars)
  freq  <- tiny_gibbs_frequencies(n_sweeps = 30000, thin = 3L, pars = pars)

  # align empirical frequencies with the 13 enumerated configurations
  cfgs <- tiny_configs()
  keys <- vapply(cfgs, function(cfg) {
    if (length(cfg[[1]]) == 0) return("0")
    paste(apply(cbind(cfg[[1]], cfg[[2]]), 1, paste, collapse = ">"),
          collapse = ";")
  }, character(1))
  emp <- freq[match(keys, names(freq))]
  emp[is.na(emp)] <- 0

  expect_length(exact, 13L)
  expect_equal(sum(exact), 1, tolerance = 1e-12)
  diffs <- abs(emp - exact)
  expect_lt(max(diffs), 0.025)                     # per-configuration deviation
  expect_lt(sum(diffs) / 2, 0.05)                  # total variation distance
})
