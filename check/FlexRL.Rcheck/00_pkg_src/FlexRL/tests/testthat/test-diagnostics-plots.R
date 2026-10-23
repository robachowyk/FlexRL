# RL_diagnostics, print/plot methods, and the standalone plots: smoke tests
# plus basic structural assertions.

std <- std_generate(seed = 11)
prep <- std$prep
PIVs <- std$PIVs
PIVs_type <- list(V1 = FALSE, V2 = FALSE, V3 = FALSE, V4 = TRUE)

fit <- suppressWarnings(do.call(StEM, c(list(data = prep), std_stem_args())))

diag <- suppressWarnings(RL_diagnostics(
  fit = fit, encodedA = prep$encodedA, encodedB = prep$encodedB,
  compare_vars = PIVs, vars_type_cont = PIVs_type, PIVs = PIVs,
  true_pairs = prep$true_pairs, FDP_estimation = TRUE,
  RL_method = "FlexRL", data = prep,
  maxIter4CV = 1, n_repeats = 2,
  StEM_iter = 6, StEM_burnin = 3, gibbs_iter = 6, gibbs_burnin = 3,
  n_post_sample = 30))

test_that("RL_diagnostics returns the expected components", {
  expect_s3_class(diag, "RL_diagnostics")
  expect_true(!is.null(diag$FDP_measures$curve))
  expect_true(all(c("RL_FDP_score", "augmRL_FDP_score", "augmRL_FDP_synth")
                  %in% names(diag$FDP_measures$curve)))
  expect_true(!is.null(diag$true_performance))
  expect_true(!is.null(diag$discrepancy_measures))
})

test_that("print.RL_diagnostics runs", {
  expect_output(print(diag), regexp = ".")
})

test_that("plot.RL_diagnostics works for all four types", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  expect_no_error(plot(diag, type = "scores"))
  expect_no_error(plot(diag, type = "distributions", threshold = 0.5))
  expect_no_error(plot(diag, type = "FDP"))
  expect_no_error(plot(diag, type = "discrepancy"))
  grDevices::dev.off()
  unlink(f)
  expect_error(plot(diag, type = "nonsense"), "type")
})

test_that("plot_StEM_convergence plots the chains", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  expect_no_error(plot_StEM_convergence(fit))
  grDevices::dev.off()
  unlink(f)
})

test_that("plot_linkage_scores and plot_distributions run", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  expect_no_error(plot_linkage_scores(nrow(prep$encodedA) * nrow(prep$encodedB),
                                      fit$Delta$x))
  dl <- list(data_baseline = prep$encodedA, data_select = prep$encodedA[1:50, ])
  expect_no_error(plot_distributions(dl, PIVs))
  grDevices::dev.off()
  unlink(f)
})
