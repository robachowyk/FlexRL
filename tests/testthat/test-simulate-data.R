# simulate_data(): the data generator must respect the configured dynamics.

test_that("REGRESSION (FINDINGS #B): changes hit the structured PIV, not the last column", {
  # Only V2 is structured; V1 and V3 are stable; zero mistakes/missingness.
  # Before the fix, all changes leaked into the LAST column (V3).
  cfg <- list(
    V1 = list(dynamics = "stable",     bound_mistakes = c(NA, NA), fix_mistakes = c(NA, NA)),
    V2 = list(dynamics = "structured", bound_mistakes = c(NA, NA), fix_mistakes = c(0.02, 0.02),
              cond_hazard_cov = list(cov1 = c(), cov2 = c())),
    V3 = list(dynamics = "stable",     bound_mistakes = c(NA, NA), fix_mistakes = c(NA, NA)))
  n_values  <- c(5, 6, 7)
  p_mistake <- list(V1 = c(0, 0), V2 = c(0, 0), V3 = c(0, 0))
  p_missing <- list(V1 = c(0, 0), V2 = c(0, 0), V3 = c(0, 0))
  cond      <- list(V1 = c(), V2 = log(0.5), V3 = c())
  set.seed(71)
  gen <- simulate_data(cfg, n_values, c(200, 250), 150, p_mistake, p_missing,
                       cond, TRUE, survival_model("exponential"))
  d1 <- gen$data1[1:150, ]
  d2 <- gen$data2[1:150, ]
  expect_equal(mean(d1$V1 != d2$V1), 0)          # stable: identical
  expect_gt(mean(d1$V2 != d2$V2), 0.05)          # structured: changes happen here
  expect_equal(mean(d1$V3 != d2$V3), 0)          # stable: must NOT change
  expect_true(all(d2$change[1:150] == (d1$V2 != d2$V2)))
})

test_that("simulate_data validates its inputs", {
  cfg <- list(V1 = list(dynamics = "stable", bound_mistakes = c(NA, NA),
                        fix_mistakes = c(NA, NA)))
  expect_error(simulate_data(cfg, 4, c(50, 40), 10, list(V1 = c(0, 0)),
                             list(V1 = c(0, 0)), list(V1 = c()), TRUE),
               "n_records\\[2\\]")
  expect_error(simulate_data(cfg, 4, c(40, 50), 60, list(V1 = c(0, 0)),
                             list(V1 = c(0, 0)), list(V1 = c()), TRUE),
               "n_links")
})

test_that("simulate_data output structure: entity ids, change flags, survival probs", {
  cfg <- list(
    V1 = list(dynamics = "stable", bound_mistakes = c(0.1, 0.1), fix_mistakes = c(NA, NA)))
  set.seed(72)
  gen <- simulate_data(cfg, 4, c(60, 80), 40, list(V1 = c(0.02, 0.02)),
                       list(V1 = c(0.01, 0.01)), list(V1 = c()), TRUE)
  expect_equal(nrow(gen$data1), 60L)
  expect_equal(nrow(gen$data2), 80L)
  # entity_id: shared ids for linked pairs, disjoint otherwise
  expect_false(any(duplicated(gen$data1$entity_id)))
  expect_false(any(duplicated(gen$data2$entity_id)))
  shared <- intersect(gen$data1$entity_id, gen$data2$entity_id)
  expect_equal(length(shared), 40L)
  # linked pairs occupy the first n_links rows of both files
  expect_setequal(gen$data1$entity_id[1:40], shared)
  expect_setequal(gen$data2$entity_id[1:40], shared)
  # proba_same_H: for a stable-only config there is no structured PIV, so no
  # survival probabilities are simulated (the field is NA by contract).
  expect_true(is.na(gen$proba_same_H))
})
