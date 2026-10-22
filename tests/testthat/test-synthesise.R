# synthesise(): augmentation with synthetic records.

std <- std_generate(seed = 11)
PIVs <- std$PIVs
encA <- std$prep$encodedA[, c(PIVs, "local_id", "source")]
encB <- std$prep$encodedB[, c(PIVs, "local_id", "source")]

test_that("synthesise keeps real rows (up to support restriction) and marks synthetic ones", {
  set.seed(91)
  nd <- suppressWarnings(synthesise("arf", encA, encB, PIVs, 50))
  # Real rows are never altered or invented: what survives is a subset of the
  # input (prepare_data's support restriction may drop a few rows).
  real <- nd$dataB[nd$dataB$source != "synthetic", ]
  expect_true(all(real$local_id %in% encB$local_id))
  expect_true(nrow(nd$dataA) <= nrow(encA))
  expect_gt(sum(nd$dataB$source == "synthetic"), 0)
  synth <- nd$dataB[nd$dataB$source == "synthetic", ]
  expect_true(all(synth$local_id > nrow(encB)))
})

test_that("restrict_support_intersection = FALSE never drops rows", {
  set.seed(92)
  nd <- suppressWarnings(synthesise("arf", encA, encB, PIVs, 200,
                                    restrict_support_intersection = FALSE))
  expect_equal(nrow(nd$dataB), nrow(encB) + 200L)
})

test_that("synthesise reports how many synthetic rows survived (post-fix contract)", {
  set.seed(93)
  nd <- suppressWarnings(synthesise("arf", encA, encB, PIVs, 50))
  expect_true("n_synth_actual" %in% names(nd))
  expect_equal(nd$n_synth_actual, sum(nd$dataB$source == "synthetic"))
})

test_that("synthpop and mice paths run on small data", {
  set.seed(94)
  nd1 <- suppressWarnings(synthesise("synthpop", encA, encB, PIVs, 30))
  expect_gt(sum(nd1$dataB$source == "synthetic"), 0)
  set.seed(95)
  nd2 <- suppressWarnings(synthesise("mice", encA, encB, PIVs, 30))
  expect_gt(sum(nd2$dataB$source == "synthetic"), 0)
})

test_that("high-cardinality PIVs take the numeric-conversion path without error", {
  set.seed(96)
  # synthpop::syn needs at least two columns, so use two PIVs.
  A <- data.frame(V1 = sample.int(65, 100, TRUE), V2 = sample.int(65, 100, TRUE),
                  local_id = as.character(1:100), source = "A")
  B <- data.frame(V1 = sample.int(65, 150, TRUE), V2 = sample.int(65, 150, TRUE),
                  local_id = as.character(1:150), source = "B")
  nd <- suppressWarnings(synthesise("synthpop", A, B, c("V1", "V2"), 20))
  expect_gt(sum(nd$dataB$source == "synthetic"), 0)
})
