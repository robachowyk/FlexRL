# prepare_data(): encoding, swap logic, support restriction, validation.

test_that("encoding maps sorted unique values to 1..n_k and NA to 0", {
  set.seed(61)
  d1 <- data.frame(V = factor(c("b", "a", NA, "c")), entity_id = 1:4)
  d2 <- data.frame(V = factor(c("a", "c", "b", "b")), entity_id = 5:8)
  cfg <- list(V = list(dynamics = "stable", bound_mistakes = c(NA, NA),
                       fix_mistakes = c(NA, NA)))
  prep <- suppressWarnings(prepare_data(d1, d2, "f1", "f2", cfg))
  expect_equal(prep$n_values[["V"]], 3L)
  expect_setequal(prep$encodedA$V, c(1, 2, 0, 3))    # a=1, b=2, c=3, NA=0
})

test_that("files are swapped so that encodedA is the smaller one, labels travel", {
  set.seed(62)
  big <- data.frame(V = sample(1:4, 60, TRUE), entity_id = 1:60)
  small <- data.frame(V = sample(1:4, 40, TRUE), entity_id = 1:40)
  cfg <- list(V = list(dynamics = "stable", bound_mistakes = c(NA, NA),
                       fix_mistakes = c(NA, NA)))
  # pass the LARGER file first
  prep <- suppressWarnings(prepare_data(big, small, "big", "small", cfg))
  expect_equal(nrow(prep$encodedA), 40L)   # small file became A
  expect_equal(nrow(prep$encodedB), 60L)
  expect_true(all(prep$encodedA$source == "small"))
  expect_true(all(prep$encodedB$source == "big"))
})

test_that("true_pairs reconstruct entity_id matches", {
  std <- std_generate(seed = 11)
  prep <- std$prep
  tp <- prep$true_pairs
  labA <- colnames(tp)[1]
  labB <- colnames(tp)[2]
  idA <- prep$encodedA[tp[[1]], "entity_id"]
  idB <- prep$encodedB[tp[[2]], "entity_id"]
  expect_true(all(idA == idB))
  expect_equal(length(idA), sum(table(c(idA, idB)) == 2))  # all shared entities
})

test_that("out-of-support records are dropped with a warning, or kept when restriction is off", {
  d1 <- data.frame(V = c(1, 2, 3), entity_id = 1:3)
  d2 <- data.frame(V = c(1, 2, 9), entity_id = 4:6)
  cfg <- list(V = list(dynamics = "stable", bound_mistakes = c(NA, NA),
                       fix_mistakes = c(NA, NA)))
  expect_warning(prep <- prepare_data(d1, d2, "f1", "f2", cfg), "common support")
  expect_false(9 %in% prep$encodedB$V)
  expect_equal(nrow(prep$encodedB), 2L)
  prep2 <- suppressWarnings(prepare_data(d1, d2, "f1", "f2", cfg,
                                         restrict_support_intersection = FALSE))
  # with restriction off, all rows are kept; note the encoding is fitted on the
  # pooled values, so 9 is remapped to the highest level (4)
  expect_equal(nrow(prep2$encodedB), 3L)
  expect_equal(max(prep2$encodedB$V), 4L)
})

test_that("more than 100 unique values in a PIV stops with an error", {
  d1 <- data.frame(V = 1:150, entity_id = 1:150)
  d2 <- data.frame(V = 1:150, entity_id = 1:150)
  cfg <- list(V = list(dynamics = "stable", bound_mistakes = c(NA, NA),
                       fix_mistakes = c(NA, NA)))
  expect_error(prepare_data(d1, d2, "f1", "f2", cfg), "more than 100 unique values")
})

test_that("PIV config validation rejects malformed entries", {
  d1 <- data.frame(V = 1:5, entity_id = 1:5)
  d2 <- data.frame(V = 1:5, entity_id = 1:5)
  base <- function(x) list(V = x)
  bad1 <- base(list(dynamics = "stable", bound_mistakes = c(0.1, 0.1, 0.1),
                    fix_mistakes = c(NA, NA)))
  expect_error(prepare_data(d1, d2, "f1", "f2", bad1), "bound_mistakes")
  bad2 <- base(list(dynamics = "stable", bound_mistakes = c(0.1, 0.1, 0.1),
                    fix_mistakes = c(NA, NA)))
  expect_error(prepare_data(d1, d2, "f1", "f2", bad2), "length-2")
  bad3 <- base(list(dynamics = "stable", bound_mistakes = c(0.1, 0.1),
                    fix_mistakes = c(0.5, 0.5, 0.5)))
  expect_error(prepare_data(d1, d2, "f1", "f2", bad3), "length-2")
})

test_that("character numerics are rejected (FINDINGS #G, fixed)", {
  # `is_numeric_or_na()` used to grep "^[0-9.]+$" on the coerced string, so
  # c("0.02","0.02") passed validation and later broke arithmetic in StEM
  # ("non-numeric argument to binary operator" in `1 - PIVs_fixMistakes[, k]`).
  d1 <- data.frame(V = 1:5, entity_id = 1:5)
  d2 <- data.frame(V = 1:5, entity_id = 1:5)
  cfg <- list(V = list(dynamics = "stable", bound_mistakes = c(NA, NA),
                       fix_mistakes = c("0.02", "0.02")))
  expect_error(prepare_data(d1, d2, "f1", "f2", cfg), "numeric/NA")
  # plain NA vectors and numerics (with NA entries) stay valid
  ok1 <- list(V = list(dynamics = "stable", bound_mistakes = c(NA, NA),
                       fix_mistakes = c(NA, NA)))
  expect_no_error(suppressWarnings(prepare_data(d1, d2, "f1", "f2", ok1)))
  ok2 <- list(V = list(dynamics = "stable", bound_mistakes = c(0.1, NA),
                       fix_mistakes = c(0.02, NA)))
  expect_no_error(suppressWarnings(prepare_data(d1, d2, "f1", "f2", ok2,
                                                same_mistakes = FALSE)))
})

test_that("structured PIV config: cov1/cov2 become covA/covB and swap with the files", {
  set.seed(63)
  d_one <- data.frame(V = sample(1:4, 40, TRUE), Xe1 = rnorm(40), date = 0,
                      entity_id = 1:40)
  d_two <- data.frame(V = sample(1:4, 60, TRUE), Xe2 = rnorm(60), date = 1,
                      entity_id = 1:60)
  cfg <- list(V = list(dynamics = "structured", bound_mistakes = c(NA, NA),
                       fix_mistakes = c(0.02, 0.01),
                       # cov1/cov2 refer to the first/second argument: data1 = d_two
                       # (covariate Xe2), data2 = d_one (covariate Xe1).
                       cond_hazard_cov = list(cov1 = "Xe2", cov2 = "Xe1")))
  # pass the LARGER file first ("two", 60 rows) to force a swap:
  # new A = "one", new B = "two". same_mistakes = FALSE is required to have
  # different per-file mistake rates at all.
  prep <- suppressWarnings(prepare_data(d_two, d_one, "two", "one", cfg,
                                        same_mistakes = FALSE))
  cc <- prep$PIVs_config$V$cond_hazard_cov
  expect_true(all(c("covA", "covB") %in% names(cc)))
  expect_identical(cc$covA, "Xe1")   # new A is "one" -> its covariate Xe1
  expect_identical(cc$covB, "Xe2")
  # fix_mistakes must be reversed with the swap: new A ("one") had 0.01
  expect_equal(prep$PIVs_config$V$fix_mistakes, c(0.01, 0.02))
})

test_that("conditional-independence warning fires for associated PIVs", {
  std <- std_generate(seed = 11)
  expect_warning(prepare_data(std$gen$data1, std$gen$data2, "f1", "f2",
                              std$PIVs_config),
                 "associated")
})
