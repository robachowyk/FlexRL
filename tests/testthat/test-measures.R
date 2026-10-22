# Discrepancy/agreement measures: mmd, support_iou, prop_level, smd, RL_agreement.

test_that("mmd is ~0 for identical samples, symmetric, and grows with separation", {
  set.seed(71)
  x <- as.data.frame(matrix(rnorm(200), 100, 2))
  expect_lt(abs(mmd(x, x)), 0.05)      # biased V-statistic: small positive bias
  y <- as.data.frame(matrix(rnorm(200), 100, 2))
  expect_equal(mmd(x, y), mmd(y, x))
  expect_gt(mmd(x, x + 3), mmd(x, x + 0.1))
})

test_that("KNOWN LIMITATION: mmd returns NaN for constant data (sigma = 0)", {
  x <- data.frame(V = rep(1, 10))
  expect_true(is.nan(suppressWarnings(mmd(x, x))))
})

test_that("support_iou is 1 for identical supports, 0 for disjoint", {
  h1 <- graphics::hist(1:10, breaks = 0:10, plot = FALSE)
  h2 <- graphics::hist(1:10, breaks = 0:10, plot = FALSE)
  expect_equal(support_iou(h1, h2), 1)
  h3 <- graphics::hist(20:30, breaks = 19:31, plot = FALSE)
  expect_equal(support_iou(h1, h3), 0)
  h4 <- graphics::hist(5:14, breaks = 4:15, plot = FALSE)
  expect_equal(support_iou(h1, h4), 6 / 15)
})

test_that("prop_level ignores NAs", {
  df <- data.frame(v = c("a", "a", NA, "b"))
  expect_equal(prop_level(df, "v", "a"), 2 / 3)
  expect_equal(prop_level(df, "v", "b"), 1 / 3)
})

test_that("smd (continuous) is (mean_select - mean_base) / sd_base", {
  set.seed(72)
  base <- data.frame(age = rnorm(200, 40, 10))
  sel  <- data.frame(age = rnorm(80, 43, 10))
  res <- smd(sel, base, "age")
  expected <- (mean(sel$age) - mean(base$age)) / sd(base$age)
  expect_equal(res[[1]], expected, tolerance = 1e-12)
})

test_that("smd (categorical) standardises level-proportion differences", {
  base <- data.frame(sex = factor(c(rep("M", 60), rep("F", 40))))
  sel  <- data.frame(sex = factor(c(rep("M", 30), rep("F", 30))))
  res <- smd(sel, base, "sex", continuous = FALSE)
  expect_equal(res[["sex_M"]], (0.5 - 0.6) / sd(base$sex == "M"))
  expect_equal(res[["sex_F"]], (0.5 - 0.4) / sd(base$sex == "F"))
})

test_that("RL_agreement computes per-variable agreement rates", {
  d1 <- data.frame(V1 = c(1, 2, 3, NA), V2 = c(1, 1, 2, 2))
  d2 <- data.frame(V1 = c(1, 2, 9, 5),  V2 = c(1, 2, 2, 2))
  # pair (4, 3) compares d1 row 4 (V1 = NA) with d2 row 3
  pairs <- data.frame(idx1 = 1:4, idx2 = c(1, 2, 3, 3))
  res <- RL_agreement(d1, d2, c("V1", "V2"), pairs)
  expect_equal(res$agreements[["V1"]], 2 / 3)   # NA pair excluded (na.rm = TRUE)
  expect_equal(res$agreements[["V2"]], 3 / 4)
  res2 <- RL_agreement(d1, d2, c("V1"), pairs, na.rm = FALSE, na_match = TRUE)
  expect_equal(res2$agreements[["V1"]], 3 / 4)  # NA counted as a match
  res3 <- RL_agreement(d1, d2, c("V1"), pairs, na.rm = FALSE, na_match = FALSE)
  expect_equal(res3$agreements[["V1"]], 2 / 4)
})
