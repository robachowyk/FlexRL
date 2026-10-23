# naive_linkage(): deterministic exact-match linkage.

test_that("naive_linkage matches records agreeing on all non-missing PIVs", {
  A <- data.frame(V1 = c(1, 2, 0), V2 = c(1, 3, 2))    # 0 = missing
  B <- data.frame(V1 = c(1, 3, 2), V2 = c(1, 2, 3))
  res <- naive_linkage(c("V1", "V2"), A, B)
  got <- paste(res$idxA, res$idxB, sep = "_")
  # B = (1,1), (3,2), (2,3); A2=(2,3) matches B3; A3 (V1 missing, V2=2)
  # matches B2 via na_match = TRUE (default)
  expect_setequal(got, c("1_1", "2_3", "3_2"))
})

test_that("na_match = TRUE lets a missing PIV match any value", {
  A <- data.frame(V1 = c(1, 0), V2 = c(1, 2))
  B <- data.frame(V1 = c(1, 2, 3), V2 = c(1, 2, 1))
  res <- naive_linkage(c("V1", "V2"), A, B, na_match = TRUE)
  got <- paste(res$idxA, res$idxB, sep = "_")
  expect_setequal(got, c("1_1", "2_2"))   # A2 has V1 missing -> matches B2 (V2=2)
  res2 <- naive_linkage(c("V1", "V2"), A, B, na_match = FALSE)
  got2 <- paste(res2$idxA, res2$idxB, sep = "_")
  expect_setequal(got2, "1_1")            # A2's V1=0 must match V1=0 only: none
})

test_that("na_match is forced to FALSE (with a warning) for a single PIV", {
  A <- data.frame(V1 = c(1, 0))
  B <- data.frame(V1 = c(1, 2))
  expect_warning(res <- naive_linkage(c("V1"), A, B, na_match = TRUE), "na_match")
  got <- paste(res$idxA, res$idxB, sep = "_")
  expect_setequal(got, "1_1")
})

test_that("na_is_zero = FALSE treats NA as missing instead of 0", {
  A <- data.frame(V1 = c(1, NA))
  B <- data.frame(V1 = c(1, 2))
  res <- naive_linkage(c("V1"), A, B, na_is_zero = FALSE)
  got <- paste(res$idxA, res$idxB, sep = "_")
  expect_setequal(got, "1_1")
})
