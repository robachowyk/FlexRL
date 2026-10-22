# link_with_* wrappers: uniform interface, index conventions, known truth.

std <- std_generate(seed = 11)
prep <- std$prep
PIVs <- std$PIVs
tp <- paste(prep$true_pairs[[1]], prep$true_pairs[[2]], sep = "_")

test_that("link_with_FlexRL returns posterior probabilities; > 0.5 is one-to-one", {
  args <- c(list(data = prep), std_stem_args())
  set.seed(81)
  res <- link_with_FlexRL(prep$encodedA, prep$encodedB, args)
  expect_true(all(res$LinkScore > 0 & res$LinkScore <= 1))
  hi <- res$LinkScore > 0.5
  expect_gt(sum(hi), 0)
  expect_false(any(duplicated(res$idxA[hi])))
  expect_false(any(duplicated(res$idxB[hi])))
})

test_that("BRL wrapper: idxA indexes dataA, idxB indexes dataB, links are mostly true", {
  set.seed(82)
  # V1 is unique per record so the true pairs are unambiguous for BRL's Gibbs;
  # V2 adds a weak secondary field. (With only a handful of levels per field,
  # ~20 distinct patterns for 120 records, the posterior is too diffuse for a
  # short chain and BRL legitimately returns almost no links.)
  B <- data.frame(V1 = 1:120, V2 = rep(1:30, 4),
                  local_id = as.character(1:120), source = "B")
  A <- B[1:60, ]
  A$V1[sample(60, 6)] <- sample(1:60, 6, TRUE)     # a few registration mistakes
  A$local_id <- as.character(1:60)
  A$source <- "A"
  args <- list(flds = c("V1", "V2"), types = c("bi", "bi"), nIter = 300)
  res <- link_with_BRL(A, B, args)
  expect_true(all(res$idxA >= 1 & res$idxA <= 60))
  expect_true(all(res$idxB >= 1 & res$idxB <= 120))
  got <- paste(res$idxA, res$idxB, sep = "_")
  true <- paste(1:60, 1:60, sep = "_")
  expect_gt(sum(got %in% true) / max(length(got), 1), 0.5)
})

test_that("fastLink wrapper returns pairs with probability scores", {
  args <- list(varnames = PIVs, tol.em = 1e-6, threshold.match = 0.5,
               return.all = TRUE, n.cores = 1)
  set.seed(83)
  res <- link_with_fastLink(prep$encodedA, prep$encodedB, args)
  expect_gt(length(res$idxA), 0)
  expect_true(all(res$LinkScore >= 0 & res$LinkScore <= 1))
})

test_that("reclin2 wrapper returns threshold-selected pairs", {
  args <- list(on = PIVs,
               formula = stats::as.formula(paste("~", paste(PIVs, collapse = "+"))),
               type = "mpost", add = TRUE, variable = "selected",
               score = "mpost", threshold = 0.5)
  set.seed(84)
  res <- link_with_reclin2(prep$encodedA, prep$encodedB, args)
  expect_gt(length(res$idxA), 0)
  expect_true(all(res$LinkScore >= 0.5))
})

test_that("fedmatch wrapper returns matches with multivar scores", {
  args <- list(by = PIVs, match_type = "multivar",
               unique_key_1 = "unique_key_A", unique_key_2 = "unique_key_B",
               multivar_settings = fedmatch::build_multivar_settings(
                 compare_type = rep("indicator", length(PIVs)),
                 wgts = rep(1 / length(PIVs), length(PIVs))))
  set.seed(85)
  res <- link_with_fedmatch(prep$encodedA, prep$encodedB, args)
  expect_gt(length(res$idxA), 0)
  expect_true(all(res$idxA >= 1 & res$idxA <= nrow(prep$encodedA)))
  expect_true(all(res$idxB >= 1 & res$idxB <= nrow(prep$encodedB)))
})

test_that("multilink wrapper (no duplicates) returns in-range one-to-one pairs", {
  args <- list(types = rep("bi", length(PIVs)),
               breaks = setNames(rep(list(NA), length(PIVs)), PIVs),
               duplicates = c(0, 0), n_iter = 200, max_cc_size = 50)
  set.seed(86)
  res <- link_with_multilink(prep$encodedA, prep$encodedB, args)
  expect_gt(length(res$idxA), 0)
  expect_true(all(res$idxA >= 1 & res$idxA <= nrow(prep$encodedA)))
  expect_true(all(res$idxB >= 1 & res$idxB <= nrow(prep$encodedB)))
  expect_false(any(duplicated(res$idxA)))
  expect_false(any(duplicated(res$idxB)))
  expect_null(res$LinkScore)
})

test_that("KNOWN BUG: multilink wrapper errors when duplicates are allowed (FINDINGS #F)", {
  args <- list(types = rep("bi", length(PIVs)),
               breaks = setNames(rep(list(NA), length(PIVs)), PIVs),
               duplicates = c(1, 0), n_iter = 200, max_cc_size = 50)
  expect_error(link_with_multilink(prep$encodedA, prep$encodedB, args),
               "pairs_to_keep")
})

test_that("diyar wrapper (vignette arguments, PIVs only) links the std data", {
  args <- list(probabilistic = TRUE, return_weights = TRUE)
  set.seed(87)
  res <- link_with_diyar(prep$encodedA[, PIVs], prep$encodedB[, PIVs], args)
  expect_null(res$LinkScore)
  expect_true(all(res$idxA >= 1 & res$idxA <= nrow(prep$encodedA)))
  expect_true(all(res$idxB >= 1 & res$idxB <= nrow(prep$encodedB)))
})

test_that("KNOWN TRAP: diyar wrapper silently uses entity_id when present (FINDINGS #E)", {
  args <- list(probabilistic = TRUE, return_weights = TRUE)
  set.seed(87)
  res <- link_with_diyar(prep$encodedA[, c(PIVs, "entity_id")],
                         prep$encodedB[, c(PIVs, "entity_id")], args)
  got <- paste(res$idxA, res$idxB, sep = "_")
  # with the ground-truth identifier treated as a PIV, essentially every
  # returned pair is a true pair (precision inflation; see FINDINGS #E)
  expect_gt(length(got), 0)
  expect_gt(sum(got %in% tp) / length(got), 0.95)
})

test_that("KNOWN BUG: diyar wrapper errors on prepare_data output containing 'change' (FINDINGS #E)", {
  expect_error(
    link_with_diyar(prep$encodedA, prep$encodedB,
                    list(probabilistic = TRUE, return_weights = TRUE)))
})
