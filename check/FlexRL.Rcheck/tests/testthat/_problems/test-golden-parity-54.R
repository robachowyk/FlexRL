# Extracted from test-golden-parity.R:54

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "FlexRL", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
gold <- test_path("golden")

# test -------------------------------------------------------------------------
for (id in GOLDEN_STATE_IDS) {
    st <- readRDS(file.path(gold, sprintf("state_%s.rds", id)))
    rebuilt <- golden_build_state(id, st$data)
    expect_identical(rebuilt$truepivsA, st$truepivsA, info = id)
    expect_identical(rebuilt$truepivsB, st$truepivsB, info = id)
    expect_identical(rebuilt$linksR, st$linksR, info = id)
    expect_identical(rebuilt$sumRowD, st$sumRowD, info = id)
    expect_identical(rebuilt$sumColD, st$sumColD, info = id)
    expect_identical(rebuilt$eta, st$eta, info = id)
    expect_identical(rebuilt$phi, st$phi, info = id)
    expect_identical(rebuilt$alpha, st$alpha, info = id)
  }
