library(testthat)
context("perfFishTestHPOSingle")

test_that("perfFishTestHPOSingle returns a tibble dataframe", {
expect_equal(class(perfFishTestHPOSingle(HPO = "Rickets",enrich_1S = TRUE)), c("tbl_df", "tbl", "data.frame"))
})

test_that("perfFishTestHPOSingle returns a tibble dataframe of 10 columns", {
  expect_equal(ncol(perfFishTestHPOSingle(HPO = "Rickets",enrich_1S = TRUE)), 10)
})

test_that("perfFishTestHPOSingle returns error for HPO not supported", {
  expect_error(perfFishTestHPOSingle(HPO = "Rockets",enrich_1S = TRUE))
})


test_that("perfFishTestHPOSingle returns error for enricH_1S isn't true", {
  expect_error(perfFishTestHPOSingle(HPO = "Rickets",enrich_1S = "ERROR"))
})

test_that("perfFishTestHPOSingle returns same values for two different capitalisations", {
  expect_true(all.equal(perfFishTestHPOSingle(HPO = "Rickets",enrich_1S = TRUE),perfFishTestHPOSingle(HPO = "RICKETS",enrich_1S = TRUE)))
})

test_that("perfFishTestHPOSingle returns different values for two-sided", {
  expect_false(isTRUE(all.equal(perfFishTestHPOSingle(HPO = "Rickets",enrich_1S = TRUE),perfFishTestHPOSingle(HPO = "Rickets",enrich_1S = FALSE))))
})

