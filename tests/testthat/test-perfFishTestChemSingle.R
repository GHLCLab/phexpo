library(testthat)
context("perfFishTestChemSingle")

test_that("perfFishTestChemSingle returns a tibble dataframe", {
expect_equal(class(perfFishTestChemSingle(cname = "Iodine",enrich_1S = TRUE)), c("tbl_df", "tbl", "data.frame"))
})

test_that("perfFishTestChemSingle returns a tibble dataframe of 10 columns", {
  expect_equal(ncol(perfFishTestChemSingle(cname = "Iodine",enrich_1S = TRUE)), 10)
})

test_that("perfFishTestChemSingle returns all 8301 HPO terms", {
  expect_equal(nrow(perfFishTestChemSingle(cname = "Iodine",enrich_1S = TRUE)), 8301)
})

test_that("perfFishTestChemSingle returns error for chemical not supported", {
  expect_error(perfFishTestChemSingle(cname = "Iodin",enrich_1S = TRUE))
})


test_that("perfFishTestChemSingle returns error for enricH_1S isn't true", {
  expect_error(perfFishTestChemSingle(cname = "Iodine",enrich_1S = "ERROR"))
})

test_that("perfFishTestChemSingle returns same values for two different capitalisations", {
  expect_true(all.equal(perfFishTestChemSingle(cname = "Iodine",enrich_1S = TRUE),perfFishTestChemSingle(cname = "IODINE",enrich_1S = TRUE)))
})

test_that("perfFishTestChemSingle returns different values for two-sided", {
  expect_false(isTRUE(all.equal(perfFishTestChemSingle(cname = "Iodine",enrich_1S = TRUE),perfFishTestChemSingle(cname = "Iodine",enrich_1S = FALSE))))
})

