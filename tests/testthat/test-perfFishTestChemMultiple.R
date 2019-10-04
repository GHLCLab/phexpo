library(testthat)
context("perfFishTestChemMultiple")

test_that("perfFishTestChemMultiple returns a tibble data frame", {
expect_equal(class(perfFishTestChemMultiple(cname = list("Iodine","Iodine","Iodine"),enrich_1S = TRUE)), c("tbl_df", "tbl", "data.frame"))
})

test_that("perfFishTestChemMultiple returns a tibble data frame of 10 columns", {
  expect_equal(ncol(perfFishTestChemMultiple(cname = list("Iodine","Iodine","Iodine"),enrich_1S = TRUE)), 10)
})

test_that("perfFishTestChemMultiple returns all 8301 HPO terms", {
  expect_equal(nrow(perfFishTestChemMultiple(cname = list("Iodine","Iodine","Iodine"),enrich_1S = TRUE)), 8301)
})

test_that("perfFishTestChemMultiple returns error for chemical not supported", {
  expect_error(perfFishTestChemMultiple(cname = list("Iodin","Iodine","Iodine"),enrich_1S = TRUE))
})


test_that("perfFishTestChemMultiple returns error for enrich_1S isn't true", {
  expect_error(perfFishTestChemMultiple(cname = list("Iodine","Iodine","Iodine"),enrich_1S = "ERROR"))
})

test_that("perfFishTestChemMultiple returns same values for two different capitalisations", {
  expect_true(all.equal(perfFishTestChemMultiple(cname = list("Iodine","Iodine","Iodine"),enrich_1S = TRUE),
                        perfFishTestChemMultiple(cname = list("IODINE","IODINE","IODINE"),enrich_1S = TRUE)))
})

test_that("perfFishTestChemMultiple returns different values for two-sided", {
  expect_false(isTRUE(all.equal(perfFishTestChemMultiple(cname = list("Iodine","Iodine","Iodine"),enrich_1S = TRUE),perfFishTestChemMultiple(cname = "IODINE",enrich_1S = FALSE))))
})

test_that("perfFishTestChemMultiple returns error for Nafronyl", {
  expect_error(perfFishTestChemMultiple(cname = list("Iodine", "Nafronyl"),enrich_1S = TRUE))
})
