library(testthat)
context("perfFishTestHPOMultiple")

test_that("perfFishTestHPOMultiple returns a tibble data frame", {
expect_equal(class(perfFishTestHPOMultiple(HPO = list("Rickets","Rickets","Rickets"),enrich_1S = TRUE)), c("tbl_df", "tbl", "data.frame"))
})

test_that("perfFishTestHPOMultiple returns a tibble data frame of 10 columns", {
  expect_equal(ncol(perfFishTestHPOMultiple(HPO = list("Rickets","Rickets","Rickets"),enrich_1S = TRUE)), 10)
})

test_that("perfFishTestHPOMultiple returns error for HPO not supported", {
  expect_error(perfFishTestHPOMultiple(HPO = list("Rockets","Rickets","Rickets"),enrich_1S = TRUE))
})


test_that("perfFishTestHPOMultiple returns error for enrich_1S isn't true", {
  expect_error(perfFishTestHPOMultiple(HPO = list("Rickets","Rickets","Rickets"),enrich_1S = "ERROR"))
})

test_that("perfFishTestHPOMultiple returns same values for two different capitalisations", {
  expect_true(all.equal(perfFishTestHPOMultiple(HPO = list("Rickets","Rickets","Rickets"),enrich_1S = TRUE),
                        perfFishTestHPOMultiple(HPO = list("RICKETS","RICKETS","RICKETS"),enrich_1S = TRUE)))
})

test_that("perfFishTestHPOMultiple returns different values for two-sided", {
  expect_false(isTRUE(all.equal(perfFishTestHPOMultiple(HPO = list("Rickets","Rickets","Rickets"),enrich_1S = TRUE),
                                perfFishTestHPOMultiple(HPO = list("Rickets","Rickets","Rickets"),enrich_1S = FALSE))))
})

