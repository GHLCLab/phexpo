library(testthat)
context("checkChemical")

test_that("checkChemical returns a tibble data frame", {
expect_equal(class(checkChemical(list("Ethanol","Iodine","Zinc"))), c("tbl_df", "tbl", "data.frame"))
})

test_that("checkChemical returns a tibble data frame of 2 columns", {
  expect_equal(ncol(checkChemical(list("Ethanol","Iodine","Zinc"))), 2)
})

test_that("checkChemical column names are Chemical_Name and Supported", {
  expect_equal(colnames(checkChemical(list("Ethanol","Iodine","Zinc"))), c("Chemical_Name", "Supported"))
})

test_that("checkChemical returns yes for known supported chemicals", {
  expect_equal(checkChemical(list("Ethanol"))$Supported, "Yes")
})

test_that("checkChemical returns no for known not supported chemicals - Nafronyl", {
  expect_equal(checkChemical(list("Nafronyl"))$Supported, "No")
})
