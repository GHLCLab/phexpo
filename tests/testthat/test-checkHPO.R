library(testthat)
context("checkHPO")

test_that("checkHPO returns a tibble data frame", {
expect_equal(class(checkHPO(list("Rickets","Breast carcinoma","Preeclampsia"))), c("tbl_df", "tbl", "data.frame"))
})

test_that("checkHPO returns a tibble data frame of 2 columns", {
  expect_equal(ncol(checkHPO(list("Rickets","Breast carcinoma","Preeclampsia"))), 2)
})

test_that("checkHPO column names are HPO_Name and Supported", {
  expect_equal(colnames(checkHPO(list("Rickets","Breast carcinoma","Preeclampsia"))), c("HPO_Name","Supported"))
})

test_that("checkHPO returns yes for known supported HPO term", {
  expect_equal(checkHPO(list("Rickets"))$Supported, "Yes")
})

test_that("checkHPO returns no for known not supported HPO term - Social and occupational deterioration", {
  expect_equal(checkHPO(list("Social and occupational deterioration"))$Supported, "No")
})
