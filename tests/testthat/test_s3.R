library(DivNet)
context("Test output")

data(Lee)
divnet_phylum_char <- phyloseq::tax_glom(Lee, taxrank="Phylum") %>%
  divnet(X = "char", tuning = "test")

test_that("hypothesis testing works", {
  expect_is(testDiversity(divnet_phylum_char), "matrix")
})


test_that("plotting works", {
  expect_is(plot(divnet_phylum_char), "ggplot")
})

test_that("printing works", {
  expect_null(print(divnet_phylum_char))
})
