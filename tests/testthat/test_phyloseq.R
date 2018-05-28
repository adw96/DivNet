library(DivNet)
library(phyloseq)
context("Test phyloseq")

test_that("phyloseq integration", {
  data(Lee)
  expect_is(phyloseq::tax_glom(Lee, taxrank="Phylum") %>% 
              divnet(tuning = "test"), "diversityEstimates")
  expect_is(phyloseq::tax_glom(Lee, taxrank="Phylum") %>% 
              divnet(X = "char", tuning = "test"), "diversityEstimates")
})
