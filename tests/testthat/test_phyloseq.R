library(DivNet)
library(phyloseq)
context("Test phyloseq")

test_that("phyloseq integration", {
  data(Lee)
  expect_is(phyloseq::tax_glom(Lee, taxrank="Phylum") %>% 
              divnet(tuning = "test"), "diversityEstimates")
  expect_is(phyloseq::tax_glom(Lee, taxrank="Phylum") %>% 
              divnet(X = "char", tuning = "test"), "diversityEstimates")
  
  
  top_five_taxa = names(sort(taxa_sums(Lee), decreasing = TRUE)[1:5])
  
  dv <- Lee %>% prune_taxa(top_five_taxa, .) %>% divnet(base = 1, tuning = "test")
  dv2 <- Lee %>% prune_taxa(top_five_taxa, .) %>% otu_table %>% divnet(base = 1, tuning = "test")
  
  expect_equal(ncol(dv$fitted_z), 5)
  expect_equal(ncol(dv2$fitted_z), 5) # should be correct dimension with 5 columns
  # changed with new X default, sample specific covariates
  expect_false(isTRUE(dv2$fitted_z[1,1] == dv2$fitted_z[2,1]))
  expect_false(isTRUE(dv2$fitted_z[1,5] == dv2$fitted_z[2,5]))
  
})
