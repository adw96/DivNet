data(Lee)
dn <- phyloseq::tax_glom(Lee, taxrank="Phylum") %>% 
  divnet(tuning = "test")

#test_that("make alpha estimates works", {
#  expect_type(make_alpha_estimates(dn), "list")
#})
