library(DivNet)
context("Test simplification of beta diversities")

data(Lee)
Lee_phylum <- phyloseq::tax_glom(Lee, taxrank="Phylum")
divnet_phylum_char <- divnet(Lee_phylum, X = "char", tuning = "test")

test_that("simplifyBeta returns a dataframe",
          {
expect_is(simplifyBeta(dv = divnet_phylum_char,
             physeq = Lee_phylum,
             measure = "bray-curtis",
             x = "char"),
      "data.frame")
          })
