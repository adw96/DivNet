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

test_that("error when sample names have non-standard characters", {
  upd_divnet_phylum_char <- divnet_phylum_char
  colnames(upd_divnet_phylum_char$`bray-curtis`)[1] <- "error!"
  expect_error(simplifyBeta(dv = upd_divnet_phylum_char,
                            physeq = Lee_phylum,
                            measure = "bray-curtis",
                            x = "char"),
               "Your sample names contain non-standard characters, please change that :)")
})

test_that("Aitchison distance works", {
  expect_is(simplifyBeta(dv = divnet_phylum_char,
                         physeq = Lee_phylum,
                         measure = "Aitchison",
                         x = "char"),
            "data.frame")
})
