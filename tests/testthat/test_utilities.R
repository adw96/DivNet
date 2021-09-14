library(DivNet)
context("Test utility functions")


data(Lee)
Lee_phylum <- phyloseq::tax_glom(Lee, taxrank="Phylum")
divnet_phylum_char <- divnet(Lee_phylum, X = "char", tuning = "test")


test_that("OLS function returns something reasonable",
          {
expect_is(OLS(X = model.matrix(lm(as.numeric(sample_data(Lee_phylum)$type)~
                          sample_data(Lee_phylum)$char)),
    Y = cbind(sapply(divnet_phylum_char$shannon,
                     function(x) x$estimate),
              sapply(divnet_phylum_char$simpson,
                     function(x) x$estimate))),
    "matrix")
          })

test_that("pickbase() picks a base in Lee phylum data",
          {
expect_is(pick_base(t(otu_table(Lee_phylum))),
          "integer")
          })
