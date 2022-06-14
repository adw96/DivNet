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

test_that("errors occur in OLS", {
  X <- matrix(rnorm(30), nrow = 5)
  Y <- rnorm(rnorm(6))
  #expect_output(OLS(X, Y), "Tried to regress X on Y where X is")
  #expect_output(OLS(X, Y), print(X))
  expect_error(OLS(X, Y), "Non-conformable multiplication")
})

test_that("errors occur in pick_base", {
  expect_error(pick_base(W = matrix(0, nrow = 3, ncol = 3)),
               "Yikes! No taxa observed in all samples!\n Pick which taxon is to be the base")
})
