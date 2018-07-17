library(DivNet)
context("Test output")

data(Lee)
lp <- phyloseq::tax_glom(Lee, taxrank="Phylum")
# set.seed(1)
# t1 <- phy_tree(ape::rtree(20, tip.label = rownames(tax_table(lp)), br = runif))
# t2 <- phy_tree(ape::rtree(20, tip.label = rownames(tax_table(lp)), br = runif))
# ape::write.tree(t1, "t1.txt")
# ape::write.tree(t2, "t2.txt")
# 
# 
# test_that("phylodivnet runs", {
#   expect_true(TRUE)  
#   
#   expect_is(phylodivnet(lp, 
#                         "type", 
#                         c("tests/testthat/t1.txt", "tests/testthat/testt2.txt"), 
#                         ncores = 1, 
#                         tuning = "test",
#                         B = 2), 
#             "list")
# })


