library(DivNet)
context("Test output")

data(Lee)
lp <- phyloseq::tax_glom(Lee, taxrank="Phylum")
t1 <- phy_tree(ape::rtree(20, tip.label = rownames(tax_table(lp)), br = runif))
t2 <- phy_tree(ape::rtree(20, tip.label = rownames(tax_table(lp)), br = runif))
ape::write.tree(t1, "t1.txt")
ape::write.tree(t2, "t2.txt")


test_that("UniFrac works", {
  expect_true(TRUE)  
  
  
  
  expect_is(phylodivnet(lp, 
                        "type", 
                        c("t1.txt", "t2.txt"), 
                        ncores = 1, 
                        tuning = "test"), 
            "list")
})
# 
# test_that("divnet works", {
#   expect_is(divnet(my_counts, variance = 0, tuning="test"), "diversityEstimates")
#   expect_is(divnet(my_counts, my_discrete_covariate, variance = 0, tuning="test"), "diversityEstimates")
#   expect_is(divnet(my_counts, X = my_continuous_covariate, variance = 0, tuning="test"), "diversityEstimates")
# })
# 
# test_that("arguments are fine", {
#   expect_is(divnet(my_counts, my_discrete_covariate,
#                    base = 2, perturbation = 0.01, 
#                    tuning="test"), "diversityEstimates")
#   expect_is(divnet(my_counts, my_discrete_covariate,
#                    base = 1, perturbation = 0.1, 
#                    tuning=list(EMiter = 4, EMburn = 2, MCiter = 5, MCburn = 2)), "diversityEstimates")
# })
# 
# 
# test_that("ordering with base works as expected", {
#   for (i in 1:5) {
#     perturbed_counts <- my_counts
#     perturbed_counts[,i] <- rpois(n, lambda=200)
#     perturbed_counts
#     fa1 <- fit_aitchison(perturbed_counts, tuning = "test", base=i)
#     estimated_compositions <- fa1$fitted_z[1, ]
#     expect_equal(which.max(estimated_compositions), i)
#   }
# })
# 
# test_that("diversity indices work as expected", {
#   expect_equal(bray_curtis_true(c(0.4, 0.6), c(0.2, 0.8)),
#                bc_fast(c(0.4, 0.6), c(0.2, 0.8)))
#   expect_equal(euclidean_true(c(0.4, 0.6), c(0.2, 0.8)),
#                euc_fast(c(0.4, 0.6), c(0.2, 0.8)))
#   
#   expect_equal(shannon_true(rep(1/2, 2)), log(2))
#   expect_equal(shannon_true(rep(1/5, 5)), log(5))
#   
#   expect_equal(simpson_true(rep(1/2, 2)), 2*0.5^2)
#   expect_equal(simpson_true(rep(1/5, 5)), 5*(0.2^2))
#   
#   # check errors are thrown when not given a proportion
#   expect_error(shannon_true(rep(1/5, 4)))
#   expect_error(simpson_true(rep(1/5, 4)))
#   expect_error(bc_fast(c(0.4, 0.5), c(0.2, 0.8)))
#   expect_error(euc_fast(c(0.4, 0.5), c(0.2, 0.8)))
#   
# })
