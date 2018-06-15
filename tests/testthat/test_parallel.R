library(DivNet)
context("Test parallelisation")

set.seed(1)
my_counts <- matrix(rpois(30, lambda=10), nrow = 6)
my_counts
my_covariate <- cbind(1, rep(c(0,1), each = 3), rep(c(0,1), 3))
my_covariate

# This cannot be tested with Travis or CRAN
# 
test_that("parallel works", {
  expect_is(divnet(my_counts, my_covariate,
                   variance="parametric",
                   nsub = 3, B = 2, ncores = 4,
                   tuning="test"), "list")
  expect_is(divnet(my_counts, my_covariate,
                   variance="nonparametric", ncores = 4,
                   nsub = 3, B = 2, tuning="test"), "list")

  expect_is(phylodivnet(lp, 
                        "type", 
                        c("t1.txt", "t2.txt"), 
                        ncores = 4, 
                        tuning = "test",
                        B = 2), 
            "list")
})
