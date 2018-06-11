library(DivNet)
context("Test variances")

set.seed(1)
my_counts <- matrix(rpois(30, lambda=50), nrow = 6)
my_counts
my_covariate <- cbind(1, rep(c(0,1), each = 3), rep(c(0,1), 3))
my_covariate

test_that("parametric variances", {
  
  expect_is(divnet(my_counts, 
                   variance="parametric",
                   nsub = 3, B = 2,
                   tuning="test"), "diversityEstimates")
  expect_is(divnet(my_counts, 
                   my_covariate, 
                   variance="parametric",
                   nsub = 3, B = 2, tuning="test"), "diversityEstimates")
  
})

test_that("nonparametric variances", {
  expect_is(divnet(my_counts, 
                   variance="nonparametric",
                   nsub = 3, B = 2,
                   tuning="test"), "diversityEstimates")
  expect_is(divnet(my_counts, 
                   my_covariate, 
                   variance="nonparametric",
                   nsub = 3, B = 2, tuning="test"), "diversityEstimates")
  
})

