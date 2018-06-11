library(DivNet)

set.seed(1)
my_counts <- matrix(rpois(16, lambda=10), nrow = 4)

test_that("network options works", {
  
  expect_error(divnet(my_counts, network="nonsense", variance = 0, tuning="test"))
  expect_is(divnet(my_counts, network="default", variance = 0, tuning="test"), "diversityEstimates")
  expect_is(divnet(my_counts, network="diagonal", variance = 0, tuning="test"), "diversityEstimates")
  
  # removed to speed up
  # expect_is(divnet(my_counts, network="stars", variance = 0, 
  #                              tuning="test", ncores = 1,
  #                              method='mb', lambda.min.ratio=1e-1, nlambda=2, 
  #                              icov.select.params=list(rep.num=3)), 
  #                       "diversityEstimates")
  
  
})
