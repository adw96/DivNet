library(DivNet)

set.seed(1)
my_counts <- matrix(rpois(16, lambda=10), nrow = 4)

test_that("network options works", {
  
  expect_error(divnet(my_counts, network="nonsense", variance = 0, tuning="test"))
  expect_is(divnet(my_counts, network="default", variance = 0, tuning="test"), "diversityEstimates")
  expect_is(divnet(my_counts, network="diagonal", variance = 0, tuning="test"), "diversityEstimates")
  expect_is(divnet(my_counts, network="stars", variance = 0, tuning="test"), "diversityEstimates")
  
  # # 15.578 s
  # t1 <- proc.time()
  # expect_is(divnet(my_counts, network="stars", variance = 0, tuning="test",
  #                  lambda.min.ratio=1e-2,
  #                  nlambda=2, 
  #                  icov.select.params=list(rep.num=5)), 
  #           "diversityEstimates")
  # proc.time() - t1
  # 
  # # 15.596 s
  # t1 <- proc.time()
  # expect_is(divnet(my_counts, network="stars", variance = 0, tuning="test",
  #                  lambda.min.ratio=1e-2,
  #                  nlambda=2, icov.select.params=list(rep.num=3)), 
  #           "diversityEstimates")
  # proc.time() - t1
  # 
  # # 15.586 s
  # t1 <- proc.time()
  # expect_is(divnet(my_counts, network="stars", variance = 0, tuning="test",
  #                  lambda.min.ratio=1e-3,
  #                  nlambda=2, icov.select.params=list(rep.num=3)), 
  #           "diversityEstimates")
  # proc.time() - t1
  # # 15.586 s
  # t1 <- proc.time()
  # expect_is(divnet(my_counts, network="stars", variance = 0, tuning="test",
  #                  lambda.min.ratio=1e-3,
  #                  nlambda=4, icov.select.params=list(rep.num=3)), 
  #           "diversityEstimates")
  # proc.time() - t1
  # 
})
