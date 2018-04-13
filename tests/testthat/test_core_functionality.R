set.seed(1)
my_counts <- matrix(rpois(30, lambda=10), nrow = 6)
my_counts
my_covariate <- cbind(1, rep(c(0,1), each = 3), rep(c(0,1), 3))
my_covariate

test_that("no covariates works", {
  
  expect_equal(divnet(my_counts, 
                      tuning=list(EMiter = 4, EMburn = 2, MCiter = 50, MCburn = 20)) %>% class, "list")
  
})

test_that("covariates work", {
  expect_equal(divnet(my_counts, my_covariate, 
                      tuning=list(EMiter = 4, EMburn = 2, MCiter = 50, MCburn = 20)) %>% class, "list")
  expect_equal(divnet(my_counts, my_covariate, 
                      variance="nonparametric",
                      nsub = 3,
                      B = 2,
                      tuning=list(EMiter = 4, EMburn = 2, MCiter = 50, MCburn = 20)) %>% class, "list")
  
})

test_that("arguments are fine", {
  expect_equal(divnet(my_counts, my_covariate,
                      base = 2,
                      tuning=list(EMiter = 4, EMburn = 2, MCiter = 50, MCburn = 20)) %>% class, "list")
  
})