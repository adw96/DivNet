set.seed(1)
my_counts <- matrix(rpois(30, lambda=10), nrow = 6)
my_counts
my_covariate <- cbind(1, rep(c(0,1), each = 3), rep(c(0,1), 3))
my_covariate

test_that("fitting the model is okay", {
  expect_true(TRUE)  
  expect_is(fit_aitchison(my_counts, tuning = "test"), "list")
  expect_is(fit_aitchison(my_counts, my_covariate, tuning = "test"), "list")
})

test_that("no covariates works", {
  expect_is(divnet(my_counts, variance = 0, tuning="test"), "list")
})

test_that("covariates work", {
  expect_is(divnet(my_counts, my_covariate, variance = 0, tuning="test"), "list")
  expect_is(divnet(my_counts, X = rnorm(6), variance = 0, tuning="test"), "list")
})

test_that("parametric variances", {
  expect_is(divnet(my_counts, 
                   variance="parametric",
                   nsub = 3, B = 2,
                   tuning="test"), "list")
  expect_is(divnet(my_counts, my_covariate, 
                   variance="parametric",
                   nsub = 3, B = 2, tuning="test"), "list")
  
})

test_that("nonparametric variances", {
  expect_is(divnet(my_counts, 
                   variance="nonparametric",
                   nsub = 3, B = 2,
                   tuning="test"), "list")
  expect_is(divnet(my_counts, my_covariate, 
                   variance="nonparametric",
                   nsub = 3, B = 2, tuning="test"), "list")
  
})

test_that("arguments are fine", {
  expect_is(divnet(my_counts, my_covariate,
                      base = 2, tuning="test"), "list")
})