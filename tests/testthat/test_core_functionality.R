library(DivNet)
context("Test output")

set.seed(1)
n <- 6
n_taxa <- 5
my_counts <- matrix(rpois(n*n_taxa, lambda=10), nrow = n)
my_discrete_covariate <- cbind(1, rep(c(0,1), each = n/2), rep(c(0,1), n/2))
my_continuous_covariate <- rnorm(n)

test_that("errors are thrown", {
  expect_error(divnet(matrix(c(10, 20, 10, 1), nrow= 1), tuning="test"))  
  expect_error(divnet(matrix(c(10, 20, 10, 1), ncol=2), tuning="test"))
  expect_error(divnet(matrix(c(10, 20, 10, 1, 50, 0), ncol=2), tuning="test"))
  expect_is(divnet(matrix(c(10, 20, 10, 1, 50, 0), nrow=2), tuning="test"), "diversityEstimates")
})

test_that("fit_aitchison works", {
  expect_true(TRUE)  
  expect_is(fit_aitchison(my_counts, tuning = "test"), "list")
  expect_is(fit_aitchison(my_counts, my_discrete_covariate, tuning = "test"), "list")
  expect_is(fit_aitchison(my_counts, my_continuous_covariate, tuning = "test"), "list")
})

test_that("divnet works", {
  expect_is(divnet(my_counts, variance = 0, tuning="test"), "diversityEstimates")
  expect_is(divnet(my_counts, my_discrete_covariate, variance = 0, tuning="test"), "diversityEstimates")
  expect_is(divnet(my_counts, X = my_continuous_covariate, variance = 0, tuning="test"), "diversityEstimates")
})

test_that("arguments are fine", {
  expect_is(divnet(my_counts, my_discrete_covariate,
                   base = 2, perturbation = 0.01, 
                   tuning="test"), "diversityEstimates")
  expect_is(divnet(my_counts, my_discrete_covariate,
                   base = 1, perturbation = 0.1, 
                   tuning=list(EMiter = 4, EMburn = 2, MCiter = 5, MCburn = 2)), "diversityEstimates")
})


test_that("ordering with base works as expected", {
  for (i in 1:5) {
    perturbed_counts <- my_counts
    perturbed_counts[,i] <- rpois(n, lambda=200)
    perturbed_counts
    fa1 <- fit_aitchison(perturbed_counts, tuning = "test", base=i)
    estimated_compositions <- fa1$fitted_z[1, ]
    expect_equal(which.max(estimated_compositions), i)
  }
})

test_that("diversity indices work as expected", {
  expect_equal(bray_curtis_true(c(0.4, 0.6), c(0.2, 0.8)),
               bc_fast(c(0.4, 0.6), c(0.2, 0.8)))
  expect_equal(euclidean_true(c(0.4, 0.6), c(0.2, 0.8)),
               euc_fast(c(0.4, 0.6), c(0.2, 0.8)))
  
  expect_equal(shannon_true(rep(1/2, 2)), log(2))
  expect_equal(shannon_true(rep(1/5, 5)), log(5))
  
  expect_equal(simpson_true(rep(1/2, 2)), 2*0.5^2)
  expect_equal(simpson_true(rep(1/5, 5)), 5*(0.2^2))
  
  # check errors are thrown when not given a proportion
  expect_error(shannon_true(rep(1/5, 4)))
  expect_error(simpson_true(rep(1/5, 4)))
  expect_error(bc_fast(c(0.4, 0.5), c(0.2, 0.8)))
  expect_error(euc_fast(c(0.4, 0.5), c(0.2, 0.8)))
  
})

# test_that("conversions work as expected", {
#   perturbed_counts <- my_counts
#   perturbed_counts[,3] <- rpois(n, lambda=200)
#   perturbed_counts
#   
#   i <- 2
#   for (i in 1:5) {
#     fa1 <- fit_aitchison(perturbed_counts, tuning = "test", base=i)
#     fa1
#     expect_equal(fa$fitted_z, to_composition_matrix(fa$fitted_y, base = i))
#     expect_equal(fa$fitted_z[1,], to_composition(fa$fitted_y[1,], base = i))
#   }
#   # works: 1
#   # doesn't work: 2, 3, 4, 5
# })

## TODO check that changing the base works as expected

