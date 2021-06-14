library(DivNet)
library(magrittr)
context("Test zeroes")

set.seed(1)
n <- 6
n_taxa <- 6
my_counts <- matrix(c(rep(c(rep(100, n_taxa/2), rep(0, n_taxa/2)), n/2),
                      rep(c(rep(0, n_taxa/2), rep(100, n_taxa/2)), n/2)),
                    nrow = n, ncol = n_taxa, byrow=T)
my_counts
my_discrete_covariate <- cbind(1, rep(c(0,1), each = n/2))


test_that("DivNet estimates are basically correct", {
  
  dv <- divnet(my_counts[1:3, ], X=my_discrete_covariate[1:3, ], base = 1)
  estimates <- dv$shannon %>% summary %$% estimate 
  expect_true(max(estimates) == min(estimates))
  expect_equal(max(estimates), 
               shannon_true(c(1/3, 1/3, 1/3)), tolerance = 1e-2)
  
  my_counts_2 <- my_counts
  my_counts_2[1:3, ] <- 100
  my_counts_2
  dv2 <- divnet(my_counts_2, X=my_discrete_covariate, base = 1)
  estimates2 <- dv2$shannon %>% summary %$% estimate
  estimates <- unname(estimates)
  estimates2 <- unname(estimates2)
  expect_equal(estimates2[6], estimates[1], tolerance = 1e-2)
  expect_equal(estimates2[1], 
               shannon_true(rep(1/6, 6)), tolerance = 1e-2)
  expect_equal(estimates2[6], 
               shannon_true(rep(1/3, 3)), tolerance = 1e-2)
  
})



test_that("Identical observations give identical estimates", {
  
  dv <- divnet(my_counts, X=my_discrete_covariate, base = 1)
  
  estimates <- dv$shannon %>% summary %$% estimate 
  expect_lt(max(estimates) - min(estimates), 0.01)
  
  W <- matrix(pmax(1, my_counts), nrow = 6)
  X = my_discrete_covariate
  tuning = NULL
  perturbation = NULL 
  network = NULL
  base = 1
  ncores = NULL
  
  fa <- fit_aitchison(W, X=my_discrete_covariate, base = 6)
  expect_equal(fa$fitted_z[1,1], 1/3, tolerance = 1e-2)
  expect_equal(fa$fitted_z[6,6], 1/3, tolerance = 1e-2)
  
})