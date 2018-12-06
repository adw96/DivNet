library(DivNet)
library(magrittr)
context("Test zeroes")

test_that("Identical observations give identical estimates", {
  
  n <- 6
  n_taxa <- 6
  my_counts <- matrix(c(rep(c(rep(100, n_taxa/2), rep(0, n_taxa/2)), n/2),
                        rep(c(rep(0, n_taxa/2), rep(100, n_taxa/2)), n/2)),
                      nrow = n, ncol = n_taxa, byrow=T)
  my_discrete_covariate <- cbind(1, rep(c(0,1), each = n/2))
  
  my_counts_2 <- my_counts
  my_counts_2[1:3, ] <- 100
  dv <- divnet(my_counts[1:3, ], X=my_discrete_covariate[1:3, ], base = 1)
  estimates <- dv$shannon %>% summary %$% estimate 
  
  expect_true(max(estimates) == min(estimates))
  expect_equal(max(estimates), 
               shannon_true(c(1/3, 1/3, 1/3)), tolerance = 1e-2)
  
  dv2 <- divnet(my_counts_2, X=my_discrete_covariate, base = 1)
  estimates2 <- dv2$shannon %>% summary %$% estimate 
  expect_equal(estimates2[6], estimates[1], tolerance = 1e-2)
  expect_equal(estimates2[1], 
               shannon_true(rep(1/6, 6)), tolerance = 1e-2)
  expect_equal(estimates2[6], 
               shannon_true(rep(1/3, 3)), tolerance = 1e-2)
  
})



test_that("DivNet estimates are basically correct", {
  my_counts <- matrix(c(rep(c(rep(100, n_taxa/2), rep(0, n_taxa/2)), n/2),
                        rep(c(rep(0, n_taxa/2), rep(100, n_taxa/2)), n/2)),
                      nrow = n, ncol = n_taxa, byrow=T)
  my_discrete_covariate <- cbind(1, rep(c(0,1), each = n/2))
  
  dv <- divnet(my_counts, X=my_discrete_covariate, base = 1)
  dv$shannon
  estimates <- dv$shannon %>% summary %$% estimate 
  expect_lt(max(estimates) - min(estimates), 0.01)
  
  my_counts <- matrix(c(rep(c(rep(100, n_taxa/2), rep(1, n_taxa/2)), n/2),
                        rep(c(rep(1, n_taxa/2), rep(100, n_taxa/2)), n/2)),
                      nrow = n, ncol = n_taxa, byrow=T)
  my_counts
  my_discrete_covariate <- cbind(1, rep(c(0,1), each = n/2))
  dv <- divnet(my_counts, X=my_discrete_covariate, base = 1)
  dv$shannon
  estimates <- dv$shannon %>% summary %$% estimate 
  expect_lt(max(estimates) - min(estimates), 0.01)
  
  fa <- fit_aitchison(my_counts, X=my_discrete_covariate, base = 1)
  expect_equal(fa$fitted_z[1,1], 1/3, tol = 1e-2)
  expect_equal(fa$fitted_z[1,4], 0, tol = 1e-2)
  expect_equal(fa$fitted_z[4,1], 0, tol = 1e-2)
  expect_equal(fa$fitted_z[4,4], 1/3, tol = 1e-2)
  
  
  n <- 100
  my_counts <- matrix(c(replicate(n/2, c(rpois(n_taxa/2, 100), rpois(n_taxa/2, 1))),
                        replicate(n/2, c(rpois(n_taxa/2, 1), rpois(n_taxa/2, 100)))),
                      nrow = n, ncol = n_taxa, byrow=T)
  my_counts
  fa <- fit_aitchison(my_counts, X=rep(c(0,1), each = n/2), base = 6)
  expect_equal(fa$fitted_z[1, ] %>% shannon_true, 
               shannon_true(rep(1/3, 3)), tolerance = 5e-2)
  expect_equal(fa$fitted_z[n, ] %>% shannon_true, 
               shannon_true(rep(1/3, 3)), tolerance = 5e-2)
  
  n <- 9
  my_counts <- matrix(c(rep(c(rep(100, n_taxa/3), rep(0, n_taxa/3*2)), n/3),
                        rep(c(rep(0, n_taxa/3), rep(100, n_taxa/3), rep(0, n_taxa/3)), n/3),
                        rep(c(rep(0, n_taxa/3*2), rep(100, n_taxa/3)), n/3)),
                      nrow = n, ncol = n_taxa, byrow=T)
  my_counts
  my_discrete_covariate <- cbind(1, rep(c(0,1,0), each = n/3), rep(c(0,0,1), each = n/3))
  dv <- fit_aitchison(my_counts, X=my_discrete_covariate, base = 6)
  shannon = dv$fitted_z %>% apply(1, true_shannon)
  expect_lt(max(shannon) - min(shannon), 0.01)
  
})
