library(DivNet)
library(magrittr)
context("Test formula")

test_that("Test formula implementation", {
  my_counts <- matrix(c(rep(c(rep(100, n_taxa/2), rep(1, n_taxa/2)), n/2),
                        rep(c(rep(1, n_taxa/2), rep(100, n_taxa/2)), n/2)),
                      nrow = n, ncol = n_taxa, byrow=T)
  my_counts
  my_discrete_covariate <- cbind(1, rep(c(0,1), each = n/2))
  colnames(my_discrete_covariate) <- c("X0", "X1")
  dv <- divnet(my_counts, X = my_discrete_covariate,formula = ~ X1, base = 1)
  dv$shannon
  estimates <- dv$shannon %>% summary %$% estimate 
  expect_lt(max(estimates) - min(estimates), 0.01)
  
  multiple_covs <- data.frame("X1" = rnorm(n), 
                              "X2" = rnorm(n), 
                              "X3" = rnorm(n))
  dv <- divnet(my_counts, X = multiple_covs, 
               formula = ~ X1 + X2 + X3, base = 1)
  expect_is(dv, "diversityEstimates")
  dv2 <- divnet(my_counts, X = multiple_covs, 
               formula = ~ X1*X3, base = 1)
  expect_is(dv2, "diversityEstimates")
  dv3 <- divnet(my_counts, X = multiple_covs, 
                formula = ~ X1 - 1, base = 1)
  expect_is(dv3, "diversityEstimates")
  
}
)