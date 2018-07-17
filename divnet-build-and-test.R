rm(list = ls())
pkg <- "/Users/adw96/Documents/software/DivNet/"
setwd(pkg)
# setwd("/Users/adw96/Documents/software/")
library(devtools)
library(roxygen2)
library(magrittr)
library(phyloseq)

# install_github("adw96/breakaway")
library(breakaway)

## For confirming build
roxygenise(pkg)
document(pkg)
build(pkg, vignettes = F)
install(pkg)
library(DivNet)
test(pkg)

set.seed(1)
mu = c(0.5, -5, 2) %>% to_composition()
A <- matrix(runif(16, -1, 1), ncol = 4)
D <- diag(seq(from = 10, to = 0.01, length.out = 4))
Sigma <- A %*% D %*% t(A)
mu
Sigma
make_w(mu = rbind(mu, mu, mu), 
       Sigma = Sigma, mm = 1e5)

build(pkg, vignettes = F)
install(pkg)
library(DivNet)
test_file("DivNet/tests/testthat/test_core_functionality.R")
