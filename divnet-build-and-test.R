library(devtools)
library(magrittr)
library(phyloseq)

# install_github("adw96/breakaway")
library(breakaway)

devtools::test(pkg)
dv = divnet(Lee %>% tax_glom("Phylum"))
dv
dv$shannon[[1]]
dv$simpson
dv$`shannon-variance`
dv$`simpson-variance`

## For confirming build
roxygenise(pkg)
document(pkg)
build(pkg, vignettes = F)
install(pkg)

library(DivNet)
test(pkg)


load_all()

set.seed(1)
nn <- 44 
qq <- 
mu = c(0.5, -5, 2) %>% to_composition_matrix()
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


