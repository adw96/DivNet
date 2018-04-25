
DivNet
=========


[![Build Status](https://travis-ci.org/adw96/DivNet.svg?branch=master)](https://travis-ci.org/adw96/DivNet.svg?branch=master)


DivNet is a method to estimate diversity when taxa in the community cooccur via a ecological network.

DivNet is under development with release (i.e. preprint explaining the methodology) planned for mid-April 2018. Please check back again soon!
  
Willis, A. and Martin, B. (2018+) *Estimating diversity in networked ecological communities*. In Preparation.

## Installation ##


```r
library(devtools)
install_github("adw96/DivNet")
library(DivNet)
```

## Basic Usage ##

Let's simulate and run an example with no covariates (biological replicates)

```r
set.seed(1)
my_counts <- matrix(rpois(30, lambda=10), nrow = 6)
rownames(my_counts) <- paste("Sample", 1:6, sep = "")
colnames(my_counts) <- paste("Taxon", 1:5, sep = "")

divnet(my_counts)
```

Now let's add some covariates

```r
my_covariate <- cbind(1, rep(c(0,1), each = 3), rep(c(0,1), 3))

divnet(my_counts, my_covariate)
```
