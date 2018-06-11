#' OLS function
#'
#' Function to simulate raw counts W from model fit
#'
#' @param X explanatory variables
#' @param Y response variables
#'
OLS <- function(X, Y) {
  
  p <- ncol(X)
  centY <- scale(Y, center=TRUE, scale = FALSE)
  
  # quick function to aid apply
  aFun <- function(ycol) {
    
    outcome <- try(tcrossprod(ginv(crossprod(X)), X) %*% ycol, silent = T)
    if (class(outcome) == "try-error") {
      cat("Tried to regress X on Y where X is\n")
      print(X)
      cat("and Y is \n")
      print(Y)
      stop("Non-conformable multiplication")
    } else {
      return(outcome)
    }
  }
  b <- apply(centY, 2, aFun)
  matrix(b, nrow = p)
}


#' acomb3
#'
#' This function works like cbind, but along the third dimension of an array
#'
#' @param ... Array
#' 
acomb3 <- function(...) abind(..., along = 3)

#' pick_base
pick_base <- function(W) {
  taxa_sums <- colSums(W)
  taxa_unobserved <- apply(W, 2, function(x) ifelse(any(x == 0), 0, 1))
  if (all(taxa_unobserved == 0)) stop("Yikes! No taxa observed in all samples!\n Pick which taxon is to be the base")
  which.max(taxa_sums*taxa_unobserved)
}


#' Make design matrix
#' 
#' @param phyloseq_object A phyloseq object
#' @param variables variable names
#' 
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq get_variable
#' 
#' @export
make_design_matrix <- function(phyloseq_object, variables) {
  predictors <- phyloseq_object %>% sample_data %>% get_variable(variables)
  model.matrix( ~predictors, data = predictors %>% as.data.frame)
}