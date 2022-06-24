#' OLS function
#'
#' Function to simulate raw counts W from model fit
#'
#' @param X explanatory variables
#' @param Y response variables
#'
OLS <- function(X, Y) {
  
  p <- ncol(X)
  centY <- scale(Y, center = TRUE, scale = FALSE)
  
  # quick function to aid apply
  aFun <- function(ycol) {
    
    outcome <- try(tcrossprod(ginv(crossprod(X)), X) %*% ycol, silent = T)
    if ("try-error" %in% class(outcome)) {
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
#' 
#' Picks the base taxon to be used in divnet fit. If no taxon is detected
#' in all samples, returns error; in this case, we recommend manually choosing 
#' a few taxa fairly abundant across samples and for each such taxon
#' fitting divnet specifying this taxon as the base taxon.
#' 
#' @param W A taxon abundance matrix (taxa as columns)
#' @param detection_cutoff The proportion of samples base taxon must be detected
#' in. Default is NULL in which case base taxon must be detected in all samples
#' unless automatic_cutoff is set to TRUE.
#' @param automatic_cutoff Choose detection cutoff automatically? Default is 
#' FALSE. If TRUE, detection_cutoff will be set equal to the maximum proportion
#' of samples any taxon is detected in.
#' @value Index corresponding to taxon chosen as base taxon
#' @author Amy Willis
#' @export
pick_base <- function(W,
                      detection_cutoff = NULL,
                      automatic_cutoff = FALSE) {
  if(is.null(detection_cutoff)& !automatic_cutoff ){
  taxa_sums <- colSums(W)
  taxa_unobserved <- apply(W, 2, function(x) ifelse(any(x == 0), 0, 1))
  if (all(taxa_unobserved == 0)) stop("Yikes! No taxa observed in all samples!\n Pick which taxon is to be the base")
  return(which.max(taxa_sums*taxa_unobserved))
  } else{
    taxa_sums <- colSums(W)
    taxa_detection_proportion <- apply(W, 2, function(x) mean(x >0))
    if(automatic_cutoff){
      detection_cutoff <- max(taxa_detection_proportion)
    }
    acceptable_taxa <- as.numeric(taxa_detection_proportion >= detection_cutoff)
    if (all(acceptable_taxa == 0)) stop(paste("Yikes! No taxa observed in at least", 100*detection_cutoff,
                                              "percent of samples.",
                                              "Set automatic_cutoff = TRUE or decrease detection_cutoff and try again."))
    return(which.max(taxa_sums*acceptable_taxa))
  }
  
}


