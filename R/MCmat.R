#' MCmat
#'
#' This function simulates MC step for an entire matrix. Should not need to be used by user directly. NOTE: currently written in parallel
#'
#' @param Y logratio matrix
#' @param W corresponding count matrix
#' @param eY current expected value of logratio matrix
#' @param N number of samples, or nrow of Y
#' @param Q number of OTUs minus the base, or ncol of Y
#' @param base OTU index used for base
#' @param sigma current estimate of sigma
#' @param MCiter number of MC samples to generate
#' @param stepsize variance used for MH samples, defaults to 1. Tweak to adjust acceptance ratio
#' @param poorman boolean of whether to just use simple diagonal inverse to calculate sigma inverse. Defaults to FALSE
#'
#' @export
MCmat <- function(Y, W, eY, N, Q, base, sigma, MCiter, stepsize = 1, poorman = FALSE, in_parallel = TRUE) {
  
  # sigInv <- solve(sigma)
  if (poorman == TRUE) {
    # take diagonal vec, make it a diagonal matrix, solve
    sigInv <- diag(1/diag(sigma))
    
    
  } else {
    test <- try(chol(sigma), silent = T)
    if (class(test) == "try-error") {
      sigInv <- MASS::ginv(sigma)
    } else {
      sigInv <- chol2inv(chol(sigma))
    }
  }
  MH_path <- function(i) {
    MCrow(Yi = Y[i, ], Wi = W[i, ], eYi = eY[i, ], Q = Q, base = base, sigInv = sigInv, MCiter = MCiter, 
          stepsize = stepsize)
  }
  
  if (in_parallel) {
    ####################
    # Parallel option ##
    ####################
    registerDoParallel(detectCores())
    Y.MH <-  foreach(i = 1:N, .combine = "acomb3", .multicombine = TRUE) %dopar% {
      MH_path(i)
    }
    stopImplicitCluster()
  } else {
    ####################
    ## Series option ###
    ####################
    Y.MH <-  foreach(i=1:N, .combine='acomb3', .multicombine=TRUE) %do% {
      MH_path(i)
    }
  }
  
  # Should be (MCiter x Q x N) Dont forget, first column is acceptance
  return(Y.MH)
}
