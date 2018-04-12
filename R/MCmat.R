#' MCmat
#'
#' This function simulates MC step for an entire matrix. Should not need to be used by user directly. NOTE: currently written in parallel
#'
#' @author Bryan Martin
#' @author Amy Willis
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
#' @param perturbation size of purturbation used for toLogRatios, defaults to 0.05
#' @param network How to estimate network. Defaults to "default" (generalised inverse, aka naive). Other options include "diagonal", or a function that takes a sample covariance matrix and returns an estimate of the inverse covariance matrix (eg glasso or SpiecEasi)
#' @param ncores number of cores to use, defaults to 1
#' @param ... additional arguments to be supplied to the network function
#'
#' @export
MCmat <- function(Y, W, eY, N, Q, base, sigma, MCiter, stepsize = 1, perturbation = 0.05, network = "default", ncores = 1, ...) {
  
  # sigInv <- solve(sigma)
  if (network == "diagonal") {
    # take diagonal vec, make it a diagonal matrix, solve
    sigInv <- diag(1/diag(sigma))
    
    
  } else if (network == "default") {
    test <- try(chol(sigma), silent = T)
    if (class(test) == "try-error") {
      sigInv <- MASS::ginv(sigma)
    } else {
      sigInv <- chol2inv(chol(sigma))
    }
  } else if (network == "stars") {
    sigInv <- stars(sigma, W, base = base, perturbation = perturbation, ncores = ncores)
  } else { #if (class(network) %in% c("function", "methods", "standardGeneric")) {
    sigInv <- try(network(sigma, ...), silent = T)
    if (class(sigInv) == "try-error") {
      stop("Cannot use supplied network option")
    }
  }
  
  
  MH_path <- function(i) {
    MCrow(Yi = Y[i, ], Wi = W[i, ], eYi = eY[i, ], Q = Q, base = base, sigInv = sigInv, MCiter = MCiter, 
          stepsize = stepsize)
  }
  
  if (ncores > 1) {
    ####################
    # Parallel option ##
    ####################
    warning("I'm going to use all your cores!! Tell Amy to fix this")
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


stars <- function(sigma, W, base, perturbation, ncores) {
  Y.p <- toLogRatios(W, base = base, perturbation = perturbation)
  inverse_covariance_estimate <- SpiecEasi::sparseiCov(data=Y.p, method = "glasso")
  selected <- SpiecEasi::icov.select(inverse_covariance_estimate, rep.num=5, ncores = ncores)
  gl <- glasso::glasso(sigma, rho = selected$opt.lambda)
  gl$wi
}