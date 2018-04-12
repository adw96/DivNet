#' LNM.EM
#'
#' This function estimates the LNM model fit from Xia et al.
#' 
#' @author Bryan Martin
#' @author Amy Willis
#' 
#' @param W count matrix, with OTUs as columns
#' @param X covariate matrix (optional)
#' @param EMiter number of EM iterations, defaults to 10
#' @param EMburn number of EM iterations to burn, defaults to 5
#' @param MCiter number MC iterations, defaults to 1000
#' @param MCburn number of MC iterations to burn, defaults to 500
#' @param stepsize variance used for MH samples, defaults to 0.01. Tweak to adjust acceptance ratio
#' @param perturbation size of purturbation used for toLogRatios, defaults to 0.05
#' @param network How to estimate network. Defaults to "default" (generalised inverse, aka naive). Other options include "diagonal", "stars" (requires glasso and SpiecEasi to be installed), or a function that you want to use to estimate the network
#' @param base OTU index to be used for base. if NULL, will use most common taxon
#' @param ncores number of cores to use, defaults to 1
#' @param ... additional arguments to be supplied to the network function
#' 
#' @importFrom magrittr "%>%"
#' @export
fit_aitchison <- function(W, 
                          X = NULL, 
                          EMiter = 10, EMburn = 5, MCiter = 1000, MCburn = 500, 
                          stepsize = 0.01, 
                          perturbation = 0.05, 
                          network = "default",
                          base = NULL,
                          ncores = 1,
                          ...) {
  X_original <- X
  W <- as.matrix(W)
  
  ### ROWS COLUMNS N is samples, Q is OTUs
  N <- nrow(W)
  Q <- ncol(W)
  if (is.null(base)) {
    taxa_sums <- apply(W, 2, sum)
    base <- which(rank(taxa_sums) == max(rank(taxa_sums)))[1]
  }

  output_list <- list()
  
  if (is.null(X)) {
    X <- matrix(1, nrow=N, ncol=1)
  }
  # if an intercept to column is included in the covariate matrix, remove it
  intercept_columns <- apply(X, 2, function(x) max(x) == min(x))
  X <- X[ , -intercept_columns] %>% as.matrix
  
  output_list$X <- X
  
  
  no_covariates <- ncol(X) == 0
  
  # purturbed Y (N x Q-1) - function in getData.R
  Y.p <- toLogRatios(W, base = base, p = perturbation)
  # this is just column means (of OTUs)
  b0 <- attr(Y.p, "center")
  
  eY <- tcrossprod(rep(1, N), b0) 
  
  if (!no_covariates)  {
    b <- OLS(X, Y.p)
    eY <- eY + X %*% b
  }
  
  # (Q-1) x (Q-1)
  sigma <- var(Y.p - eY)
  
  
  #### STORING RESULTS b0's by row of a matrix, rbind
  b0.list <- b0
  if (!no_covariates)   b.list <- b
  # sigma's by arraw, acomb3
  sigma.list <- sigma
  accept.list <- c()
  
  # Should be (MCiter x Q x N) (1000 x 75 x 119) Dont forget, first column is acceptance (ie want 74 x
  # 119 for data)
  pb <- txtProgressBar(min = 0, max = EMiter, style = 3)
  for (em in 1:EMiter) {
    setTxtProgressBar(pb, em-1)
    start <- proc.time()
    MCarray <- MCmat(Y = Y.p, W = W, eY = eY, N = N, Q = Q, base = base, sigma = sigma, MCiter = MCiter, 
                     stepsize = stepsize, network = network, ncores = ncores, ...)
    
    # should call 119 apply functions, each time, get 75 means. each of iteration values for OTU.  ORIGINAL
    # for(i in 1:119) { Y.new[i,] <- apply(MCarray[(MCburn+1):MCiter,,i],2,mean) } ALTERNATIVE, no for loop
    # if large N, but transpose: seems faster by system.time
    burnt <- MCarray[(MCburn + 1):MCiter, , ]
    
    Y.new <- t(apply(burnt, 3, colMeans))
    
    # recall first column
    accepts <- Y.new[, 1]
    Y.new <- Y.new[, 2:Q]
    
    # update beta, sigma Recall: b0 is means across OTUs
    b0 <- apply(Y.new, 2, mean)
    
    # sigSums <- matrix(0,Q-1,Q-1) # Get all sample matrices for(i in (MCburn+1):MCiter) { # first part is
    # Y sample Eps.samp <- t(MCarray[i,2:Q,1:N]) - eY sigSums <- sigSums + crossprod(Eps.samp) } CAN EASILY
    # MAKE ABOVE APPLY: next just take mean
    sigSumFun <- function(i) {
      return(crossprod(t(MCarray[i, 2:Q, 1:N]) - eY))
    }
    
    sigSum <- foreach(i = (MCburn + 1):MCiter, .combine = "+") %do% sigSumFun(i)
    
    sigma <- sigSum/(N * (MCiter - MCburn))
    
    eY <- tcrossprod(rep(1, N), b0) 
    
    if (!no_covariates)  {
      b <- OLS(X, Y.new)
      eY <- eY + X %*% b
    }
    
    
    ### STORE after updating
    ## Amy@BryanTODO: can this be cleaned up?
    accept.list <- rbind(accept.list, accepts)
    b0.list <- rbind(b0.list, b0)
    if (!no_covariates) {
      if (is.matrix(b)) {
        b.list <- acomb3(b.list, b)
      } else {
        b.list <- cbind(b.list, b)
      }
    }
    sigma.list <- acomb3(sigma.list, sigma)
    end <- proc.time()
  }
  setTxtProgressBar(pb, EMiter)
  cat("\n")
  
  ## Next: take average of EM samples past burn. Included init, so have EMiter+1 total Note, below doesn't
  ## have any inital input
  b0.EM <- colMeans(b0.list[(EMburn + 1):(EMiter + 1), ])
  
  output_list$beta0 <- b0.EM
  
  ### estimated value of Y
  if (no_covariates) {
    fitted_y <- matrix(b0.EM, ncol = Q-1, nrow = N, byrow=T)
  } else {
    b.list.reduced <- b.list[, , (EMburn + 1):(EMiter + 1)] 
    
    if (length(dim(b.list.reduced)) == 2) {
      b.list.reduced <- b.list.reduced %>% array(c(1, dim(b.list.reduced)))
    }
    
    b.EM <- apply(b.list.reduced, c(1, 2), mean)
    output_list$beta <- b.EM
    
    fitted_y <- X %*% b.EM + matrix(b0.EM, ncol = Q-1, nrow = N, byrow=T)
    
  }
  
  # burn-in sigma and take average 
  # it's fine to do this because class of positive definite matrices is closed under addition
  sigma.EM <- apply(sigma.list[, , (EMburn + 1):(EMiter + 1)], c(1, 2), mean)
  
  output_list$sigma <- sigma.EM
  output_list$fitted_y <- fitted_y
  output_list$fitted_z <- toCompositionMatrix(fitted_y, base=base)
  
  output_list
}
