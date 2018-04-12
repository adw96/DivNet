#' LNM.EM
#'
#' This function estimates the LNM model fit from Xia et al.
#' 
#' @author Bryan Martin
#' @author Amy Willis
#' 
#' @param W count matrix, with columns listing the OTUs and rows listing the samples
#' @param X covariate matrix (optional)
#' @param tuning settings for tuning the MC-MH algorithm. Options include NULL (defaults to "fast"), "fast", "careful" or a named list with components EMiter (number of EM iterations; 6 for fast, 10 for careful), EMburn (number of EM iterations to burn; 3 for fast, 5 for careful), MCiter (number of MC iterations; 500 for fast, 1000 for careful), MCburn (number of MC iterations to burn; 250 for fast, 500 for careful) and stepsize (variance used for MH samples; 0.01 for both fast and careful)
#' @param perturbation size of purturbation for zero counts, defaults to NULL (0.05)
#' @param network How to estimate network. Defaults to NULL (the default), "default" (generalised inverse, aka naive). Other options include "diagonal", "stars" (requires glasso and SpiecEasi to be installed), or a function that you want to use to estimate the network
#' @param base OTU index to be used for base. if NULL, will use most common taxon
#' @param ncores number of cores to use for MH sampling, defaults to 1
#' @param ... additional arguments to be supplied to the network function
#' 
#' @importFrom magrittr "%>%"
#' @export
fit_aitchison <- function(W, 
                          X = NULL, 
                          tuning = NULL,
                          perturbation = NULL, 
                          network = NULL,
                          base = NULL,
                          ncores = NULL,
                          ...) {
  W <- as.matrix(W)
  output_list <- list()
  
  ### ROWS COLUMNS N is samples, Q is OTUs
  N <- nrow(W)
  Q <- ncol(W)
  
  # covariate matrix: confirm full rank
  if (is.null(X)) {
    X <- matrix(1, nrow=N, ncol=1)
  }
  # if an intercept to column is included in the covariate matrix, remove it
  intercept_columns <- apply(X, 2, function(x) max(x) == min(x))
  X <- X[ , -intercept_columns] %>% as.matrix
  
  output_list$X <- X
  no_covariates <- ncol(X) == 0
  
  # base taxon
  if (is.null(base)) {
    taxa_sums <- apply(W, 2, sum)
    base <- which(rank(taxa_sums) == max(rank(taxa_sums)))[1]
  }
  
  # take log ratio, adding perturbation term to zero counts 
  if (is.null(perturbation)) perturbation <- 0.05
  Y.p <- toLogRatios(W = W, base = base, perturbation = perturbation) # (N x Q-1) 
  
  # initialize EM algorithm
  b0 <- attr(Y.p, "center")  # this is just column means (of OTUs)
  eY <- tcrossprod(rep(1, N), b0) 
  if (!no_covariates)  {
    b <- OLS(X, Y.p)
    eY <- eY + X %*% b
  }
  sigma <- var(Y.p - eY) # (Q-1) x (Q-1)
  
  

  
  ## set up tuning parameters for EM-MH algorithm
  if (is.null(tuning)) tuning <- "fast"
  if (tuning == 0) tuning <- "fast"
  if (tuning == "fast") {
    EMiter <- 6
    EMburn <- 3
    MCiter <- 500
    MCburn <- 250
    stepsize <- 0.01
  } else if (tuning == "careful") {
    EMiter <- 10
    EMburn <- 5
    MCiter <- 1000
    MCburn <- 500
    stepsize <- 0.01
  } else {
    if (class(tuning) != "list") stop("Input parameter tuning should be `fast`, `careful`, or a list")
    EMiter <- ifelse(is.null(tuning$EMiter), 6, tuning$EMiter)
    EMburn <- ifelse(is.null(tuning$EMburn), 3, tuning$EMburn)
    MCiter <- ifelse(is.null(tuning$MCiter), 500, tuning$MCiter)
    MCburn <- ifelse(is.null(tuning$MCburn), 250, tuning$MCburn)
    stepsize <- ifelse(is.null(tuning$stepsize), 0.01, tuning$stepsize)
  }
  
  if (is.null(network)) network <- "default"
  if (is.null(ncores)) ncores <- 1
  
  
  ## set up storage for EM algorithm
  #b0.list <- b0
  b0.list <- matrix(0, nrow = EMiter + 1, ncol = length(b0))
  b0.list[1,] <- b0
  if (!no_covariates) {
    if (is.matrix(b)) {
      b.list <- array(0, dim = c(dim(b), EMiter + 1))
      b.list[,, 1] <- b
    } else {
      b.list <- matrix(0, nrow = length(b), ncol = EMiter + 1)
      b.list[, 1] <- b
    }
  }
  
  
  # sigma's by arraw, acomb3
  sigma.list <- array(0, dim = c(dim(sigma), EMiter + 1))
  sigma.list[,, 1] <- sigma
  accept.list <- matrix(0, nrow = EMiter, ncol = N)
  
  # start EM algorithm
  pb <- txtProgressBar(min = 0, max = EMiter, style = 3)
  for (em in 1:EMiter) {
    setTxtProgressBar(pb, em-1)
    start <- proc.time()
    
    # MC step
    MCarray <- MCmat(Y = Y.p, W = W, eY = eY, N = N, Q = Q, base = base, sigma = sigma, MCiter = MCiter, 
                     stepsize = stepsize, network = network, ncores = ncores, ...)
    
    # MC burn-in
    burnt <- MCarray[(MCburn + 1):MCiter, , ]
    
    Y.new <- t(apply(burnt, 3, colMeans))
    accepts <- Y.new[, 1] #  first column is acceptance ratio
    Y.new <- Y.new[, 2:Q]
    
    # update b0, means across OTUs
    b0 <- apply(Y.new, 2, mean)
    
    # update sigma
    sigSumFun <- function(i) {
      return(crossprod(t(MCarray[i, 2:Q, 1:N]) - eY))
    }
    sigSum <- foreach(i = (MCburn + 1):MCiter, .combine = "+") %do% sigSumFun(i)
    sigma <- sigSum/(N * (MCiter - MCburn))
    
    # update b
    eY <- tcrossprod(rep(1, N), b0) 
    if (!no_covariates)  {
      b <- OLS(X, Y.new)
      eY <- eY + X %*% b
    }
    
    ### STORE after updating
    ## @Bryan TODO: can this be cleaned up?
    accept.list[em,] <- accepts
    #b0.list <- rbind(b0.list, b0)
    b0.list[em + 1,] <- b0
    if (!no_covariates) {
      if (is.matrix(b)) {
        b.list[,, em + 1] <- b
      } else {
        b.list[, em + 1] <- b
      }
    }
    sigma.list[,, em + 1] <- sigma
    end <- proc.time()
  }
  setTxtProgressBar(pb, EMiter)
  cat("\n")
  
  ## Next: take average of EM samples past burn. Included init, so have EMiter+1 total 
  b0.EM <- colMeans(b0.list[(EMburn + 1):(EMiter + 1), ])
  
  output_list$beta0 <- b0.EM
  
  ### fitted value of Y
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
  
  # burn-in sigma and take average; fine b/c class of positive def matrices is closed under addition
  sigma.EM <- apply(sigma.list[, , (EMburn + 1):(EMiter + 1)], c(1, 2), mean)
  
  # store parameters
  output_list$sigma <- sigma.EM
  output_list$fitted_y <- fitted_y
  output_list$fitted_z <- toCompositionMatrix(fitted_y, base=base)
  
  output_list
}
