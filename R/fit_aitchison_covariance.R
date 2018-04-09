#' @export
fit_aitchison_covariance <- function(W, X = NULL, 
                                     EMiter = 10, EMburn = 5, MCiter = 1000, MCburn = 500, 
                                     stepsize = 0.01, 
                                     perturbation = 0.05, 
                                     poorman = FALSE, 
                                     interval = 0.95,
                                     base = NULL,
                                     in_parallel = TRUE,
                                     sigInvFn) {
  X_original <- X
  W <- as.matrix(W)
  
  ### ROWS COLUMNS N is samples, Q is OTUs
  N <- nrow(W)
  Q <- ncol(W)
  if (is.null(base)) {
    base <- Q ## this is compatible with later
  }
  
  if (is.null(X)) {
    X <- matrix(1, nrow=N, ncol=1)
  }
  intercept_columns <- apply(X, 2, function(x) max(x) == min(x))
  X <- X[ , -intercept_columns] %>% as.matrix
  
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
    MCarray <- MCmat_covariance(Y = Y.p, W = W, eY = eY, N = N, Q = Q, base = base, sigma = sigma, MCiter = MCiter, 
                                stepsize = stepsize, poorman = poorman, in_parallel = in_parallel, sigInvFn=sigInvFn)
    
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
  # accept.EM <- colMeans(accept.list[(EMburn):(EMiter), ])
  b0.EM <- colMeans(b0.list[(EMburn + 1):(EMiter + 1), ])
  
  output_list <- list()
  output_list$beta0 <- b0.EM
  
  ### estimated value of Y
  if (no_covariates) {
    fitted_y <- matrix(b0.EM, ncol = Q-1, nrow = N, byrow=T)
  } else {
    b.list.reduced <- b.list[, , (EMburn + 1):(EMiter + 1)] 
    
    if (length(dim(b.list.reduced)) == 2) {
      b.list.reduced %>% array(c(1, dim(b.list.reduced))) -> b.list.reduced
    }
    
    b.EM <- apply(b.list.reduced, c(1, 2), mean)
    output_list$beta <- b.EM
    
    fitted_y <- X %*% b.EM + matrix(b0.EM, ncol = Q-1, nrow = N, byrow=T)
    
  }
  
  sigma.EM <- apply(sigma.list[, , (EMburn + 1):(EMiter + 1)], c(1, 2), mean)
  
  output_list$X <- X
  output_list$sigma <- sigma.EM
  output_list$fitted_y <- fitted_y
  
  
  # compute estimate  
  exp_b <- fitted_y %>% exp
  denominator <- 1 + apply(exp_b, 1, sum)
  zz <- cbind(exp_b, 1) / matrix(denominator, nrow=N, ncol = Q, byrow=F)
  
  output_list$fitted_z <- zz
  
  
  # estimate Var(Z) = Var(h(Y))
  h_function <- function(zz) exp(zz)/(exp(zz) %>% sum + 1)
  delta_h <-  function(zz) (-(exp(zz))^2 + exp(zz)*(exp(zz) %>% sum + 1))/(exp(zz) %>% sum + 1)^2 
  
  empirical_delta_h <- apply(exp_b / denominator, 2, delta_h)
  variance_y <- sigma.EM
  variance_z <- crossprod(empirical_delta_h %*% variance_y, empirical_delta_h) 
  output_list$variance_z <- variance_z
  
  output_list
}

#' @export
MCmat_covariance <- function(Y, W, eY, N, Q, base, sigma, MCiter, stepsize = 1, poorman = FALSE, in_parallel = TRUE, sigInvFn) {
  
  sigInv <- sigInvFn(sigma, W)
  
  
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