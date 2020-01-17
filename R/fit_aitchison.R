#' fit_aitchison
#'
#' This function estimates the parameters of the log-ratio model
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
#' @export
fit_aitchison <- function(W, 
                          X = NULL, 
                          tuning = NULL,
                          perturbation = NULL, 
                          network = NULL,
                          base = NULL,
                          ncores = NULL,
                          ...) {
  output_list <- list()
  
  ### ROWS COLUMNS N is samples, Q is OTUs
  N <- nrow(W)
  Q <- ncol(W)
  
  # covariate matrix: confirm full rank
  if (is.null(X)) {
    X <- matrix(1, nrow=N, ncol=1)
  }
  if (length(X) == N) {
    X <- matrix(X, nrow=N, ncol=1)
  }
  # if an intercept to column is included in the covariate matrix, remove it
  intercept_columns <- apply(X, 2, function(x) max(x) != min(x))
  X <- X[ , intercept_columns] %>% as.matrix
  
  no_covariates <- ncol(X) == 0
  X_c <- scale(X, scale = FALSE)
  output_list$X <- X_c
  
  X_c <- scale(X, scale = FALSE)
  
  # base taxon
  if (is.null(base)) { 
    base <- pick_base(W)
  }
  
  # take log ratio, adding perturbation term to zero counts 
  if (is.null(perturbation)) perturbation <- 0.05
  Y_p <- to_log_ratios(W = W, base = base, perturbation = perturbation) # (N x Q-1) 
  
  # initialize EM algorithm
  b0 <- colMeans(Y_p)
  eY <- tcrossprod(rep(1, N), b0) 
  if (!no_covariates)  {
    b <- OLS(X_c, Y_p)
    eY <- eY + X_c %*% b
  }
  sigma <- var(Y_p - eY) # (Q-1) x (Q-1)
  

  
  ## set up tuning parameters for EM-MH algorithm
  if (is.null(tuning)) tuning <- "fast"
  if (class(tuning) == "list") {
    EMiter <- ifelse(is.null(tuning$EMiter), 6, tuning$EMiter)
    EMburn <- ifelse(is.null(tuning$EMburn), 3, tuning$EMburn)
    MCiter <- ifelse(is.null(tuning$MCiter), 500, tuning$MCiter)
    MCburn <- ifelse(is.null(tuning$MCburn), 250, tuning$MCburn)
    stepsize <- ifelse(is.null(tuning$stepsize), 0.01, tuning$stepsize)
  } else if (tuning == "careful") {
    EMiter <- 10
    EMburn <- 5
    MCiter <- 1000
    MCburn <- 500
    stepsize <- 0.01
  } else if (tuning == "test") {
    EMiter <- 4
    EMburn <- 2
    MCiter <- 5
    MCburn <- 2
    stepsize <- 0.05
  } else {
    EMiter <- 6
    EMburn <- 3
    MCiter <- 500
    MCburn <- 250
    stepsize <- 0.01
  }
  
  if (is.null(network)) network <- "default"
  if (is.null(ncores)) ncores <- 1
  
  
  ## set up storage for EM algorithm
  #b0_list <- b0
  b0_list <- matrix(0, nrow = EMiter + 1, ncol = length(b0))
  b0_list[1,] <- b0
  if (!no_covariates) {
    if (is.matrix(b)) {
      b_list <- array(0, dim = c(dim(b), EMiter + 1))
      b_list[,, 1] <- b
    } else {
      b_list <- matrix(0, nrow = length(b), ncol = EMiter + 1)
      b_list[, 1] <- b
    }
  }
  
  
  # sigma's by arraw, acomb3
  sigma_list <- array(0, dim = c(dim(sigma), EMiter + 1))
  sigma_list[, , 1] <- sigma
  # accept_list <- matrix(0, nrow = EMiter, ncol = N)
  
  # start EM algorithm
  pb <- utils::txtProgressBar(min = 0, max = EMiter, style = 3)
  for (em in 1:EMiter) {
    utils::setTxtProgressBar(pb, em-1)
    #start <- proc.time()

    if (network == "diagonal") {
      sigInv <- diagonal_network(sigma)
    } else if (network == "default") {
      sigInv <- default_network(sigma)
    } else if (network == "stars") {
      sigInv <- stars(sigma, W, base = base, perturbation = perturbation, ncores = ncores, ...)
    } else {
      sigInv <- try(network(sigma, ...), silent = T)
      if (class(sigInv) == "try-error") {
        stop("Cannot use supplied network option?")
      }
    }
    
    ret_val <- get_Y_new_and_sigSum(N, Y_p, W, eY, base, sigInv, MCiter, MCburn, stepsize)
    Y_new <- ret_val$Y_new
    sigSum <- ret_val$sigSum
    
    # update b0, means across OTUs
    b0 <- colMeans(Y_new)
    
    sigma <- sigSum/(N * (MCiter - MCburn))
    
    # update b
    eY <- tcrossprod(rep(1, N), b0) 
    if (!no_covariates)  {
      b <- OLS(X_c, Y_new)
      eY <- eY + X_c %*% b
    }
    
    ### STORE after updating
    ## @Bryan TODO: can this be cleaned up?
    # accept_list[em,] <- accepts
    #b0_list <- rbind(b0_list, b0)
    b0_list[em + 1,] <- b0
    if (!no_covariates) {
      if (is.matrix(b)) {
        b_list[,, em + 1] <- b
      } else {
        b_list[, em + 1] <- b
      }
    }
    sigma_list[,, em + 1] <- sigma
    #end <- proc.time()
  }
  utils::setTxtProgressBar(pb, EMiter)
  cat("\n")
  
  ## Next: take average of EM samples past burn. Included init, so have EMiter+1 total 
  b0_EM <- colMeans(b0_list[(EMburn + 1):(EMiter + 1), ])
  
  output_list$beta0 <- b0_EM
  
  ### fitted value of Y
  if (no_covariates) {
    fitted_y <- matrix(b0_EM, ncol = Q-1, nrow = N, byrow=T)
  } else {
    b_list_reduced <- b_list[, , (EMburn + 1):(EMiter + 1)] 
    
    if (length(dim(b_list_reduced)) == 2) {
      b_list_reduced <- b_list_reduced %>% array(c(1, dim(b_list_reduced)))
    }
    
    b_EM <- apply(b_list_reduced, c(1, 2), mean)
    output_list$beta <- b_EM
    
    fitted_y <- X_c %*% b_EM + matrix(b0_EM, ncol = Q-1, nrow = N, byrow=T)
    
  }
  
  # burn-in sigma and take average; fine b/c class of positive def matrices is closed under addition
  sigma_em <- apply(sigma_list[, , (EMburn + 1):(EMiter + 1)], c(1, 2), mean)
  
  # store parameters
  output_list$sigma <- sigma_em
  output_list$fitted_y <- fitted_y
  output_list$fitted_z <- to_composition_matrix(fitted_y, base=base)
  output_list$base <- base
  output_list
}

diagonal_network <- function(sigma) {
  # take diagonal vec, make it a diagonal matrix, solve
  diag(1/diag(sigma))
}


default_network <- function(sigma) {
  test <- try(chol(sigma), silent = T)
  if (class(test) == "try-error") {
    test2 <- try(svd(sigma), silent = T)
    if (class(test2) == "try-error") {
      message("SVD failed; sigma is")
      print(sigma)
      stop()
    } else {
      sigInv <- MASS::ginv(sigma)
    }
  } else {
    sigInv <- chol2inv(chol(sigma))
  }
  sigInv
}

#' stars
#'
#' Estimate the network using the package SpiecEasi
#'
#' @param sigma current estimate of sigma
#' @param W corresponding count matrix
#' @param base OTU index used for base
#' @param perturbation size of purturbation used for to_log_ratios, defaults to 0.05
#' @param ncores number of cores to use, defaults to 1
#' @param ... other arguments to pass
#'
#'
stars <- function(sigma, W, base, perturbation, ncores, ...) {
  
  if (!requireNamespace("glasso", quietly = TRUE) | !requireNamespace("SpiecEasi", quietly = TRUE) |
      !requireNamespace("pulsar", quietly = TRUE) | !requireNamespace("huge", quietly = TRUE)) {
    stop("Packages glasso, pulsar and SpiecEasi are needed for this function to work. \n
         Please install them.",
         call. = FALSE)
  }
  
  Y_p <- to_log_ratios(W, base = base, perturbation = perturbation)
  print("pulsar::getMaxCov(Y_p)")
  print(pulsar::getMaxCov(Y_p))
  hugeargs <- list(lambda=seq(from=ifelse(is.nan(pulsar::getMaxCov(Y_p)), 50, pulsar::getMaxCov(Y_p)),
                              to = 0.01, length.out=5),
                   verbose=FALSE)
  
  out.p <- pulsar(Y_p, fun=huge::huge, fargs=hugeargs,
                  rep.num=10, criterion='stars', ncores=1)
  gl <- glasso::glasso(sigma, rho = lams[out.p$stars$opt.index])
  gl$wi
  
}

