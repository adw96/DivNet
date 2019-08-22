diagonal_network <- function(sigma) {
  # take diagonal vec, make it a diagonal matrix, solve
  diag(1/diag(sigma))
}


default_network <- function(sigma) {
  # If possible, use the Cholesky decomp as it is faster than svd
  test <- try(chol(sigma), silent = T)
  if (class(test) == "try-error") { # cholesky decomp failed
    # If svd will fail, then MASS::ginv will fail as it calls svd.  So no need to test svd first.
    test2 <- try(sigInv <- MASS::ginv(sigma))
    if (class(test2) == "try-error") {
      message("MASS::ginv failed (probably because svd failed); sigma is")
      print(sigma)
      stop()
    }
  } else {
    sigInv <- chol2inv(chol(sigma))
  }
  sigInv
}

#' stars
#'
#' @param sigma current estimate of sigma
#' @param W corresponding count matrix
#' @param base OTU index used for base
#' @param perturbation size of purturbation used for to_log_ratios, defaults to 0.05
#' @param ncores number of cores to use, defaults to 1
#' @param ... other arguments to pass
#'
#' Estimate the network using the package SpiecEasi
stars <- function(sigma, W, base, perturbation, ncores, ...) {
  
  if (!requireNamespace("glasso", quietly = TRUE) | !requireNamespace("SpiecEasi", quietly = TRUE)) {
    stop("Packages glasso and SpiecEasi are needed for this function to work. \n
         Please install them.",
         call. = FALSE)
  }
  Y_p <- to_log_ratios(W, base = base, perturbation = perturbation)
  inverse_covariance_estimate <- SpiecEasi::sparseiCov(data=Y_p, method = "glasso")
  selected <- SpiecEasi::icov.select(inverse_covariance_estimate, rep.num=5, ncores = ncores)
  gl <- glasso::glasso(sigma, rho = selected$opt.lambda)
  gl$wi
}
