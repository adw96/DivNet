#' Simulate counts
#' 
#' @param mu Mean matrix. Rows are samples
#' @param Sigma Covariance matrix 
#' @param mm Sampling depth (= number of individuals observed)
#' @param base The base taxon for logratio transformations
#' 
#' @export
make_w <- function(mu, Sigma, mm, base = NULL) {
  Q <- ncol(Sigma)
  N <- nrow(mu)
  
  if (ncol(Sigma) != ncol(mu)) {
    stop("mu and Sigma are not the right dimensions")
  }
  # Global variable \fix
  i <- NULL
  if (length(mm) == 1) mm <- rep(mm, N)
  
  Y <- foreach(i = 1:N, .combine = rbind) %do% {
    MASS::mvrnorm(n = 1, mu = mu[i, ], Sigma)
  }
  
  compositions <- to_composition_matrix(Y, base = base) 
  
  my_w <- foreach(i = 1:N, .combine = rbind) %do% {
    stats::rmultinom(n = 1, size = mm[i], prob = compositions[i,]) %>% c
  }
  my_w
}
