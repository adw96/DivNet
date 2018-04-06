#' Xsim
#'
#' Function to simulate raw compositions X from model fit
#'
#' @param out model fit from LNM.EM
#' @param W raw count matrix used in model fit
#' @param niter number of simulations, defaults to 1000
#' @param fast whether or not to simulate from diagonal sigma with rmvn
#'
#' @export
Xsim <- function(out, W, X = NULL, niter = 1000, fast = FALSE) {
    N <- nrow(out$Y)
    Q <- ncol(out$Y) + 1
    base <- out$base
    X.m <- array(0, dim = c(N, Q, niter))
    M <- apply(W, 1, sum)
    if (is.null(X)) {
        mu <- get_mu(out)
    } else {
        mu <- get_mu(out, X)
    }
    # if can sample from the same mu every time (no covariates)
    if (is.vector(mu)) {
        for (i in 1:niter) {
            # set up Y take mu, sigma, simulate N Y_i's
          if (fast == TRUE) {
            Y.m <- rmvn(n = N, mu = mu, sigma = diag(diag(out$sigma)))
          } else {
            Y.m <- mvrnorm(n = N, mu = mu, Sigma = out$sigma)
          }
            W.m <- YtoW(Y = Y.m, M = M, base = base)
            X.m[, , i] <- makeComp(W.m)
        }
    }
    if (is.matrix(mu)) {
        for (i in 1:niter) {
          if (fast == TRUE) {
            # apply out as vector stores as columns, transpose
            Y.m <- t(apply(mu, 1, function(x) rmvn(n = 1, mu = x, sigma = diag(diag(out$sigma)))))
          } else {
            Y.m <- t(apply(mu, 1, function(x) mvrnorm(n = 1, mu = x, sigma = out$sigma)))
          }
            W.m <- YtoW(Y = Y.m, M = M, base = base)
            X.m[, , i] <- makeComp(W.m)
        }
    }
    return(X.m)
}
