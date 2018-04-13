#' @export
make_w <- function(mu, Sigma, mm) {
  q <- ncol(Sigma)
  n <- nrow(mu)
  if (length(mm) == 1) mm <- rep(mm, n)
  
  Y <- matrix(NA, nrow=n, ncol=q)
  for (i in 1:n) {
    Y[i, ] <- MASS::mvrnorm(n=1, mu=mu[i, ], Sigma=Sigma)
  }
  exp_Y <- exp(Y)
  denominators <- (exp_Y %>% apply(1, sum) + 1)
  compositions <- cbind(exp_Y, 1) / denominators # n x (q+1)
  my_w <- matrix(NA, nrow=n, ncol=q+1)
  for (i in 1:n) {
    my_w[i, ] <- rmultinom(n=1, size=mm[i], prob=compositions[i,])
  }
  my_w
}