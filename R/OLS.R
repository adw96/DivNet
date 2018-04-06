#' OLS function
#'
#' Function to simulate raw counts W from model fit
#'
#' @param X explanatory variables
#' @param Y response variables
#'
#' @export
OLS <- function(X, Y) {
  p <- ncol(X)
  centY <- scale(Y, center=TRUE, scale = FALSE)
  # quick function to aid apply
  aFun <- function(ycol) {
    return(tcrossprod(ginv(crossprod(X)), X) %*% ycol)
  }
  b <- apply(centY, 2, aFun)
  matrix(b, nrow = p)
}
