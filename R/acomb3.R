#' acomb3
#'
#' This function works like cbind, but along the third dimension of an array
#'
#' @param ... Array
#' @export
acomb3 <- function(...) abind(..., along = 3)
