#' get_mu
#'
#' Function to get mu from model output
#'
#' @param out model fit from LNM.EM or LNM.EM.nocov
#'
#' @export
get_mu <- function(out, X = NULL) {
    N <- nrow(out$Y)
    if (!is.null(out$b)) {
        return(X %*% out$b + tcrossprod(rep(1, N), out$b0))
    } else {
        return(out$b0)
    }
}
