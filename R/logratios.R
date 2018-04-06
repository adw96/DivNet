#' logratios
#'
#' This function transforms from count data to logratios
#'
#' @param W raw count data, with OTUs as columns
#' @param base base OTU value
#' @param p how much to purturb zero counts, defaults to 0.05
#'
#' @export
logratios <- function(W, base, p = 0.05) {
    W <- as.matrix(W)
    # get purturbed Y, apply returns arguments as columns CHECK: Apply forces transpose. Worth it?
    Y.purt <- t(apply(W, 1, getPurt, base = base, p = p))
    attr(Y.purt, "center") = apply(Y.purt, 2, mean)
    return(Y.purt)
}
