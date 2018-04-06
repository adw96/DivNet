#' extractOTU
#'
#' This function extracts the desired combination of random and most prevalent OTUs, including 0.
#'
#' @param W raw count matrix, with OTUs as columns
#' @param nrand number of random OTUs desired
#' @param nmost number of most prevalent OTUs desired
#'
#' @export
extractOTU <- function(W, nrand, nmost) {
    if (nmost >= 1) {
        mostI <- whichpart(colSums(W), n = nmost)
        leftO <- seq(dim(W)[2])[-mostI]
        randI <- sample(leftO, nrand)
        return(W[, c(mostI, randI)])
    } else {
        randI <- sample(seq(dim(W)[2]), nrand)
        return(W[, c(randI)])
    }
}
