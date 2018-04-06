#' MCrow
#'
#' This function simulates MC step for a single row. Should not need to be used by user directly.
#'
#' @param Yi row of logratio matrix
#' @param Wi corresponding row of count matrix
#' @param eYi current expected value of logratio matrix
#' @param Q number of OTUs minus the base, or length of Yi
#' @param base OTU index used for base
#' @param sigInv current estimate of sigma inverse
#' @param MCiter number of MC samples to generate
#' @param stepsize variance used for MH samples, defaults to 1. Tweak to adjust acceptance ratio
#'
#' @export
MCrow <- function(Yi, Wi, eYi, Q, base, sigInv, MCiter, stepsize = 1) {
    # extra column for acceptance indicator
    Yi.MH <- matrix(0, MCiter, Q)
    for (i in 1:MCiter) {
        # cat(i,'\n') proposal
        Yi.star <- Yi + rnorm(Q - 1, 0, stepsize)
        # denominator
        Eq5pt1 <- sum(Wi) * (log(sum(exp(Yi)) + 1) - log(sum(exp(Yi.star)) + 1))
        # numerator
        Eq5pt2 <- sum(Wi[-base] * (Yi.star - Yi))
        Eq5pt3 <- -0.5 * crossprod((Yi.star - eYi), sigInv) %*% (Yi.star - eYi)
        Eq5pt4 <- -0.5 * crossprod((Yi - eYi), sigInv) %*% (Yi - eYi)
        fullRat <- Eq5pt1 + Eq5pt2 + Eq5pt3 - Eq5pt4
        acceptance <- min(1, exp(fullRat))
        temp <- runif(1)
        aVal <- 0
        if (temp < acceptance | is.nan(acceptance)) {
            Yi <- Yi.star
            aVal <- 1
        }
        # each row is a resampling of Yi, first col is whether accepted or not
        Yi.MH[i, ] <- c(aVal, Yi)
    }
    
    ## returns a matrix Yi.MH: MH samples of row, first column ind if accepted (MCiter rows, Q cols)
    return(Yi.MH)
}
