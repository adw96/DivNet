#' getPurt
#'
#' This function purturbates a row of the raw count matrix. Likely will not be used directly.
#'
#' @param Wi raw count data row
#' @param base base OTU value
#' @param p how much to purturb zero counts, defaults to 0.05
#'
#' @export
getPurt <- function(Wi, base, p = 0.05) {
    zeros <- which(Wi == 0)
    nz <- length(zeros)
    Q <- length(Wi)
    Z.purt <- rep(0, Q)
    if (nz != 0) {
        Z.purt[zeros] <- (Wi[zeros] + p)/(sum(Wi) + p * nz)
        Z.purt[-zeros] <- (Wi[-zeros])/(sum(Wi) + p * nz)
    } else {
        Z.purt <- Wi/sum(Wi)
    }
    
    # Return Y.purt
    return(log(Z.purt[-base]/Z.purt[base]))
}
