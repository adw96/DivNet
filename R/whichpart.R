#' whichpart
#'
#' This function selects the n most prevalent OTUs
#'
#' @param x OTU data
#' @param n number of OTUs to extract, defaults to 75
#'
#' @export
whichpart <- function(x, n = 75) {
    nx <- length(x)
    p <- nx - n
    xp <- sort(x, partial = p)[p]
    return(which(x > xp))
}
