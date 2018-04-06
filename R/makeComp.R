#' makeComp
#'
#' This function calculates the observed composition from raw counts
#'
#' @param W raw count matrix, with OTUs as columns
#'
#' @export
makeComp <- function(W) {
    return(W/rowSums(W))
}
