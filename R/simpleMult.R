#' simpleMult
#'
#' Function to plot simple multimial model fit, with error bars
#'
#' @param W observed count matrix
#' @param main desired title for plot, defaults to Simple Multinomial Variance
#'
#' @export
simpleMult <- function(W, main = "Simple Multinomial Variance") {
    W <- as.matrix(W)
    # W_ij / \sum_j W_{ij} summing all column entries together? That's a rowSum
    Xaxis <- W/rowSums(W)
    # phat_j = \sum_i W_ij / \sum_ij W_ij
    phatj <- colSums(W)/sum(W)
    Yaxis <- matrix(phatj, nrow = nrow(W), ncol = ncol(W), byrow = TRUE)
    par(mar = c(9, 9, 4, 3) + 0.1, cex.main = 2)
    plot(Xaxis, Yaxis, pch = 20, xlab = "", ylab = "", main = main)
    mtext(expression(frac(W[ij], sum(W[ij], j))), side = 1, line = 8, cex = 2)
    mtext(expression(frac(sum(W[ij], i), sum(W[ij], ij))), side = 2, las = 1, line = 3, cex = 2)
    xx <- seq(0, max(makeComp(W)) * 1.1, length = 100)
    lyy <- xx - 2 * sqrt(xx * (1 - xx)/sum(W))
    uyy <- xx + 2 * sqrt(xx * (1 - xx)/sum(W))
    lines(xx, lyy, col = "red")
    lines(xx, uyy, col = "red")
}
