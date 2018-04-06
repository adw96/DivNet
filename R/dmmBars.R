#' dmmBars
#'
#' Function to plot DMM model fit, with empirical 95% bars
#'
#' @param out model fit from package DirichletMultinomial
#' @param W raw count matrix used in model fit
#' @param main desired title for plot, defaults to DMM Model Fit
#'
#' @export
dmmBars <- function(out, W, main = "DMM Model Fit") {
    mu <- out@fit$Estimate/sum(out@fit$Estimate)
    eZ <- matrix(mu, byrow = TRUE, nrow = nrow(W), ncol = ncol(W))
    Z <- makeComp(W)
    par(mar = c(9, 6, 4, 3) + 0.1, cex.main = 2)
    plot(as.vector(Z), as.vector(eZ), xlab = "", ylab = "", main = main, pch = 20)
    # abline(a=0,b=1)
    mtext(expression(frac(W[ij], sum(W[ij], j))), side = 1, line = 7.5, cex = 2)
    mtext("Scaled Dirichlet Component Estimates", side = 2, line = 3, cex = 1.5)
    eBarY <- eZ[1, ]
    eBarX1 <- out@fit$Lower/sum(out@fit$Lower)
    eBarX2 <- out@fit$Upper/sum(out@fit$Upper)
    arrows(eBarX1, eBarY, eBarX2, eBarY, col = "red", code = 3, angle = 90, length = 0.01, lwd = 2)
}
#' dmmBarsT
#'
#' Function to plot DMM model fit, with empirical 95% bars, transposed
#'
#' @param out model fit from package DirichletMultinomial
#' @param W raw count matrix used in model fit
#' @param main desired title for plot, defaults to DMM Model Fit
#'
#' @export
dmmBarsT <- function(out, W, main = "DMM Model Fit") {
    eZ <- makeComp(W)
    Q <- ncol(W)
    N <- nrow(W)
    
    Z <- rep(1:Q, each = N)
    par(mar = c(4.5, 9, 4, 3) + 0.1, cex.main = 2)
    plot(as.vector(Z), as.vector(eZ), xlab = "", ylab = "", main = main, pch = 20)
    # abline(a=0,b=1)
    mtext("OTU Index", side = 1, line = 3, cex = 2)
    mtext(expression(frac(W[ij], sum(W[ij], j))), side = 2, las = 1, line = 3, cex = 2)
    eBarX <- 1:Q
    eBarY1 <- out@fit$Lower/sum(out@fit$Lower)
    eBarY2 <- out@fit$Upper/sum(out@fit$Upper)
    arrows(eBarX, eBarY1, eBarX, eBarY2, col = "red", code = 3, angle = 90, length = 0.01, lwd = 2)
}
