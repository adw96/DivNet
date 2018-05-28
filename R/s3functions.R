

#' Print function 
#' 
#' @param x An object of class diversityEstimates
#' @param ... other arguments to be passed to print
#' @return NULL
#' 
#' @export
print.diversityEstimates <- function(x, ...) {
  dv <- x
  cat("An object of class `diversityEstimates` with the following elements:\n")
  sapply(1:length(names(dv)), function(i) { cat("  - ", names(dv)[i], "\n")})
  cat("Access individual components with, e.g., object$shannon and object$`shannon-variance`\n")
  cat("Use function testDiversity() to test hypotheses about diversity")
}



#' Plot function
#' 
#' TODO make more like the phyloseq plot richness
#' 
#' @param x An object of class diversityEstimates
#' @param ... other arguments to be passed to plot
#' @return An object of class ggplot
#' @export
plot.diversityEstimates <- function(x, ...) {
  dv <- x
  args <- match.call(expand.dots = TRUE)
  if (is.null(args$xx)) {
    args$xx <- "samples"
  }
  if (is.null(args$h0)) {
    args$h0 <- "shannon"
  }
  xx <- args$xx
  h0 <- args$h0
  
  lci <- dv[[h0]] - 2*sqrt(dv[[paste(h0, "-variance", sep = "")]])
  uci <- dv[[h0]] + 2*sqrt(dv[[paste(h0, "-variance", sep = "")]])
  df <- data.frame("names" = names(dv[[h0]]), 
                   "h0" = dv[[h0]], lci, uci, dv$X)
  df$names <- factor(df$names, levels = df$names)
  
  ggplot2::ggplot(df, ggplot2::aes(x = names, xend = names)) +
    ggplot2::geom_point(ggplot2::aes(x = names, y = h0)) +
    ggplot2::geom_segment(ggplot2::aes(y = lci, yend = uci)) +
    ggplot2::ylab(paste(h0, "estimate")) +
    ggplot2::xlab(xx) +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Test diversity
#' 
#' Hypothesis testing for alpha-diversity. 
#' 
#' 
#' @references Willis, A., Bunge, J., and Whitman, T. (2017). Improved detection of changes in species richness in high diversity microbial communities. \emph{JRSS-C.} 
#' 
#' @param dv An object of class diversityEstimates. The variable `X` used for the construction
#' @param h0 The alpha-diversity index to be tested for equality
#' @return A data frame similar to the output of `lm`
#' 
#' @export
testDiversity <- function(dv, h0 = "shannon") {
  cat("Hypothesis testing:\n")
  bt <- breakaway::betta(dv[[h0]], dv[[paste(h0, "-variance", sep="")]], X = dv[["X"]])
  cat(paste("  p-value for global test:", bt$global[2], "\n"))
  bt$table
}