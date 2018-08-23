#' to_composition_matrix
#'
#' This function transforms a logratio matrix Y to a composition matrix X
#'
#' @param Y A vector or matrix of logratios. If a matrix, the samples are listed across the rows and the taxa are listed across the columns 
#' @param base base taxon used to calculate logratios
#'
#' @export
to_composition_matrix <- function(Y, base = NULL) {
  
  if (class(Y) == "list")  {
    stop("to_composition_matrix was passed a list and doesn't know what to do")
  }
  
  if (is.null(dim(Y))) { # if a vector
    qq <- length(Y)
    nn <- 1
  } else if (length(dim(Y)) == 2) {
    qq <- ncol(Y)
    nn <- nrow(Y)
  } else {
    stop("to_composition_matrix generated an error")
  }

  if (is.null(base)) {
    warning("base is NULL, setting to last taxon")
    base <- qq + 1
  }
  
  out <- matrix(nrow = nn, ncol = qq + 1)
  out[,-base] <- exp(Y)
  out[,base] <- 1
  out / rowSums(out)
}

#' to_log_ratios
#'
#' This function transforms from original scale to logratios
#'
#' @param W Any compositional (count or proportion) data to be transformed, with OTUs as columns
#' @param base base OTU value
#' @param perturbation how much to purturb zero counts, defaults to 0.05
#'
#' @export
to_log_ratios <- function(W, base, perturbation = 0.05) {
  
  if (is.null(base)) {
    stop("to_log_ratios needs which taxon is base taxon")
  }
  
  if (is.null(dim(W))) { # if a vector
    W <- matrix(W, nrow = 1)
  } 
  
  if (all(W %% 1 == 0)) {
    W <- pmax(W, perturbation)
    W <- W/rowSums(W)
  }
  
  Y_purterbed <- log(W[,-base]/W[,base])
  return(Y_purterbed)
}

