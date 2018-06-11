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
    error("to_log_ratios needs which taxon is base taxon")
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


#' #' to_counts
#' #'
#' #' This function is currently not used by any function in DivNet and is slated for removal
#' #' This function transforms logratio matrix Y to counts W
#' #'
#' #' @param Y logratio matrix
#' #' @param M vector of observed counts by row
#' #' @param base base value used to calculate logratios
#' #'
#' #' @export
#' to_counts <- function(Y, M, base) {
#'   N <- nrow(Y)
#'   Q <- ncol(Y) + 1
#'   exp_Y <- exp(Y)
#'   sum_exp_Y <- rowSums(exp_Y)
#'   X <- matrix(0, N, Q)
#'   X[,-base] <- exp_Y/(sum_exp_Y + 1)
#'   X[,base] <- 1 - rowSums(X)
#'   W <- matrix(rmultinom(1, M, prob = X), nrow = N, ncol = Q)
#'   ## @Bryan TODO: can this be cleaned up?
#'   
#'   # for (i in 1:N) {
#'   #   X[i, -base] <- exp_Y[i, ]/(sum_exp_Y[i] + 1)
#'   #   X[i, base] <- 1/(sum_exp_Y[i] + 1)
#'   #   W[i, ] <- rmultinom(n = 1, size = M[i], prob = X[i, ])
#'   # }
#'   return(W)
#' }

#' #' make_composition
#' #'
#' #' This function is currently not used by any function in DivNet and is slated for removal
#' #' This function calculates the observed composition from raw counts
#' #'
#' #' @param W raw count matrix, with OTUs as columns
#' #'
#' #' @export
#' make_composition <- function(W) {
#'   return(W/rowSums(W))
#' }
