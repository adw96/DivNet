#' toComposition
#'
#' This function transforms a vector of logratios Y to a composition matrix X
#'
#' @param Y A vector of logratios
#' @param base base taxon used to calculate logratios
#' 
#' @importFrom R.utils insert
#'
#' @export
toComposition <- function(Y, base = NULL) {
  
  stopifnot(!is.array(Y))
  
  if (is.null(base)) base <- (length(Y) + 1)
  
  exp_Y <- exp(Y)
  denominator <- sum(exp_Y) + 1
  # c(exp_Y, 1) / denominator
  
  numerators <- R.utils::insert(exp_Y, ats = base, values = 1)
  numerators / denominator
}

#' toCompositionMatrix
#'
#' This function transforms a logratio matrix Y to a composition matrix X
#'
#' @param Y A matrix of logratios. The samples are listed across the rows and the taxa are listed across the columns 
#' @param base base taxon used to calculate logratios
#'
#' @export
toCompositionMatrix <- function(Y, base = NULL) {
  
  apply(X=Y, MARGIN=1, FUN=toComposition, base = base) %>% t
  
}

#' toLogRatios
#'
#' This function transforms from count data to logratios
#'
#' @param W raw count data, with OTUs as columns
#' @param base base OTU value
#' @param perturbation how much to purturb zero counts, defaults to 0.05
#'
#' @export
toLogRatios <- function(W, base, perturbation = 0.05) {
  W <- as.matrix(W)
  # get purturbed Y, apply returns arguments as columns CHECK: Apply forces transpose. Worth it?
  Y.purt <- t(apply(W, 1, getPurt, base = base, perturbation = perturbation))
  attr(Y.purt, "center") = apply(Y.purt, 2, mean)
  return(Y.purt)
}


#' toCounts
#'
#' This function transforms logratio matrix Y to counts W
#'
#' @param Y logratio matrix
#' @param M vector of observed counts by row
#' @param base base value used to calculate logratios
#'
#' @export
toCounts <- function(Y, M, base) {
  N <- nrow(Y)
  Q <- ncol(Y) + 1
  exp_Y <- exp(Y)
  sum_exp_Y <- apply(exp_Y, 1, sum)
  X <- W <- matrix(0, N, Q)
  
  ## @Bryan TODO: can this be cleaned up?
  for (i in 1:N) {
    X[i, -base] <- exp_Y[i, ]/(sum_exp_Y[i] + 1)
    X[i, base] <- 1/(sum_exp_Y[i] + 1)
    W[i, ] <- rmultinom(n = 1, size = M[i], prob = X[i, ])
  }
  return(W)
}

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
