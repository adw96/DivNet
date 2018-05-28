#' #' to_composition
#' #'
#' #' This function is currently not used by any function in DivNet and is slated for removal
#' 
#' #' This function transforms a vector of logratios Y to a composition matrix X
#' #'
#' #' @param Y A vector of logratios
#' #' @param base base taxon used to calculate logratios
#' #' 
#' #' @importFrom R.utils insert
#' #'
#' #' @export
#' to_composition <- function(Y, base = NULL) {
#'   
#'   stopifnot(!is.array(Y))
#'   
#'   if (is.null(base)) base <- (length(Y) + 1)
#'   
#'   exp_Y <- exp(Y)
#'   denominator <- sum(exp_Y) + 1
#'   # c(exp_Y, 1) / denominator
#'   
#'   numerators <- R.utils::insert(exp_Y, ats = base, values = 1)
#'   numerators / denominator
#' }

#' to_composition_matrix
#'
#' This function transforms a logratio matrix Y to a composition matrix X
#'
#' @param Y A matrix of logratios. The samples are listed across the rows and the taxa are listed across the columns 
#' @param base base taxon used to calculate logratios
#'
#' @export
to_composition_matrix <- function(Y, base = NULL) {
  if (is.null(base)) base <- (ncol(Y) + 1)
  
  
  out <- matrix(0, nrow = nrow(Y), ncol = ncol(Y) + 1)
  out[,-base] <- exp(Y)
  out[,base] <- 1
  out / rowSums(out)
}

#' to_log_ratios
#'
#' This function transforms from count data to logratios
#'
#' @param W raw count data, with OTUs as columns
#' @param base base OTU value
#' @param perturbation how much to purturb zero counts, defaults to 0.05
#'
#' @export
to_log_ratios <- function(W, base, perturbation = 0.05) {
  stopifnot(is.matrix(W))
  
  tmp <- pmax(W, perturbation)
  Ztmp <- tmp/rowSums(tmp)
  Y_purterbed <- log(Ztmp[,-base]/Ztmp[,base])
  return(as.matrix(Y_purterbed))
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
