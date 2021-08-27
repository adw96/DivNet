#' Shannon Index
#'
#' @param proportions proportions
#'
#' @export
shannon_true <- function(proportions) {

  if (all(is.na(proportions)) | (sum(proportions) - 1)^2 < 1e-10) {
    input <- proportions[proportions > 0]
    output <- -sum(input*log(input, base=exp(1)))
  } else {
    stop("shannon_true needs a vector of proportions that sum to 1")
  }
  output
}
#' Simpson Index
#'
#' @param proportions proportions
#'
#' @export
simpson_true <- function(proportions) {

  if (all(is.na(proportions)) | (sum(proportions) - 1)^2 < 1e-8) {
    output <- sum(proportions^2)
  } else {
    stop("simpson_true needs a vector of proportions that sum to 1")
  }
  output
}



#' Bray-Curtis Distance (Fast)
#'
#' @param p1s First set of proportions
#' @param p2s Second set of proportions
#'
#' @export
bc_fast <- function(p1s, p2s,
                    sum_to_one = TRUE) {

  if ( (!sum_to_one ) | all(is.na(p1s)) | all(is.na(p2s)) | ((sum(p1s) - 1)^2 < 1e-8 & (sum(p2s) - 1)^2 < 1e-8)) {
    output <- 1 - (pmin(p1s, p2s) %>% sum)
  } else {
    stop("bc_fast needs a vector of proportions that sum to 1")
  }
  output
}


#' Euclidean Distance (Fast)
#'
#' @param p1s First set of proportions
#' @param p2s Second set of proportions
#'
#' @export
euc_fast <- function(p1s, p2s) {

  if (all(is.na(p1s)) | all(is.na(p2s)) | ((sum(p1s) - 1)^2 < 1e-8 & (sum(p2s) - 1)^2 < 1e-8)) {
    output <-  (p1s - p2s)^2 %>% sum %>% sqrt
  } else {
    stop("euc_fast needs a vector of proportions that sum to 1")
  }
  output
}


#' Bray-Curtis Distance (True)
#'
#' @param p1s First set of proportions
#' @param p2s Second set of proportions, defaults to NULL
#'
#' @export
bray_curtis_true <- function(p1s, p2s = NULL,
                             sum_to_one = TRUE) {
  if (p2s %>% is.null & !(dim(p1s)[1] %>% is.null)) {
    bc <- outer(1:dim(p1s)[1], 1:dim(p1s)[1], FUN = Vectorize( function(i,j) bc_fast(p1s[i,], p1s[j,],
                                                                                     sum_to_one) ))
  } else {
    if (!((p1s %>% length) == (p2s %>% length))) {
      stop("Attempted to compute the Bray-Curtis distance between communities of different sizes.")
    } else {
      bc <- bc_fast(p1s, p2s)
    }
  }
  bc
}



#' Euclidean Distance (True)
#'
#' @param p1s First set of proportions
#' @param p2s Second set of proportions, defaults to NULL
#'
#' @export
euclidean_true <- function(p1s, p2s = NULL) {
  if (!(dim(p1s) %>% is.null)) {
    euc <- outer(1:dim(p1s)[1], 1:dim(p1s)[1], FUN = Vectorize( function(i,j) euc_fast(p1s[i,], p1s[j,]) ))
  } else {
    if (!((p1s %>% length) == (p2s %>% length))) {
      stop("Attempted to compute the Euclidean distance between communities of different sizes.")
    } else {
      euc <- euc_fast(p1s, p2s)
    }
  }
  euc
}

