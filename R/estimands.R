#' Shannon Index 
#' 
#' @param proportions proportions
#' 
#' @export
shannon_true <- function(proportions) {
  proportions <- proportions/sum(proportions)
  -sum(proportions * log(proportions))
}

#' Simpson Index 
#' 
#' @param proportions proportions
#' 
#' @export
simpson_true <- function(proportions) {
  proportions <- proportions/sum(proportions)
  sum(proportions^2)
}



#' Bray-Curtis Distance (Fast)
#' 
#' @param p1s First set of proportions
#' @param p2s Second set of proportions
#' 
#' @export
bc_fast <- function(p1s, p2s) {
  p1s <- p1s/sum(p1s)
  p2s <- p2s/sum(p2s)
  1 - (pmin(p1s, p2s) %>% sum)
}


#' Euclidean Distance (Fast)
#' 
#' @param p1s First set of proportions
#' @param p2s Second set of proportions
#' 
#' @export
euc_fast <- function(p1s, p2s) {
  p1s <- p1s/sum(p1s)
  p2s <- p2s/sum(p2s)
  (p1s - p2s)^2 %>% sum %>% sqrt
}


#' Bray-Curtis Distance (True)
#' 
#' @param p1s First set of proportions
#' @param p2s Second set of proportions, defaults to NULL
#' 
#' @export
bray_curtis_true <- function(p1s, p2s = NULL) {
  if (p2s %>% is.null & !(dim(p1s)[1] %>% is.null)) {
    bc <- outer(1:dim(p1s)[1], 1:dim(p1s)[1], FUN = Vectorize( function(i,j) bc_fast(p1s[i,], p1s[j,]) ))
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

#' Calculate Alpha
#' 
#' @param mu Amy TODO
#' @param my_function Amy TODO
#' 
#' @export
calculate_alpha <- function(mu, my_function) {
  apply(mu %>% toCompositionMatrix, 1, my_function)
}

