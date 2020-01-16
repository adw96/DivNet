#' MCrow
#'
#' This function simulates MC step for a single row. Should not need to be used by user directly.
#'
#' @author Bryan Martin
#'
#' @param Yi row of logratio matrix
#' @param Wi corresponding row of count matrix
#' @param eYi current expected value of logratio matrix
#' @param Q number of OTUs, or length of Wi
#' @param base OTU index used for base
#' @param sigInv current estimate of sigma inverse
#' @param MCiter number of MC samples to generate
#' @param stepsize variance used for MH samples, defaults to 1. Tweak to adjust acceptance ratio
MCrow <- function(Yi, Wi, eYi, Q, base, sigInv, MCiter, stepsize = 1) {
  # extra column for acceptance indicator
  Yi.MH <- matrix(0, MCiter, Q)
  
  # Precompute all uniform dist. vals we will need.
  unif_vals <- runif(MCiter)
  
  # Calculate this once for speed.
  Wi_no_base <- Wi[-base]
  
  Yi.MH <- mcrow_mc_iteration(MCiter,
                              stepsize,
                              unif_vals,
                              Wi,
                              Wi_no_base,
                              Yi,
                              eYi,
                              sigInv)
  
  ## returns a matrix Yi.MH: MH samples of row, first column ind if accepted (MCiter rows, Q cols)
  return(Yi.MH)
}