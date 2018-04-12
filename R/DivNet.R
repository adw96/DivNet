#' DivNet package documentation.
#'
#' This package provides various methods and tools to analyze microbiome data
#'
#' A complete description of all package features, along
#' with examples, can be found at \url{https://github.com/adw96/DivNet}.
#'
#' @import stats
#' @import readr
#' @import extrafont
#' @import fontcm
#' @importFrom MASS mvrnorm ginv
#' @importFrom parallel detectCores
#' @import doParallel
#' @import abind
#' @import doSNOW
#' @import PDSCE
#' @import grDevices
#' @import graphics
#' @import foreach
#' @importFrom utils globalVariables
#' @importFrom mvnfast rmvn
#'
#' @name DivNet-package
#' @aliases DivNet
NULL

## Quiets global binding variables warning when using 'i' in foreach
if (getRversion() >= "2.15.1") utils::globalVariables(c("i"))

#' @export
divnet <-  function(W, 
                    X = NULL, 
                    tuning = NULL,
                    perturbation = NULL, 
                    network = NULL,
                    base = NULL,
                    ncores = NULL,
                    variance = "parametric",
                    B = 5,
                    nsub = NULL,
                    ...) {
  
  fitted_aitchison <- fit_aitchison(W, 
                                    X = X, 
                                    tuning = tuning,
                                    perturbation = perturbation, 
                                    network = network,
                                    base = base,
                                    ncores = ncores,
                                ...)
  zz <- fitted_aitchison$fitted_z
  output_list <- get_diversities(zz)
  
  # @Amy TODO implement variance estimates
  # @Amy TODO parallelise
  if (variance == "parametric") {
    
    # resample from models
    parametric_list <- replicate(B, 
                                 parametric_variance(fitted_aitchison, 
                                                     W = W,
                                                     X = X, 
                                                     tuning = tuning,
                                                     perturbation = perturbation, 
                                                     network = network,
                                                     base = base,
                                                     ncores = ncores,
                                                     ...), 
                                 simplify=F)
    
    variance_estimates <- get_diversity_variance(parametric_list)
    output_list <- c(output_list, variance_estimates)
  } else if (variance == "nonparametric") {
    if (is.null(nsub)) nsub <- ceiling(dim(W)[1]/2)
    nonparametric_list <- replicate(B, 
                                    nonparametric_variance(W = W,
                                                           X = X, 
                                                           tuning = tuning,
                                                           perturbation = perturbation, 
                                                           network = network,
                                                           base = base,
                                                           ncores = ncores,
                                                            ...), 
                                    simplify=F)
    variance_estimates <- get_diversity_variance(nonparametric_list)
    output_list <- c(output_list, variance_estimates)
  }
  output_list
}