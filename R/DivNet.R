#' DivNet package documentation.
#'
#' This package provides various methods and tools to analyze microbiome data
#'
#' A complete description of all package features, along
#' with examples, can be found at \url{https://github.com/adw96/DivNet}.
#'
#' @import stats
#' @importFrom MASS mvrnorm ginv
#' @importFrom utils globalVariables
#' @importFrom mvnfast rmvn
#'
#' @useDynLib DivNet
#' @importFrom Rcpp evalCpp
#' @import RcppEigen
#' 
#' @import doParallel
#' @import abind
#' @import foreach
#' 
#' @name DivNet-package
#' @aliases DivNet
NULL


## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))