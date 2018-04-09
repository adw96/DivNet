#' MicrobiomePack package documentation.
#'
#' This package provides various methods and tools to analyze microbiome data
#'
#' A complete description of all package features, along
#' with examples, can be found at \url{https://github.com/martinbryan/MicrobiomePack}.
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
#' @name MicrobiomePack-package
#' @aliases MicrobiomePack
NULL

## Quiets global binding variables warning when using 'i' in foreach
if (getRversion() >= "2.15.1") utils::globalVariables(c("i"))
