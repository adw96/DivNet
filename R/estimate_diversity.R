#' estimate_diversity
#' 
#' @param fitted_aitchison object produced by fit_aitchison
#' 
#' 
#' @export
estimate_diversity <- function(fitted_aitchison) {
  
  # compute estimate  
  zz <- fitted_aitchison$fitted_z

  output_list <- list()
  
  # Estimate Shannon diversity
  output_list[["shannon"]] <- apply(zz, 1, shannon_true)
  
  # Estimate Simpson
  output_list[["simpson"]] <- apply(zz, 1, simpson_true)
  
  # Estimate Bray Curtis
  output_list[["bray-curtis"]] <- bray_curtis_true(zz)
  
  # Estimate Euclidean
  output_list[["euclidean"]] <- euclidean_true(zz)
  
  output_list
}
