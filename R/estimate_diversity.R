#' @export
estimate_diversity <- function(fitted_aitchison, errors = "delta") {
  
  
  # Get estimated proportions
  X <- fitted_aitchison$X
  sigma <- fitted_aitchison$sigma
  fitted_y <- fitted_aitchison$fitted_y
  N <- dim(fitted_y)[1]
  Q <- dim(fitted_y)[2] + 1
  
  # compute estimate  
  exp_b <- fitted_y %>% exp
  denominator <- 1 + apply(exp_b, 1, sum)
  zz <- cbind(exp_b, 1) / matrix(denominator, nrow=N, ncol = Q, byrow=F)
  
  output_list <- list()
  
  # Estimate Shannon diversity
  shannon_estimates <- zz %>% apply(1, breakaway::shannon)
  output_list[["shannon"]] <- shannon_estimates
  
  # Estimate Simpson
  simpson_estimates <- zz %>% apply(1, breakaway::simpson)
  output_list[["simpson"]] <- simpson_estimates
  
  
  # Estimate Bray Curtis
  output_list[["bray-curtis"]] <- bray_curtis(zz)
  
  # Estimate Euclidean
  output_list[["euclidean"]] <- euclidean(zz)
  
  
  empirical_delta_g <- log(exp_b / denominator) - log ( 1 / denominator)
  variance_shannon <- empirical_delta_g %*% fitted_aitchison$variance_z %*% t(empirical_delta_g)
  output_list[["shannon_variance"]] <- sqrt(diag(variance_shannon))
  
  empirical_delta_g <- 2 * (exp_b / denominator) - 2 * (1 / denominator)
  variance_simpson <- empirical_delta_g %*% fitted_aitchison$variance_z %*% t(empirical_delta_g)
  output_list[["simpson_variance"]] <- sqrt(diag(variance_simpson))
  
  output_list
}