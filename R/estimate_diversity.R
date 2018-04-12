#' estimate_diversity
#' 
#' @param fitted_aitchison object produced by fit_aitchison
#' @param variance method to get variance of estimates. Current options are "parametric" for parametric bootstrap, "nonparametric" for nonparametric bootstrap, and "none" for no variance estimates
#' 
#' @author Amy Willis
#' 
#' @export
estimate_diversity <- function(fitted_aitchison, W, variance = "parametric", B = 5, nsub = NULL, ...) {
  
  # compute estimate  
  zz <- fitted_aitchison$fitted_z
  output_list <- get_diversities(zz)
  
  # @Amy TODO implement variance estimates
  # @Amy TODO parallelise
  if (variance == "parametric") {
    
    # resample from models
    parametric_list <- replicate(B, 
                                 parametric_variance(fitted_aitchison, 
                                                     W = W,
                                                     X = fitted_aitchison$X, ...), 
                                 simplify=F)
    
    variance_estimates <- get_diversity_variance(parametric_list)
    output_list <- c(output_list, variance_estimates)
  } else if (variance == "nonparametric") {
    if (is.null(nsub)) nsub <- ceiling(dim(W)[1]/2)
    nonparametric_list <- replicate(B, 
                                    nonparametric_variance(W = W,
                                                           X = fitted_aitchison$X, 
                                                           nsub = nsub, ...), 
                                    simplify=F)
    variance_estimates <- get_diversity_variance(nonparametric_list)
    output_list <- c(output_list, variance_estimates)
  }
  
  output_list
}

get_diversities <- function(zz) {
  
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

get_diversity_variance <- function(list_of_fitted_models) {
  output_list <- list()
  
  output_list[["shannon-variance"]]  <- list_of_fitted_models %>% 
    lapply(function(x) x$shannon) %>% 
    simplify2array %>% 
    apply(1, var)
  
  output_list[["simpson-variance"]]  <- list_of_fitted_models %>% 
    lapply(function(x) x$simpson) %>% 
    simplify2array %>% 
    apply(1, var)
  
  output_list[["bray-curtis-variance"]] <- list_of_fitted_models %>% 
    lapply(function(x) x[["bray-curtis"]]) %>% 
    simplify2array %>% 
    apply(1:2, var)
  
  output_list[["euclidean-variance"]] <- list_of_fitted_models %>% 
    lapply(function(x) x$euclidean) %>% 
    simplify2array %>% 
    apply(1:2, var)
  
  output_list
}

### nonparametric bootstrap
nonparametric_variance <- function(W, X, nsub, ...) {
  N <- dim(W)[1]
  curly_b <- sample(1:N, replace=T, size = nsub)
  fitted_model <- fit_aitchison(W[curly_b, ],  X[curly_b, ], ...)
  
  eY <- fitted_model$beta0 + fitted_model$X %*% fitted_model$beta
  
  dots <- list(...)
  base <- ifelse(is.null(dots$base), 
                 which(rank(apply(W, 2, sum)) == max(rank(apply(W, 2, sum))))[1], 
                 dots$base)
  toCompositionMatrix(Y=eY, base=base) %>% get_diversities
  
}


### parametric bootstrap
parametric_variance <- function(fitted_aitchison, W, X, ...) {
  
  # unfortunately in this model we condition on M_i, so no randomising here
  ms <- apply(W, 1, sum)
  
  mw <- make_w(mu=fitted_aitchison$fitted_y, 
               Sigma = fitted_aitchison$sigma, 
               mm=ms)
  
  fitted_model <- fit_aitchison(mw,  X, ...)
  get_diversities(fitted_model$fitted_z)
}
