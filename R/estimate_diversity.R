
#' divnet
#' 
#' @param fitted_aitchison object produced by fit_aitchison
#' @param variance method to get variance of estimates. Current options are "parametric" for parametric bootstrap, "nonparametric" for nonparametric bootstrap, and "none" for no variance estimates
#' 
#' @author Amy Willis
#' 
#' @export
divnet <-  function(W, 
                    X = NULL, 
                    fitted_model = NULL,
                    tuning = NULL,
                    perturbation = NULL, 
                    network = NULL,
                    base = NULL,
                    ncores = NULL,
                    variance = "parametric",
                    B = 5,
                    nsub = NULL,
                    ...) {
  
  if (is.null(fitted_model)) {
    fitted_model <- fit_aitchison(W, 
                                  X = X, 
                                  tuning = tuning,
                                  perturbation = perturbation, 
                                  network = network,
                                  base = base,
                                  ncores = ncores,
                                  ...)
  }
  zz <- fitted_model$fitted_z
  output_list <- get_diversities(zz)
  
  # @Amy TODO parallelise
  if (variance == "parametric") {
    
    # resample from models
    parametric_list <- replicate(B, 
                                 parametric_variance(fitted_model, 
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
                                                           nsub = nsub,
                                                           ...), 
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
  
  # na.rm = TRUE because nonparametric bootstrap works by subsampling
  output_list[["shannon-variance"]]  <- list_of_fitted_models %>% 
    lapply(function(x) x$shannon) %>% 
    simplify2array %>% 
    apply(1, var, na.rm = TRUE)
  
  output_list[["simpson-variance"]]  <- list_of_fitted_models %>% 
    lapply(function(x) x$simpson) %>% 
    simplify2array %>% 
    apply(1, var, na.rm = TRUE)
  
  output_list[["bray-curtis-variance"]] <- list_of_fitted_models %>% 
    lapply(function(x) x[["bray-curtis"]]) %>% 
    simplify2array %>% 
    apply(1:2, var, na.rm = TRUE)
  
  output_list[["euclidean-variance"]] <- list_of_fitted_models %>% 
    lapply(function(x) x$euclidean) %>% 
    simplify2array %>% 
    apply(1:2, var, na.rm = TRUE)
  
  output_list
}

### nonparametric bootstrap
nonparametric_variance <- function(W, 
                                   X, 
                                   tuning,
                                   perturbation, 
                                   network,
                                   base,
                                   ncores,
                                   nsub, 
                                   ...) {
  curly_b <- sample(1:(nrow(W)), size = nsub, replace=T)
  
  fitted_model <- fit_aitchison(W[curly_b, ],  
                                X[curly_b, ], 
                                tuning = tuning,
                                perturbation = perturbation, 
                                network = network,
                                base = base,
                                ncores = ncores, 
                                ...)
  eY <- matrix(NA, ncol = ncol(W)-1, nrow = nrow(W))
  eY[curly_b, ] <- fitted_model$X %*% fitted_model$beta + 
    matrix(fitted_model$beta0 , ncol = ncol(W)-1, nrow = nsub, byrow=T)
  toCompositionMatrix(Y=eY, base=base) %>% get_diversities
  
}


### parametric bootstrap
parametric_variance <- function(fitted_aitchison, 
                                W, X, 
                                tuning,
                                perturbation, 
                                network,
                                base,
                                ncores, ...) {
  
  # unfortunately in this model we condition on M_i, so no randomising here
  #ms <- apply(W, 1, sum)
  ms <- rowSums(W)
  
  mw <- make_w(mu=fitted_aitchison$fitted_y, 
               Sigma = fitted_aitchison$sigma, 
               mm=ms)
  
  fitted_model <- fit_aitchison(mw,  
                                X,
                                tuning = tuning,
                                perturbation = perturbation, 
                                network = network,
                                base = base,
                                ncores = ncores, ...)
  get_diversities(fitted_model$fitted_z)
}
