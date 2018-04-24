
#' divnet
#' 
#' @param W TODO
#' @param X TODO
#' @param fitted_model object produced by fit_aitchison
#' @param tuning TODO
#' @param perturbation Perturbation magnitude for zero values when calculating logratios.
#' @param network TODO
#' @param base Base taxon index.
#' @param ncores Number of cores to use for parallelization
#' @param variance method to get variance of estimates. Current options are "parametric" for parametric bootstrap, "nonparametric" for nonparametric bootstrap, and "none" for no variance estimates
#' @param B TODO
#' @param nsub TODO
#' @param ... TODO
#' 
#' @importFrom magrittr "%>%"
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
  
  if ("phyloseq" %in% class(W)) {
    
    input_data <- W
    
    W <- input_data %>% phyloseq::otu_table %>% as.matrix
    suppressWarnings({class(W) <- "matrix"})
    
    if (phyloseq::taxa_are_rows(input_data)) W <- W %>% t
    
    samples_names <- input_data %>% phyloseq::sample_names
    
    # make the design matrix
    if (is.character(X)) {
      X <- make_design_matrix(input_data, X)
    }
  } else {
    samples_names <- rownames(W)
  }
  
  # remove taxa that weren't observed 
  # yes, this is a good idea
  W <- W[ , which(colSums(W) > 0)]
  
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
  output_list <- get_diversities(zz, samples_names)
  
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
    
    variance_estimates <- get_diversity_variance(parametric_list, samples_names)
    
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
    variance_estimates <- get_diversity_variance(nonparametric_list, samples_names)
    output_list <- c(output_list, variance_estimates)
  }
  output_list
}

get_diversities <- function(zz, samples_names = NULL) {
  
  output_list <- list()
  
  # Estimate Shannon diversity
  output_list[["shannon"]] <- apply(zz, 1, shannon_true)
  
  # Estimate Simpson
  output_list[["simpson"]] <- apply(zz, 1, simpson_true)
  
  # Estimate Bray Curtis
  output_list[["bray-curtis"]] <- bray_curtis_true(zz)
  
  # Estimate Euclidean
  output_list[["euclidean"]] <- euclidean_true(zz)
  
  if (!is.null(samples_names)) {
    names(output_list[["shannon"]]) <- samples_names
    names(output_list[["simpson"]])  <- samples_names
    rownames(output_list[["bray-curtis"]]) <- samples_names
    rownames(output_list[["euclidean"]]) <- samples_names
    colnames(output_list[["bray-curtis"]]) <- samples_names
    colnames(output_list[["euclidean"]]) <- samples_names
    
  }
  output_list
}

get_diversity_variance <- function(list_of_fitted_models, samples_names = NULL) {
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
  
  if (!is.null(samples_names)) {
    names(output_list[["shannon-variance"]]) <- samples_names
    names(output_list[["simpson-variance"]])  <- samples_names
    rownames(output_list[["bray-curtis-variance"]]) <- samples_names
    rownames(output_list[["euclidean-variance"]]) <- samples_names
    colnames(output_list[["bray-curtis-variance"]]) <- samples_names
    colnames(output_list[["euclidean-variance"]]) <- samples_names
  }
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
  eY[curly_b, ] <- fitted_model$fitted_y
    # fitted_model$X %*% fitted_model$beta + 
    # matrix(fitted_model$beta0 , ncol = ncol(W)-1, nrow = nsub, byrow=T)
  
  to_composition_matrix(Y=eY, base=base) %>% get_diversities
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
