
#' divnet
#' 
#' @param W An abundance table with taxa as columns and samples as rows; or a phyloseq object. 
#' @param X The covariate matrix, with samples as rows and variables as columns. Defaults to NULL (sample_names are the covariates).
#' @param fitted_model object produced by fit_aitchison. Defaults to NULL.
#' @param tuning settings for tuning the MC-MH algorithm. Options include NULL (defaults to "fast"), "fast", "careful" or a named list with components EMiter (number of EM iterations; 6 for fast, 10 for careful), EMburn (number of EM iterations to burn; 3 for fast, 5 for careful), MCiter (number of MC iterations; 500 for fast, 1000 for careful), MCburn (number of MC iterations to burn; 250 for fast, 500 for careful) and stepsize (variance used for MH samples; 0.01 for both fast and careful)
#' @param perturbation Perturbation magnitude for zero values when calculating logratios.
#' @param network How to estimate network. Defaults to NULL (the default), "default" (generalised inverse, aka naive). Other options include "diagonal", "stars" (requires glasso and SpiecEasi to be installed), or a function that you want to use to estimate the network
#' @param base The column index of the base taxon in the columns of W, or the name of the taxon (must be a column name of W, or a taxon name if W is a phyloseq object). If NULL, will use `pick_base` to choose a taxon. If no taxa are observed in all samples, an error will be thrown. In that case, we recommend trying a number of different highly abundant taxa to confirm the results are robust to the taxon choice. 
#' @param ncores Number of cores to use for parallelization
#' @param variance method to get variance of estimates. Current options are "parametric" for parametric bootstrap, "nonparametric" for nonparametric bootstrap, and "none" for no variance estimates
#' @param B Number of bootstrap iterations for estimating the variance.
#' @param nsub Number of subsamples for nonparametric bootstrap. Defaults to half the number of observed samples.
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted. Optional, defaults to \code{NULL}. Formula objects must match column names found in \code{X}.
#' @param ... Additional parameters to be passed to the network function
#' 
#' @importFrom breakaway make_design_matrix
#' @importFrom magrittr "%>%"
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq sample_names
# #' @importFrom breakaway alpha_estimate
# #' @importFrom breakaway alpha_estimates
# #' @importClassesFrom phyloseq phyloseq
# #' @importClassesFrom breakaway alpha_estimate alpha_estimates
#' 
#' @import phyloseq
#' @import breakaway
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
                    formula = NULL,
                    ...) {
  
  if (!is.null(formula)) {
    # get model.matrix
    X <- data.frame(X)
    X <- stats::model.matrix(object = formula, data = X)
  }
  
  if ("phyloseq" %in% class(W)) {
    
    input_data <- W
    
    W <- input_data %>% otu_table %>% as.matrix
    suppressWarnings({class(W) <- "matrix"})
    
    if (phyloseq::taxa_are_rows(input_data)) W <- W %>% t
    
    samples_names <- input_data %>% sample_names
    
    # make the design matrix
    if (is.character(X)) {
      X <- breakaway::make_design_matrix(input_data, X)
    } else if (is.null(X)) {
      xx <- input_data %>% sample_data %>% rownames 
      X <- model.matrix(~xx)
    }
  } else if ("otu_table" %in% class(W)) {
    
    input_data <- W
    W <- input_data %>% as.matrix
    suppressWarnings({class(W) <- "matrix"})
    
    if (phyloseq::taxa_are_rows(input_data)) W <- W %>% t
    
    samples_names <- input_data %>% sample_names
    
  } else {
    samples_names <- rownames(W)
  } 
  

  if (is.null(X)) {
    X <- matrix(1, ncol=1, nrow=nrow(W))
  }
  
  # remove taxa that weren't observed 
  # yes, this is a good idea
  if (any(colSums(W) == 0)) {
    message("Removing absent taxa!")
    W <- W[ , which(colSums(W) > 0)]
  }
  
  if (nrow(W) == 1) {
    stop("DivNet requires more than 1 sample")
  }
  if (ncol(W) == 2) {
    stop("Cannot fit a network model with 2 taxa")
  }
  
  if (class(base) == "character") {
    if (base %in% colnames(W)) {
      base <- which(base == colnames(W))
    } else {
      stop("base not found in taxon names. base taxon may be unobserved in all samples, or you may have a typo.")
    }
  }
  
  
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
  
  base <- fitted_model$base
  
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
    
    for(i in 1:length(variance_estimates)) {
      output_list[[names(variance_estimates)[i]]] <-  variance_estimates[[i]]
    }
    
    # Add variance to alpha_diversity class
    for(i in 1:nrow(W)) {
      output_list$shannon[[i]]$error <- sqrt(variance_estimates$`shannon-variance`)[i]
      output_list$shannon[[i]]$interval <- c(output_list$shannon[[i]]$estimate - 2*output_list$shannon[[i]]$error,
                                             output_list$shannon[[i]]$estimate + 2*output_list$shannon[[i]]$error)
      output_list$shannon[[i]]$interval_type <- "symmetric"
      
      output_list$simpson[[i]]$error <- sqrt(variance_estimates$`simpson-variance`)[i]
      output_list$simpson[[i]]$interval <- c(output_list$simpson[[i]]$estimate - 2*output_list$simpson[[i]]$error,
                                             output_list$simpson[[i]]$estimate + 2*output_list$simpson[[i]]$error)
      output_list$simpson[[i]]$interval_type <- "symmetric"
    }
    variance_estimates$`shannon-variance` <- variance_estimates$`simpson-variance` <- NULL
    
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
    
    for(i in 1:length(variance_estimates)) {
      output_list[[names(variance_estimates)[i]]] <-  variance_estimates[[i]]
    }
  } else if (variance == 0 | variance == "none") {
    # do nothing, no warning
  } else {
    warning("No variance estimate is computed")
  }
  
  output_list[["X"]] <- X
  output_list[["fitted_z"]] <- zz
  
  class(output_list) <- c("diversityEstimates", class(output_list))
  
  output_list
}