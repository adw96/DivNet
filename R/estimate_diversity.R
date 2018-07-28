
#' divnet
#' 
#' @param W An abundance table with taxa as columns and samples as rows; or a phyloseq object
#' @param X The covariate matrix, with samples as rows and variables as columns. Defaults to NULL (samples are biological replicates).
#' @param fitted_model object produced by fit_aitchison. Defaults to NULL.
#' @param tuning settings for tuning the MC-MH algorithm. Options include NULL (defaults to "fast"), "fast", "careful" or a named list with components EMiter (number of EM iterations; 6 for fast, 10 for careful), EMburn (number of EM iterations to burn; 3 for fast, 5 for careful), MCiter (number of MC iterations; 500 for fast, 1000 for careful), MCburn (number of MC iterations to burn; 250 for fast, 500 for careful) and stepsize (variance used for MH samples; 0.01 for both fast and careful)
#' @param perturbation Perturbation magnitude for zero values when calculating logratios.
#' @param network How to estimate network. Defaults to NULL (the default), "default" (generalised inverse, aka naive). Other options include "diagonal", "stars" (requires glasso and SpiecEasi to be installed), or a function that you want to use to estimate the network
#' @param base Base taxon index.
#' @param ncores Number of cores to use for parallelization
#' @param variance method to get variance of estimates. Current options are "parametric" for parametric bootstrap, "nonparametric" for nonparametric bootstrap, and "none" for no variance estimates
#' @param B Number of bootstrap iterations for estimating the variance.
#' @param nsub Number of subsamples for nonparametric bootstrap. Defaults to half the number of observed samples.
#' @param ... Additional parameters to be passed to the network function
#' 
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
                    ...) {
  
  if ("phyloseq" %in% class(W)) {
    
    input_data <- W
    
    W <- input_data %>% otu_table %>% as.matrix
    suppressWarnings({class(W) <- "matrix"})
    
    if (phyloseq::taxa_are_rows(input_data)) W <- W %>% t
    
    samples_names <- input_data %>% sample_names
    
    # make the design matrix
    if (is.character(X)) {
      X <- make_design_matrix(input_data, X)
    } else if (is.null(X)) {
      xx <- input_data %>% sample_data %>% rownames 
      X <- model.matrix(~xx)
    }
  } else {
    samples_names <- rownames(W)
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
  }
  
  output_list[["X"]] <- X
  
  # Adding class alpha-estimates
  if (!is.null(output_list$shannon)) {
    if (!is.null(output_list$`shannon-variance`)) {
      output_list$shannon <- mapply(breakaway:::alpha_estimate, 
                                    estimate = output_list$shannon, 
                                    error = output_list$`shannon-variance`,
                                    estimand = "Shannon",
                                    name = "DivNet",
                                    #interval = NULL,
                                    #interval_type = NULL,
                                    #type = NULL,
                                    model = "Aitchison",
                                    #warnings = NULL,
                                    frequentist = TRUE,
                                    parametric = TRUE,
                                    reasonable = TRUE,
                                    other = list(fitted_model = fitted_model),
                                    SIMPLIFY = F) %>%
        breakaway:::alpha_estimates
      output_list$`shannon-variance` <- NULL
    } else {
      output_list$shannon <- mapply(breakaway:::alpha_estimate, 
                                    estimate = output_list$shannon, 
                                    #error = output_list$`shannon-variance`,
                                    estimand = "Shannon",
                                    name = "DivNet",
                                    #interval = NULL,
                                    #interval_type = NULL,
                                    #type = NULL,
                                    model = "Aitchison",
                                    #warnings = NULL,
                                    frequentist = TRUE,
                                    parametric = TRUE,
                                    reasonable = TRUE,
                                    other = list(fitted_model = fitted_model),
                                    SIMPLIFY = F) %>%
        breakaway:::alpha_estimates
    }

  }
  if (!is.null(output_list$simpson)) {
    if (!is.null(output_list$`simpson-variance`)) {
      output_list$simpson <- mapply(breakaway:::alpha_estimate, 
                                    estimate = output_list$simpson, 
                                    error = output_list$`simpson-variance`, 
                                    estimand = "Simpson",
                                    name = "DivNet",
                                    #interval = NULL,
                                    #interval_type = NULL,
                                    #type = NULL,
                                    model = "Aitchison",
                                    #warnings = NULL,
                                    frequentist = TRUE,
                                    parametric = TRUE,
                                    reasonable = TRUE,
                                    other = list(fitted_model = fitted_model),
                                    SIMPLIFY = F) %>%
        breakaway:::alpha_estimates
      output_list$`simpson-variance` <- NULL
    } else {
      output_list$simpson <- mapply(breakaway:::alpha_estimate, 
                                    estimate = output_list$simpson, 
                                    #error = output_list$`simpson-variance`, 
                                    estimand = "Simpson",
                                    name = "DivNet",
                                    #interval = NULL,
                                    #interval_type = NULL,
                                    #type = NULL,
                                    model = "Aitchison",
                                    #warnings = NULL,
                                    frequentist = TRUE,
                                    parametric = TRUE,
                                    reasonable = TRUE,
                                    other = list(fitted_model = fitted_model),
                                    SIMPLIFY = F) %>%
        breakaway:::alpha_estimates
    }
  }
  

  class(output_list) <- c("diversityEstimates", class(output_list))

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
  compositions <- to_composition_matrix(Y=eY, base=base) 
  get_diversities(compositions)
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
               mm=ms, base = base)
  
  fitted_model <- fit_aitchison(mw,  
                                X,
                                tuning = tuning,
                                perturbation = perturbation, 
                                network = network,
                                base = base,
                                ncores = ncores, ...)
  get_diversities(fitted_model$fitted_z)
}
