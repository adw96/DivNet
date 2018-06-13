
#' phylodivnet
#' 
#' @param W An abundance table with taxa as columns and samples as rows; or a phyloseq object
#' @param X The covariate matrix, with samples as rows and variables as columns. Note that the estimated UniFrac distance of biological replicates will be zero.
#' @param trees A vector of file paths to files containing the tree estimates in Newick format.
#' @param tuning settings for tuning the MC-MH algorithm. Options include NULL (defaults to "fast"), "fast", "careful" or a named list with components EMiter (number of EM iterations; 6 for fast, 10 for careful), EMburn (number of EM iterations to burn; 3 for fast, 5 for careful), MCiter (number of MC iterations; 500 for fast, 1000 for careful), MCburn (number of MC iterations to burn; 250 for fast, 500 for careful) and stepsize (variance used for MH samples; 0.01 for both fast and careful)
#' @param perturbation Perturbation magnitude for zero values when calculating logratios.
#' @param network How to estimate network. Defaults to NULL (the default), "default" (generalised inverse, aka naive). Other options include "diagonal", "stars" (requires glasso and SpiecEasi to be installed), or a function that you want to use to estimate the network
#' @param base Base taxon index.
#' @param ncores Number of cores to use for parallelization
#' @param B Number of bootstrap iterations for estimating the variance.
#' @param alpha For a 100*(1-alpha) percent confidence interval. Defaults to 0.05.
#' @param ... Additional parameters to be passed to the network function
#' 
#' @importFrom magrittr "%>%"
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq sample_names
#' @importFrom phyloseq UniFrac
#' @importFrom phyloseq read_tree
#' @importFrom phyloseq phyloseq
#' @importFrom phyloseq taxa_names
#' 
#' @author Amy Willis
#' 
#' @export
phylodivnet <-  function(W, 
                         X = NULL, 
                         trees = NULL,
                         tuning = NULL,
                         perturbation = NULL, 
                         network = NULL,
                         base = NULL,
                         ncores = NULL,
                         B = 5,
                         alpha = 0.05,
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
  
  fitted_model <- fit_aitchison(W, 
                                X = X, 
                                tuning = tuning,
                                perturbation = perturbation, 
                                network = network,
                                base = base,
                                ncores = ncores,
                                ...)
  base <- fitted_model$base
  
  ### Find unique rows of the design matrix
  ### to reduce computation
  X_unique <- X[!duplicated(X), ]
  fitted_z <- fitted_model$fitted_z[!duplicated(X), ]
  fitted_y <- fitted_model$fitted_y[!duplicated(X), ]
  fitted_otu <- round(1e10*fitted_z)
  
  if ("phyloseq" %in% class(input_data)) {
    colnames(fitted_otu) <- taxa_names(input_data)
  } else {
    colnames(fitted_otu) <- colnames(input_data)
  }
  ##################################################
  ##### Step 1: Calculate theta-hat(T_i) for all i
  ##################################################
  message("Calculating UniFrac over different trees...\n")
  
  otu_table_fitted <- otu_table(fitted_otu, taxa_are_rows = F)
  
  get_unifrac_i <- function(i) {
    new_tree <- read_tree(trees[i]) 
    ps <- phyloseq(otu_table_fitted, 
                   new_tree)
    UniFrac(ps, weighted = T) %>% as.matrix
  }
  
  if (!is.null(ncores)) {
    message("parallel not supported yet, sorry!")
    # cl <- makeCluster(ncores)  
    # clusterEvalQ(cl, {
    #   library(magrittr)
    #   library(phyloseq)
    #   library(DivNet)
    # })
    # clusterExport(cl, "trees")
    # clusterExport(cl, "otu_table_fitted")
    # # clusterExport(cl, "get_unifrac_i")
    # unifracs <- parSapply(cl, 1:length(trees), 
    #                       get_unifrac_i, 
    #                       otus = otu_table_fitted, 
    #                       simplify = F)
    # stopCluster(cl)
    # 
  } #else {
  unifracs <- sapply(1:length(trees), get_unifrac_i, simplify = F)
  # }
  
  
  #########################
  ##### Step 2: Calculate theta-hat
  #########################
  message("Calculating theta-hat...\n")
  theta_hat <- apply(simplify2array(unifracs), 1:2, mean)
  
  #########################
  ##### Step 3: Parametric boostrap
  #########################
  message("Param BS...\n")
  
  ms <- apply(W, 1, sum)
  parametric_bs <- function(b) {
    message(paste("BS", b, "...\n"))
    ww <- make_w(mu = fitted_y,
                 Sigma = fitted_model$sigma, 
                 mm = ms/2, 
                 base = base)
    fa_temp <- fit_aitchison(ww, 
                             X_unique, 
                             base = base)
    
    fitted_otu2 <- round((1e10*fa_temp$fitted_z))
    
    if ("phyloseq" %in% class(input_data)) {
      colnames(fitted_otu2) <- taxa_names(input_data)
    } else {
      colnames(fitted_otu2) <- colnames(input_data)
    }
    new_tree <- read_tree(sample(trees, 1)) 
    ps <- phyloseq(otu_table(fitted_otu2, 
                             taxa_are_rows = F), 
                   new_tree)
    UniFrac(ps, weighted = T) %>% as.matrix
  }
  # if (!is.null(ncores)) {
  #   message("parallel not supported yet, sorry!")
  #   cl <- makeCluster(ncores)  
  #   clusterEvalQ(cl, {
  #     library(magrittr)
  #     library(phyloseq)
  #     library(DivNet)
  #   })
  #   clusterExport(cl, c("fitted_y", "fitted_model", 
  #                       "ms", "X_unique",
  #                       "base", "trees"))
  #   bs_unifracs <- parSapply(cl, 1:length(trees), parametric_bs, simplify = F)
  #   stopCluster(cl)
  # } else {
    bs_unifracs <- sapply(1:B, parametric_bs, simplify = F)
  # }
  
  message("Finishing up...\n")
  bs_unifracs <- simplify2array(bs_unifracs)
  
  deltas <- sapply(1:B, function(i) {bs_unifracs[, , i] - theta_hat}, simplify = F)
  
  delta_lower <- apply(simplify2array(deltas), 1:2, quantile, probs=c(alpha/2))
  delta_upper <- apply(simplify2array(deltas), 1:2, quantile, probs=c(1-alpha/2))
  
  upper <- theta_hat - delta_lower
  lower <- theta_hat - delta_upper
  
  output_list <- list()
  
  output_list[["lower"]] <- lower
  output_list[["upper"]] <- upper
  output_list
}

