get_diversities <- function(zz, samples_names = NULL) {

  ## Note: at this point, variances have not been
  # estimated, and so the error and interval fields
  # are empty for now

  output_list <- list()

  # Estimate Shannon diversity
  shannon_tmp <- apply(zz, 1, shannon_true)
  output_list[["shannon"]] <- breakaway::alpha_estimates(mapply(breakaway::alpha_estimate,
                                                                estimate = shannon_tmp,
                                                                #error = shannon_sd,
                                                                estimand = "Shannon",
                                                                name = "DivNet",
                                                                # interval = c(shannon_tmp - 2*shannon_sd,
                                                                #             shannon_tmp + 2*shannon_sd),
                                                                # interval_type = "symmetric",
                                                                #type = NULL,
                                                                model = "Aitchison",
                                                                #warnings = NULL,
                                                                frequentist = TRUE,
                                                                parametric = TRUE,
                                                                reasonable = TRUE,
                                                                # other = list(fitted_model = fitted_model),
                                                                SIMPLIFY = F))

  # Estimate Simpson
  simpson_tmp <- apply(zz, 1, simpson_true)
  output_list[["simpson"]] <- breakaway::alpha_estimates(mapply(breakaway::alpha_estimate,
                                                                estimate = simpson_tmp,
                                                                #error = simpson_sd,
                                                                estimand = "Simpson",
                                                                name = "DivNet",
                                                                #interval = c(simpson_tmp - 2*simpson_sd,
                                                                #             simpson_tmp + 2*simpson_sd),
                                                                #interval_type = "symmetric",
                                                                #type = NULL,
                                                                model = "Aitchison",
                                                                #warnings = NULL,
                                                                frequentist = TRUE,
                                                                parametric = TRUE,
                                                                reasonable = TRUE,
                                                                # other = list(fitted_model = fitted_model),
                                                                SIMPLIFY = F))

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
    lapply(function(x) sapply(x$shannon, function(x) x$estimate)) %>%
    simplify2array %>%
    apply(1, var, na.rm = TRUE)

  output_list[["simpson-variance"]]  <- list_of_fitted_models %>%
    lapply(function(x) sapply(x$simpson, function(x) x$estimate)) %>%
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
                                   return_compositions = FALSE,
                                   cluster_on_X = FALSE,
                                   ...) {
  if(cluster_on_X){
    if(!all(unique(as.numeric(X)) %in% c(0,1))){
      stop("All values in X must be either 1 or 0 if cluster_on_X is TRUE.")
    }
    curly_b_pre <- sample(1:ncol(X),
                          size = ncol(X),
                          replace = TRUE)

    curly_b <- lapply(curly_b_pre,
                      function(x) which(X[,x]==1))

    curly_b <- do.call(c,curly_b)


  } else{
    curly_b <- sample(1:(nrow(W)), size = nsub, replace=T)
  }

  fitted_model <- fit_aitchison(W[curly_b, ],
                                X[curly_b, ],
                                tuning = tuning,
                                perturbation = perturbation,
                                network = network,
                                base = base,
                                ncores = ncores,
                                ...)


  eY <- matrix(NA, ncol = ncol(W)-1, nrow =length(curly_b))

  if(!cluster_on_X){
  eY[curly_b, ] <- fitted_model$fitted_y
  } else{
    eY[1:length(curly_b),] <- fitted_model$fitted_y
  }

  compositions <- to_composition_matrix(Y=eY, base=base)
  if(cluster_on_X){
    rownames(compositions) <- rownames(W)[curly_b] }
  if(!return_compositions){
    get_diversities(compositions)
  } else{
    return(list("diversities" = get_diversities(compositions),
                "compositions" = compositions,
                "samples" = curly_b))
  }

}


### parametric bootstrap
parametric_variance <- function(fitted_aitchison,
                                W, X,
                                tuning,
                                perturbation,
                                network,
                                base,
                                ncores,
                                return_compositions = FALSE,
                                ...) {

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

  if(!return_compositions){
    get_diversities(fitted_model$fitted_z)
  } else{
    return(list("diversities" =   get_diversities(fitted_model$fitted_z),
                "compositions" = fitted_model$fitted_z))
  }
}
