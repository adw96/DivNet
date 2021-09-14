

#' Print function
#'
#' @param x An object of class diversityEstimates
#' @param ... other arguments to be passed to print
#' @return NULL
#'
#' @export
print.diversityEstimates <- function(x, ...) {
  dv <- x
  cat("An object of class `diversityEstimates` with the following elements:\n")
  sapply(1:length(names(dv)), function(i) { cat("  - ", names(dv)[i], "\n")})
  cat("Access individual components with, e.g., object$shannon and object$`shannon-variance`\n")
  cat("Use function testDiversity() to test hypotheses about diversity")
}



#' Plot function
#'
#' TODO make more like the phyloseq plot richness
#'
#' @param x An object of class diversityEstimates
#' @param ... other arguments to be passed to plot
#' @return An object of class ggplot
#' @export
plot.diversityEstimates <- function(x, ...) {
  dv <- x
  args <- match.call(expand.dots = TRUE)
  if (is.null(args$xx)) {
    args$xx <- "samples"
  }
  if (is.null(args$h0)) {
    args$h0 <- "shannon"
  }
  xx <- args$xx
  h0 <- args$h0
  if (h0 %in% c("shannon", "simpson")) {
    ests <- sapply(dv[[h0]], function(x) x$estimate)
    # vars <- sapply(dv[[h0]], function(x) x$error)
    lci <- sapply(dv[[h0]], function(x) x$interval[1])
    uci <- sapply(dv[[h0]], function(x) x$interval[2])
    df <- data.frame("names" = names(ests),
                     "h0" = ests, lci, uci, dv$X)
  } else {
    lci <- dv[[h0]] - 2*sqrt(dv[[paste(h0, "-variance", sep = "")]])
    uci <- dv[[h0]] + 2*sqrt(dv[[paste(h0, "-variance", sep = "")]])
    df <- data.frame("names" = names(dv[[h0]]),
                     "h0" = dv[[h0]], lci, uci, dv$X)
  }

  df$names <- factor(df$names, levels = df$names)

  ggplot2::ggplot(df, ggplot2::aes(x = names, xend = names)) +
    ggplot2::geom_point(ggplot2::aes(x = names, y = h0)) +
    ggplot2::geom_segment(ggplot2::aes(y = lci, yend = uci)) +
    ggplot2::ylab(paste(h0, "estimate")) +
    ggplot2::xlab(xx) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Test diversity
#'
#' Hypothesis testing for alpha-diversity.
#'
#'
#' @references Willis, A., Bunge, J., and Whitman, T. (2017). Improved detection of changes in species richness in high diversity microbial communities. \emph{JRSS-C.}
#'
#' @param dv An object of class diversityEstimates. The variable `X` used for the construction
#' @param h0 The alpha-diversity index to be tested for equality
#' @return A data frame similar to the output of `lm`
#'
#' @export
testDiversity <- function(dv, h0 = "shannon") {
  cat("Hypothesis testing:\n")
  if (h0 %in% c("shannon", "simpson")) {
    bt <- breakaway::betta(sapply(dv[[h0]], function(x) x$estimate),
                           sapply(dv[[h0]], function(x) x$error),
                           X = dv[["X"]])
  } else {
    bt <- breakaway::betta(dv[[h0]], dv[[paste(h0, "-variance", sep="")]], X = dv[["X"]])
  }
  cat(paste("  p-value for global test:", bt$global[2], "\n"))
  bt$table
}

#' Test beta diversity
#'
#' Hypothesis testing for beta-diversity.
#'
#' This function uses output from DivNet() to estimate community centroids
#' within groups defined by the groups argument and test a null hypothesis
#' of equality of all group centroids against a general alternative. This test
#' is conducted using a pseudo-F statistic with null distribution approximated
#' via a nonparametric bootstrap.
#'
#' For more details and suggested workflow see the beta diversity vignette:
#' \code{vignette("beta_diversity", package = "DivNet")}
#'
#' @param dv An object of class diversityEstimates. The variable `X` used for the construction
#' @param h0 The beta-diversity index to be tested for equality
#' @param groups A numeric vector giving group membership of each specimen
#' @param sample_specimen_matrix A matrix with ik-th entry 1 if the i-th sequenced sample is taken from specimen k, 0 otherwise.
#' The columns of this matrix should correspond to unique specimens and must be named.
#' @param n_boot Number of (cluster) bootstrap resamples to use
#' @return A list containing the observed pseudo-F statistic, the beta diversity used, the
#' p-value returned by the bootstrapped pseudo-F test of equality of (measured) centroids,
#' a vector of computed bootstrapped test statistics, a matrix of estimated group centroids,
#' and a list of group centroids estimated from each bootstrap resample
#' #'
#' @export
testBetaDiversity <- function(dv,
                              h0,
                              groups,
                              sample_specimen_matrix,
                              n_boot = 1000){

  if(length(colnames(sample_specimen_matrix)) != ncol(sample_specimen_matrix)){
            stop("Columns of argument sample_specimen_matrix must be named.
            Recommended column names are names of unique specimens in your data.")
}


  n_groups <- length(unique(groups))
  unique_groups <- unique(groups)
  unique_specimens <- colnames(sample_specimen_matrix)
  n_specimens <- ncol(sample_specimen_matrix)
  group_specimens <- sapply(unique_groups,
                            function(x) apply(sample_specimen_matrix[groups == x,,drop = F],2,max) %>%
                              (function(y) names(y)[y==1]))


if(h0 == "bray-curtis"){
  bc_matrix <- dv$`bray-curtis`

  observed_test_statistic <- get_bc_test_statistic(bc_mat = bc_matrix, groups, unique_groups,
                                                               n_groups,
                                                               n_specimens)


  boot_test_statistics <- numeric(n_boot)

  np_boot_pulls <-replicate(n_boot,
                            sample(1:ncol(sample_specimen_matrix),
                                   ncol(sample_specimen_matrix), replace = T))


  group_centroids <- lapply(unique_groups,
                            function(gr){
                              samples <- sapply(group_specimens[[gr]],
                                                function(specname)
                                                  which(sample_specimen_matrix[,specname] ==1)[1])
                              return(apply(dv$fitted_z[samples,],2,median))}
  )

  names(group_centroids) <- unique_groups


  boot_test_statistics <- numeric(n_boot)

  centroid_matrix <- do.call(rbind, lapply(groups,
                                           function(k) group_centroids[[k]]))


  boot_centroids <- vector(n_boot, mode = "list")
  for(k in 1:n_boot){

    which_samples <- do.call(c,lapply(np_boot_pulls[,k],
                                        function(x) which(sample_specimen_matrix[,x] ==1)))

    comps <- dv$fitted_z[which_samples,]

    boot_group_specimens <-sapply(unique_groups,
           function(x) apply(sample_specimen_matrix[groups == x,np_boot_pulls[,k]],2,max) %>%
             (function(y) names(y)[y==1]))

    boot_centroids[[k]] <- lapply(unique_groups,
                                  function(gr){
                                    samples <- unlist(sapply(boot_group_specimens[[gr]],
                                                      function(specname)
                                                        which(sample_specimen_matrix[,specname] ==1)))
                                    return(apply(dv$fitted_z[samples,,drop= F],2,median))}
    )

    names(boot_centroids[[k]]) <- unique_groups



    centered_comps <- comps - centroid_matrix[which_samples,]

    boot_mat <- matrix(0,
                     ncol = nrow(centered_comps),
                     nrow = nrow(centered_comps))

    for(i in 1:(nrow(centered_comps) - 1)){
      for(j in (i + 1):nrow(centered_comps)){
        boot_mat[i,j] <- boot_mat[j,i] <- 0.5*sum(abs(centered_comps[i,] - centered_comps[j,]))
      }
    }

    boot_test_statistics[k] <- get_bc_test_statistic(bc_mat = boot_mat,groups = groups[which_samples],
                                                     unique_groups = unique_groups, n_groups = n_groups,
                                                     n_specimens = n_specimens)


  }

}

if(h0 == "euclidean"){
  euc_matrix <- dv$'euclidean'

  observed_test_statistic <- get_euc_test_statistic(euc_mat = euc_matrix, groups, unique_groups,
                                                   n_groups,
                                                   n_specimens)

  boot_test_statistics <- numeric(n_boot)

  np_boot_pulls <-replicate(n_boot,
                            sample(1:ncol(sample_specimen_matrix),
                                   ncol(sample_specimen_matrix), replace = T))


  group_centroids <- lapply(unique_groups,
                            function(gr){
                              samples <- sapply(group_specimens[[gr]],
                                                function(specname)
                                                  which(sample_specimen_matrix[,specname] ==1)[1])
                              return(apply(dv$fitted_z[samples,,drop= F],2,mean))}
  )

  names(group_centroids) <- unique_groups


  boot_test_statistics <- numeric(n_boot)

  centroid_matrix <- do.call(rbind, lapply(groups,
                                           function(k) group_centroids[[k]]))


  boot_centroids <- vector(n_boot, mode = "list")
  for(k in 1:n_boot){


    which_samples <- do.call(c,lapply(np_boot_pulls[,k],
                                      function(x) which(sample_specimen_matrix[,x] ==1)))

    comps <- dv$fitted_z[which_samples,]

    boot_group_specimens <-sapply(unique_groups,
                                  function(x) apply(sample_specimen_matrix[groups == x,np_boot_pulls[,k]],2,max) %>%
                                    (function(y) names(y)[y==1]))

    boot_centroids[[k]] <- lapply(unique_groups,
                                  function(gr){
                                    samples <- unlist(sapply(boot_group_specimens[[gr]],
                                                             function(specname)
                                                               which(sample_specimen_matrix[,specname] ==1)))
                                    return(apply(dv$fitted_z[samples,,drop = F],2,mean))})

    names(boot_centroids[[k]]) <- unique_groups


    centered_comps <- comps - centroid_matrix[which_samples,]

    boot_mat <- matrix(0,
                       ncol = nrow(centered_comps),
                       nrow = nrow(centered_comps))

    for(i in 1:(nrow(centered_comps) - 1)){
      for(j in (i + 1):nrow(centered_comps)){
        boot_mat[i,j] <- boot_mat[j,i] <- sqrt(sum((centered_comps[i,] - centered_comps[j,])^2))
      }
    }

    boot_test_statistics[k] <- get_euc_test_statistic(euc_mat = boot_mat,groups = groups[which_samples],
                                                     unique_groups = unique_groups, n_groups = n_groups,
                                                     n_specimens = n_specimens)


  }
}

if(h0 == "aitchison"){

  aitch_matrix <- get_aitchison_distance(dv$fitted_z)

  observed_test_statistic <- get_euc_test_statistic(aitch_matrix, groups, unique_groups,
                                                    n_groups,
                                                    n_specimens)

  group_centroids <- lapply(unique_groups,
                            function(gr){
                            samples <- sapply(group_specimens[[gr]],
                                              function(specname)
                                              which(sample_specimen_matrix[,specname] ==1)[1])
                            return(apply(log_ratio(dv$fitted_z[samples,]),2,mean))}
                              )

  names(group_centroids) <- unique_groups


  boot_test_statistics <- numeric(n_boot)

  np_boot_pulls <-replicate(n_boot,
                            sample(1:ncol(sample_specimen_matrix),
                                   ncol(sample_specimen_matrix), replace = T))

  centroid_matrix <- do.call(rbind, lapply(groups,
                                                function(k) group_centroids[[k]]))

  boot_groups <- groups

  boot_centroids <- vector(n_boot,mode = "list")

  for(k in 1:n_boot){

    which_samples <- do.call(c,lapply(np_boot_pulls[,k],
                              function(x) which(sample_specimen_matrix[,x] ==1)))

    comps <- log_ratio(dv$fitted_z[which_samples,])

    boot_group_specimens <-sapply(unique_groups,
                                  function(x) apply(sample_specimen_matrix[groups == x,np_boot_pulls[,k]],2,max) %>%
                                    (function(y) names(y)[y==1]))

    boot_centroids[[k]] <- lapply(unique_groups,
                                  function(gr){
                                    samples <- unlist(sapply(boot_group_specimens[[gr]],
                                                             function(specname)
                                                               which(sample_specimen_matrix[,specname] ==1)))
                                    return(apply(log_ratio(dv$fitted_z[samples,,drop = F]),2,mean))})
    names(boot_centroids[[k]]) <- unique_groups


    centered_comps <- comps - centroid_matrix[which_samples,]

    boot_mat <- matrix(0,
                       ncol = nrow(centered_comps),
                       nrow = nrow(centered_comps))

    for(i in 1:(nrow(centered_comps) - 1)){
      for(j in (i + 1):nrow(centered_comps)){
        boot_mat[i,j] <- boot_mat[j,i] <- sqrt(sum((centered_comps[i,] - centered_comps[j,])^2))
      }
    }

    boot_test_statistics[k] <- get_euc_test_statistic(euc_mat = boot_mat,groups = groups[which_samples],
                                                      unique_groups = unique_groups , n_groups = n_groups,
                                                      n_specimens = n_specimens)


  }

}

  p.val <- mean(boot_test_statistics >= observed_test_statistic)
  if(p.val == 0){
    p.val <- paste(" < ", signif(1/n_boot,2),sep = "", collapse = "")
  }

centroids <- do.call(rbind,group_centroids)
rownames(centroids) <- unique_groups


return(list("Test statistic" = observed_test_statistic,
              "h0" = h0,
              "p_value" = p.val,
              "bootstrapped_statistics" = boot_test_statistics,
              "centroids" = centroids,
              "boot_centroids" = boot_centroids
            ))


}



get_bc_test_statistic <- function(bc_mat, groups, unique_groups,
                      n_groups,
                      n_specimens){

test_statistic_numerator <- 0
test_statistic_denominator <- 0

for(group in unique_groups){
  sub_matrix <- bc_mat[groups == group,groups == group]
  test_statistic_denominator <- test_statistic_denominator + sum(sub_matrix[upper.tri(sub_matrix)])
  test_statistic_numerator <- test_statistic_numerator + sum(bc_mat[groups == group,groups != group])
}

observed_test_statistic <- (test_statistic_numerator/(n_groups - 1))/(test_statistic_denominator/(n_specimens - n_groups - 1))
return(observed_test_statistic)
}



get_euc_test_statistic <- function(euc_mat, groups, unique_groups,
                                  n_groups,
                                  n_specimens){

  euc_mat <- euc_mat^2 #squared distances for Euclidean distance test

  test_statistic_numerator <- 0
  test_statistic_denominator <- 0

  for(group in unique_groups){
    sub_matrix <- euc_mat[groups == group,groups == group]
    test_statistic_denominator <- test_statistic_denominator + sum(sub_matrix[upper.tri(sub_matrix)])
    test_statistic_numerator <- test_statistic_numerator + sum(euc_mat[groups == group,groups != group])
  }

  observed_test_statistic <- (test_statistic_numerator/(n_groups - 1))/(test_statistic_denominator/(n_specimens - n_groups - 1))
  return(observed_test_statistic)
}


get_aitchison_distance <- function(comp_matrix){
lr_matrix <- log_ratio(comp_matrix)
return(as.matrix(dist(lr_matrix)))

}

log_ratio <- function(comp_matrix){
  lr_matrix <- log(comp_matrix)
  lr_matrix <- lr_matrix -matrix(apply(lr_matrix,1, mean),ncol = 1)%*%matrix(1, ncol = ncol(lr_matrix))
  return(lr_matrix)
}
