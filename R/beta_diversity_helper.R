#' A function to facilitate viewing and manipulating beta diversity estimates
#'
#' @param dv The output of a DivNet() call
#' @param physeq The phyloseq object containing the sample data and  abundance table
#' @param measure The beta diversityindex of interest
#' @param x The covariate
#'
#' @importFrom phyloseq sample_data get_variable
#' @importFrom tibble rownames_to_column add_column
#' @importFrom tidyr gather
#' @importFrom dplyr select distinct mutate
#'
#' @return A data frame with the ecosystems, beta diversity estimates, and CIs
#'
#' @export
simplifyBeta <- function(dv,
                         physeq,
                         measure,
                         x) {

  if(measure == "Aitchison"){
    dv[[measure]] <- get_aitchison_distance(dv$fitted_z)
    rownames(dv[[measure]]) <- rownames(dv[["bray-curtis"]])
    colnames(dv[[measure]]) <- colnames(dv[["bray-curtis"]])
  }

  Covar2 <- Covar1 <- beta_est <- beta_var <- Sample2 <- Sample1 <- NULL

  in_sample_names <- physeq %>% sample_names
  changed_sample_names <- dv[[measure]] %>% data.frame %>% names

  if ( ! all(in_sample_names == changed_sample_names) ) {
    stop("Your sample names contain non-standard characters, please change that :)")
  }

  beta_var_matrix <- dv[[paste(measure, "-variance", sep = "")]]

  vars <- physeq %>% sample_data %>% get_variable(x)
  names(vars) <- physeq %>% sample_names

  if(!is.null(beta_var_matrix)){
    dv[[measure]] %>%
      data.frame %>%
      rownames_to_column("Sample1") %>%
      gather("Sample2", "beta_est", names(vars)) %>%
      add_column("beta_var" = beta_var_matrix %>%
                   data.frame %>%
                   rownames_to_column("Sample1") %>%
                   gather("Sample2", "var", names(vars)) %>%
                   select("var") %>% c %>% unlist) %>%
      mutate("Covar1" = vars[Sample1],
             "Covar2" = vars[Sample2]) %>%
      select(Covar1, Covar2, beta_est, beta_var) %>%
      dplyr::filter(beta_est > 1e-16) %>%
      unique %>%
      distinct(beta_est, .keep_all = TRUE) %>%
      mutate("lower" = pmax(0, beta_est - 2*sqrt(beta_var)),
             "upper" = pmax(0, beta_est + 2*sqrt(beta_var))) %>%
      return()
  } else{
    dv[[measure]] %>%
      data.frame %>%
      rownames_to_column("Sample1") %>%
      gather("Sample2", "beta_est", names(vars)) %>%
      mutate("Covar1" = vars[Sample1],
             "Covar2" = vars[Sample2]) %>%
      select(Covar1, Covar2, beta_est) %>%
      dplyr::filter(beta_est > 1e-16) %>%
      unique %>%
      distinct(beta_est, .keep_all = TRUE) %>%
      return()
  }
}
