#' Make alpha estimates
#' 
#' This function is a wrapper for \code{alpha_estimates} in \code{breakaway}.
#' 
#' @param dn A DivNet object. 
#' 
#' @return An object of class \code{alpha_estimate} containing alpha estimate
#' information in \code{dn}. 
#' 
#' export
make_alpha_estimates <- function(dn) {
  my_alpha <- list()
  
  unlisted_shannon <- unlist(dn$shannon)
  shannon_ests <- as.numeric(unlisted_shannon[grep("estimate", 
                                                   names(unlisted_shannon))])
  names(shannon_ests) <- names(dn$`shannon-variance`)
  
  my_alpha$shannon <- mapply(breakaway::alpha_estimate, 
                             estimate = shannon_ests, 
                             #estimate = dn$shannon, 
                             error = dn$`shannon-variance`,
                             estimand = "Shannon",
                             name = "DivNet",
                             model = "Aitchison",
                             # interval = c(dn$shannon - 1.96*sqrt(dn$`shannon-variance`),
                             #              dn$shannon + 1.96*sqrt(dn$`shannon-variance`)),
                             frequentist = TRUE,
                             parametric = TRUE,
                             reasonable = TRUE,
                             SIMPLIFY = F) %>%
    breakaway::alpha_estimates()
  
  unlisted_simpson <- unlist(dn$simpson)
  simpson_ests <- as.numeric(unlisted_simpson[grep("estimate", 
                                                   names(unlisted_simpson))])
  names(simpson_ests) <- names(dn$`simpson-variance`)
  
  my_alpha$simpson <- mapply(breakaway::alpha_estimate, 
                             estimate = simpson_ests, 
                             #estimate = dn$simpson, 
                             error = dn$`simpson-variance`,
                             estimand = "Simpson",
                             name = "DivNet",
                             model = "Aitchison",
                             frequentist = TRUE,
                             parametric = TRUE,
                             reasonable = TRUE,
                             SIMPLIFY = F) %>%
    breakaway::alpha_estimates()
  
  my_alpha
}

