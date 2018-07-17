# data(Lee)
# dn <- phyloseq::tax_glom(Lee, taxrank="Phylum") %>% 
#   divnet(tuning = "test")
# dn %>% names
make_alpha_estimates <- function(dn) {
  my_alpha <- list()
  
  my_alpha$shannon <- mapply(breakaway::alpha_estimate, 
                             estimate = dn$shannon, 
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
    alpha_estimates
  
  my_alpha$simpson <- mapply(breakaway::alpha_estimate, 
                             estimate = dn$simpson, 
                             error = dn$`simpson-variance`,
                             estimand = "Simpson",
                             name = "DivNet",
                             model = "Aitchison",
                             frequentist = TRUE,
                             parametric = TRUE,
                             reasonable = TRUE,
                             SIMPLIFY = F) %>%
    alpha_estimates
  
  my_alpha
}
# make_alpha_estimates(dn)$shannon 
