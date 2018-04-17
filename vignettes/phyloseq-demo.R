library(phyloseq)
data(GlobalPatterns)

library(DivNet)

gp <- GlobalPatterns %>% 
  subset_samples(SampleType == "Mock") %>%
  subset_taxa(Order=="Cenarchaeales") %>% 
  divnet
gp