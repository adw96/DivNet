library(DivNet)
library(magrittr)
context("Test phylodivnet")

data(GlobalPatterns)
GP = GlobalPatterns
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
GP1 = prune_taxa(wh0, GP)

#object <- DivNet::phylodivnet(GP1, X = "SampleType")
# Unable to test phylodivnet as it appears to be an incompleted function? 

