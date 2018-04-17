set.seed(1)
my_counts <- matrix(rpois(30, lambda=10), nrow = 6)
rownames(my_counts) <- paste("Sample", 1:6, sep = "")
colnames(my_counts) <- paste("Taxon", 1:5, sep = "")

# 6 Samples, 5 taxa
my_counts

my_covariate <- cbind(1, rep(c(0,1), each = 3), rep(c(0,1), 3))
my_covariate



library(DivNet)
divnet(my_counts, NULL)
divnet(my_counts, X =my_covariate[,3]) # help


fit_aitchison(my_counts)
fit_aitchison(my_counts, my_covariate)

divnet(my_counts, X =my_covariate[,3]) # help

divnet(my_counts, X =my_covariate[,2:3])
divnet(my_counts, NULL)
divnet(my_counts, X =my_covariate[,1])
divnet(my_counts, my_covariate[,3])
divnet(my_counts, my_covariate, variance="parametric")
