------------------------------------------------------------------------

DivNet <img src="docs/divnet-logo.png" align="right" width="165px"/>
====================================================================

[![Build Status](https://travis-ci.org/adw96/DivNet.svg?branch=master)](https://travis-ci.org/adw96/DivNet) <!-- [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/adw96/DivNet?branch=master&svg=true)](https://ci.appveyor.com/project/adw96/DivNet) --> [![codecov.io](https://codecov.io/gh/adw96/DivNet/coverage.svg?branch=master)](https://codecov.io/gh/adw96/DivNet?branch=master)

DivNet: an R package to estimate diversity when taxa in the community cooccur via a ecological network.

Willis, A.D. and Martin, B.D. (2018+) [*DivNet: Estimating diversity in networked communities*](https://www.biorxiv.org/content/early/2018/04/21/305045). Under review.

About
-----

If you knew the exact composition of a community (such as all of the microbes living on your skin), you could calculate the Shannon diversity of the community, or the Bray-Curtis distance between your skin's microbes and your cat's fur microbes. However, you will only ever observe a small fraction of the microbes on your skin from your experiment, so you have to estimate these diversity indices based on the data, and ideally you will get a confidence interval for those diversity indices.

Almost all ecologists use "plug-in" diversity estimates, obtained by taking the observed relative abundances of each taxon and plugging them into the formula for the true diversity. This is problematic for a number of reasons, including that there may be unobserved taxa, and the observed relative abundances probably don't exactly equal the true relative abundances due to random sampling. Furthermore, typical approaches to obtaining a confidence interval assume that all taxa behave independently: if taxon A is there, this tells you nothing about whether taxon B is likely to be there too. DivNet takes care of all of these issues, and gives you a confidence interval that incorporates all of these challenges.

See [the vignette](https://github.com/adw96/DivNet/blob/master/vignettes/getting-started.Rmd) for a full tutorial, and [the paper](https://www.biorxiv.org/content/early/2018/04/21/305045) for more details.

Want a confidence interval for a different diversity index? Let us know by posting an [issue](https://github.com/adw96/DivNet/issues)!

Installation
------------

`DivNet` is currently available via GitHub, and can be installed as follows. Please note that you will need to install the latest version of `breakaway` to install `DivNet`; this is included in the following commands:

``` r
# install.packages("devtools")   ## run this line if you do not already have devtools installed
library(devtools)
install_github("adw96/breakaway")
install_github("adw96/DivNet")
library(DivNet)
```

Basic Usage
-----------

Let's analyse the Lee *et al.* dataset using DivNet. We could analyse at the ASV level, but to keep things simple let's analyse the phylum level diversity.

``` r
library(phyloseq)
data(Lee)
divnet_phylum <-  divnet(tax_glom(Lee, taxrank="Phylum"),
                         X = "char",
                         ncores = 4)
divnet_phylum
```

To test if alpha-diversity (by default, Shannon) is equal across the values of the covariate X, we can run

``` r
testDiversity(divnet_phylum)
```

and to plot the DivNet estimates of diversity we can run

``` r
library(ggplot2)
plot(divnet_phylum)
```

Integration
-----------

DivNet integrates with - [phyloseq](https://joey711.github.io/phyloseq/) (for reproducible microbiome analysis) - [breakaway](https://github.com/adw96/breakaway) (for hypothesis testing)

Coming soon
-----------

-   DivNet will hopefully soon be available as a QIIME2 plug-in
-   The same philosophy of DivNet can be applied to estimating UniFrac and other diversity indices that are a function of phylogeny. Watch this space!
