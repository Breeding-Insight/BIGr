<!-- badges: start -->
[![R-CMD-check](https://github.com/Breeding-Insight/BIGr/workflows/R-CMD-check/badge.svg)](https://github.com/Breeding-Insight/BIGr/actions)
![GitHub Release](https://img.shields.io/github/v/release/Breeding-Insight/BIGr)
[![Development Status](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![CRAN Status Badge](https://www.r-pkg.org/badges/version/BIGr)](https://cran.r-project.org/package=BIGr)
![GitHub License](https://img.shields.io/github/license/Breeding-Insight/BIGr)
[![codecov](https://app.codecov.io/gh/Breeding-Insight/BIGr/graph/badge.svg?token=PJUZMRN1NF)](https://app.codecov.io/gh/Breeding-Insight/BIGr)


<!-- badges: end -->

=======
# BIGr <img src="https://github.com/user-attachments/assets/2168c801-fcee-4999-b04e-f7b01fed9cfa" align="right" width="250"/>


### Core Genomic Analysis Functions for Breeding Insight

</div>

BIGr is an R package developed by [Breeding Insight](https://breedinginsight.org/) that provides a robust set of functions for analyzing genomic and pedigree data in diploid and polyploid breeding programs. It's designed to streamline the analysis of breeding and genetic data, empowering researchers and breeders to make informed decisions.

## Installation

The stable version of BIGr is now available on CRAN. To install from the R terminal:
```R
install.packages("BIGr")
```


To install the development version of BIGr, you'll need to have `BiocManager` installed. Then, install from GitHub.
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    install.packages("remotes")

BiocManager::install("Breeding-Insight/BIGr", dependencies = TRUE)
library(BIGr)
```
> Note: This GitHub version of BIGr is in development. So, there could be bugs present, and the stable version of BIGr on CRAN should be viewed as more reliable.

## Funding

BIGr development is supported by [Breeding Insight](https://breedinginsight.org/), a USDA-funded initiative based at the University of Florida - IFAS.

## Citation

If you use BIGr in your research, please cite as:

Sandercock, Alexander M., Cristiane H. Taniguti, Josue Chinchilla-Vargas, Dongyan Zhao, Shufen Chen, Meng Lin, Manoj Sapkota, and Breeding Insight Team. 2025. “Breeding Insight Genomics Functions for Polypoid and Diploid Species.” https://github.com/Breeding-Insight/BIGr.

