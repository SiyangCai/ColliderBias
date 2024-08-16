# ColliderBias

This is the R package for CWLS (Corrected Weighted Least Squares) and CWBLS (Corrected Weighted Bivariate Least Squares), a desgined method to adjust for collider bias in Mendelian Randomisation using summary statistics from GWAS studies.

For detail of our methods, please refer to our paper of [CWLS](https://pubmed.ncbi.nlm.nih.gov/35583096/). The paper of CWBLS is still under review.


## Installation
From Rstudio:

```
library(devtools)
install_github("SiyangCai/ColliderBias")
```

## Usage
This package includes two main functions to adjust for collider bias:
1. `methodCB`. This function adjusts the statistics for the subsequent trait for selection bias and weak instrument bias through the index trait.
2. `CWBLS`. This function is designed for adjusting collider bias and weak instrument bias to estimate unbiased causal effect between an exposure trait and disease progression trait conditioning on disease incidence (or other scenarios that interested in estimating causal effect from one risk factor to the outcome, while outcome trait GWAS is conditioned on another risk factor), using generalised instrumental effect regression and CWBLS adjustment in bivariate Mendelian randomization. 

For both methods, we recommend using as much as independent instruments selected from original GWAS. This could be performed using PCA, or LD-pruning with high quality. For `methodCB`, one could use one-sample or two-sample design in MR, while using `CWBLS` requires at least exposure and disease GWAS are from two different samples.

A simple example and simulated dataset is included in the package to provide guidance to the users.

For more information, please refer to the help page in the R package.
