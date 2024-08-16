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
1. `methodCB`. Given effect sizes and standard errors for predictors of an index trait and a subsequent trait, this function adjusts the statistics for the subsequent trait for selection bias and weak instrument bias through the index trait.
2. `CWBLS`. This function is designed for adjusting collider bias and weak instrument bias to estimate unbiased causal effect between an exposure trait and disease progression trait conditioning on disease incidence, using generalised instrumental effect regression and CWBLS adjustment in bivariate Mendelian randomization.
