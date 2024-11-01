# ColliderBias

This is the R package for CWLS (Corrected Weighted Least Squares) and CWBLS (Corrected Weighted Bivariate Least Squares), a designed method to adjust for collider bias in Mendelian Randomisation using summary statistics from GWAS studies.

For details of our methods, please refer to our paper of [CWLS](https://pubmed.ncbi.nlm.nih.gov/35583096/) and [CWBLS](https://pubmed.ncbi.nlm.nih.gov/39445745/).


## Installation
From Rstudio:

```
library(devtools)
install_github("SiyangCai/ColliderBias")
```

## Usage
This package includes two main functions to adjust for collider bias:
1. `methodCB`. This function adjusts the statistics for the subsequent trait for selection bias and weak instrument bias through the index trait.
2. `CWBLS`. This function is designed for adjusting collider bias and weak instrument bias to estimate unbiased causal effect between an exposure trait and disease progression trait conditioning on disease incidence (or other scenarios that interested in estimating causal effect from one risk factor to the outcome, while outcome GWAS is conditioned on another risk factor), using generalised instrumental effect regression and CWBLS adjustment in bivariate Mendelian randomization. 

For both methods, we recommend using as many as independent instruments selected from the original GWAS. This could be performed using PCA, or LD-pruning with high quality. For `methodCB`, one could use one-sample or two-sample design in MR, while using `CWBLS` requires at least exposure and disease GWAS are from two different samples.

A simple example and simulated dataset is included in the package to provide guidance to the users. To run the example, 

```
library(ColliderBias)

# Load the test dataset
data(testData)

# Adjust for collider bias using instrumental effect regression,
# and weak instrument bias using CWLS.
methodCB(testData$dbeta, testData$dse, testData$ybeta, testData$yse, method = "CWLS")

# Find the true causal between an exposure of interest and disease progression.
CWBLS(testData)
```

For more information, please refer to the help page in the R package.
