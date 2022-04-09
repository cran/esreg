esreg
=====

The goal of esreg is to simultaneously model the quantile and the
Expected Shortfall of a response variable given a set of covariates.

Installation
------------

### CRAN (stable release)

You can install the released version from
[CRAN](https://cran.r-project.org/):

    install.packages("esreg")

### GitHub (development)

The latest version of the package is under development at
[GitHub](https://github.com/BayerSe/esreg). You can install the
development version using these commands:

    install.packages("devtools")
    devtools::install_github("BayerSe/esreg")

If you are using Windows, you need to install the
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) for compilation
of the codes.

Examples
--------

    # Load the esreg package
    library(esreg)

    # Simulate data from DGP-(2) in the paper
    set.seed(1)
    x <- rchisq(1000, df = 1)
    y <- -x + (1 + 0.5 * x) * rnorm(1000)

    # Estimate the model and the covariance
    fit <- esreg(y ~ x, alpha = 0.025)
    cov <- vcov(object = fit, sparsity = "nid", cond_var = "scl_sp")

References
----------

[A Joint Quantile and Expected Shortfall Regression
Framework](https://projecteuclid.org/euclid.ejs/1559872834)
