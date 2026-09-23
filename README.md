dbarts
======

[![R-CMD-check](https://github.com/vdorie/dbarts/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/vdorie/dbarts/actions/workflows/check-standard.yaml)

Discrete Bayesian Additive Regression Trees Sampler

A package for R, with C/C++.

Pre-built binaries of the package are built by [CRAN](https://cran.r-project.org/package=dbarts). These can be installed from within R using the typical `install.packages()` mechanism.

Upgrading from 0.9-x
--------------------

`bart` is now the formula-first front door (it was `bart2`, which remains as an alias for one release). The BayesTree-style interface is `bartBT`, with 0.9-x's argument names and defaults; a `bart` call spelled the old way is forwarded to it with a warning for this release. Results differ from 0.9-x. The UPGRADING section of `inst/NEWS.Rd` lists every change that can break existing code.

Features
--------

- Response families via `family =` on `bart()`: `"gaussian"`, `"probit"`, `"logistic"`, accelerated-failure-time survival (`"aft"`), discrete-time survival hazard (`"hazard"`, `"hazard.probit"`, `"hazard.logistic"`), multinomial (`"multinomial"`), ordered categorical (`"ordinal"`), negative-binomial counts (`"nbinom"`), and semicontinuous two-part (`"hurdle.lognormal"`); family objects (`?dbartsFamilies`) carry family settings, such as the residual prior `gaussian(sigma = chisq(3, 0.9))`
- Outlier-robust Student-t residuals (`family = student(df)`)
- Heteroscedastic variance forest (`variance = ~ x1 + x2`)
- Monotonicity (`monotone`), interaction (`interactions`), and block-additive (`blocks`) constraints
- Missing predictor values modeled in place (MIA)
- Priors as objects (`?dbartsPriors`): the tree prior (`tree.prior = cgm(...)`, or DART variable selection with `dart(...)`) and the leaf prior (`node.prior = normal(...)`, or linear and Gaussian-process leaves with `linear(...)` and `gp(...)`)
- Categorical predictors split on level subsets (`factors = "categorical"`); sparse `Matrix::dgCMatrix` and mixed dense/sparse predictor input
- Warm starts from a previous fit (`warm.start`) or XBART-style grow-from-root (`n.grow.sweeps`)
- Reduced-precision residual storage for large problems (`storage = "single"`)
- A sampler object (`dbarts()`) whose response, predictors, offset, and weights can be changed between draws, for use as one step of a larger Gibbs sampler; multi-forest models such as BCF through `dbarts(forests = )`

Steps to install from source:

1. Install development tools for your operating system (dbarts needs R 4.2.0 or newer and a C++20 compiler):
    1. Linux/Unix should already have this installed; if not, use your package manager to install a C/C++ compiler.
    2. macOS: [Xcode](https://developer.apple.com/xcode/resources/)
    3. Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

2. Install the `remotes` package from within R:

```R
install.packages("remotes")
```

3. Run:

```R
remotes::install_github("vdorie/dbarts")
```

For package authors
--------------------

Using dbarts from R needs no special linkage. Packages that call dbarts from compiled code use the flat C API in `inst/include/dbarts/dbarts.h` (`LinkingTo: dbarts`): the sampler is created in R and the C API works on its handle. The C++ headers of 0.9-x, `R_C_interface.hpp` among them, were removed in 1.0-0.
