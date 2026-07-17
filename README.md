dbarts
======

[![R-CMD-check](https://github.com/vdorie/dbarts/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/vdorie/dbarts/actions/workflows/check-standard.yaml)

Discrete Bayesian Additive Regression Trees Sampler

A package for R, with C/C++.

Pre-built binaries of the package are built by [CRAN](https://cran.r-project.org/package=dbarts). These can be installed from within R using the typical `install.packages()` mechanism.

Features
--------

- Gaussian, probit, logistic, accelerated-failure-time survival, and multinomial response families (`family = "gaussian"/"probit"/"logistic"/"aft"/"multinomial"`)
- Missing predictor values handled in place via MIA (`missing = "incorporate"`)
- DART variable selection prior (`dart = TRUE` or `tree.prior = dart()`)
- Linear and Gaussian-process leaf models (`node.prior = linear(...)` or `gp(...)`)
- Sparse `Matrix::dgCMatrix` predictor input
- Grouped random effects (`rbart_vi`)

See `inst/NEWS.Rd` for release notes.

Steps to install from source:

1. Install development tools for your operating system:
    1. Linux/Unix should already have this installed; if not, use your package manager to install a C/C++ compiler.
    2. OS X: [XCode](https://developer.apple.com/xcode/resources/)
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

Using dbarts from R needs no special linkage. Packages that call into dbarts at the C level should link against the flat C API in `inst/include/dbarts/dbarts.h`; the old C++ ABI, `R_C_interface.hpp`, was removed in 1.0-0.
