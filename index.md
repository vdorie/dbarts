# dbarts

[![R-CMD-check](https://github.com/vdorie/dbarts/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/vdorie/dbarts/actions/workflows/check-standard.yaml)

Discrete Bayesian Additive Regression Trees Sampler

A package for R, with C/C++.

Pre-built binaries of the package are built by
[CRAN](https://cran.r-project.org/package=dbarts). These can be
installed from within R using the typical
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html)
mechanism.

## Features

- Response families, via `family =`: `"gaussian"`, `"probit"`,
  `"logistic"`, accelerated-failure-time survival (`"aft"`),
  discrete-time survival hazard (`"hazard"`, `"hazard.logistic"`),
  multinomial (`"multinomial"`), ordered categorical (`"ordinal"`),
  negative-binomial counts (`"nbinom"`), and semicontinuous two-part
  (`"hurdle.lognormal"`)
- Outlier-robust Student-t residuals (`resid.dist = student(...)`)
- Heteroscedastic variance forest (`variance = ~ x1 + x2`)
- Monotonicity (`monotone`), interaction (`interactions`), and
  block-additive (`blocks`) constraints
- Missing predictor values handled in place via MIA
  (`missing = "incorporate"`)
- DART variable selection prior (`dart = TRUE` or `tree.prior = dart()`)
- Linear and Gaussian-process leaf models (`node.prior = linear(...)` or
  `gp(...)`)
- Categorical predictors split on level subsets
  (`factors = "categorical"`); sparse `Matrix::dgCMatrix` and mixed
  dense/sparse predictor input
- Grouped random effects (`rbart_vi`)
- Warm starts from a previous fit (`warm.start`) or XBART-style
  grow-from-root (`n.grow.sweeps`)
- Reduced-precision residual storage for large problems
  (`storage = "single"`)

See `inst/NEWS.Rd` for release notes.

Steps to install from source:

1.  Install development tools for your operating system:
    1.  Linux/Unix should already have this installed; if not, use your
        package manager to install a C/C++ compiler.
    2.  OS X: [XCode](https://developer.apple.com/xcode/resources/)
    3.  Windows:
        [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
2.  Install the `remotes` package from within R:

``` r

install.packages("remotes")
```

3.  Run:

``` r

remotes::install_github("vdorie/dbarts")
```

## For package authors

Using dbarts from R needs no special linkage. Packages that call into
dbarts at the C level should link against the flat C API in
`inst/include/dbarts/dbarts.h`; the old C++ ABI, `R_C_interface.hpp`,
was removed in 1.0-0.
