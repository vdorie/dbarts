dbarts
======

Discrete Bayesian Additive Regression Trees Sampler

A package for R, with C/C++.

Pre-built binaries of the package are built by [CRAN](https://cran.r-project.org/package=dbarts). These can be installed from within R using the typical `install.packages()` mechanism.

Steps to install from source:

1. Install development tools for your operating system:
    1. Linux/Unix should already have this installed; if not, use your package manager to install a C/C++ compiler.
    2. OS X: [XCode and gfortran](https://cran.r-project.org/bin/macosx/tools/)
       * If on an ARM processor, install gfortran from [here](https://mac.r-project.org/libs-arm64/)
       * If on an Intel processor, install gfortran from [here](https://github.com/fxcoudert/gfortran-for-macOS/releases/download/8.2/gfortran-8.2-Mojave.dmg)
    3. Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

2. Install the `remotes` package from within R:

```R
install.packages("remotes")
```

3. Run:

```R
remotes::install_github("vdorie/dbarts")
```
