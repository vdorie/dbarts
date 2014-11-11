dbarts
======

Discrete Bayesian Additive Regression Trees Sampler

A package for R, with C/C++.

Pre-built binaries of the package are available on http://cran.r-project.org/web/packages/dbarts/index.html. These can be installed from within R using the typical `install.packages()` mechanism.

Steps to install from source:

  1. Install development tools for your operating system:

    1. Linux/Unix should already have this installed
    2. OS X:
        1. Xcode (https://developer.apple.com/xcode/downloads/)
        2. gfortran (http://coudert.name/software/gfortran-4.8.2-Mavericks.dmg)
    3. Windows: Rtools (http://cran.r-project.org/bin/windows/Rtools/)

  2. Install the devtools package from within R:

    `install.packages("devtools")`

  3. Run:

    `install_github("dbarts", "vdorie")`