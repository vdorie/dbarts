# Guess Number of Cores

Attempts to guess the number of CPU ‘cores’, both physical and logical.

## Usage

``` r
guessNumCores(logical = FALSE)
```

## Arguments

- logical:

  A logical value. When `FALSE`, an estimate of the number of physical
  cores is returned. When `TRUE`, so-called “logical” cores as also
  included.

## Details

Because of different definitions of cores used by different
manufacturers, the distinction between logical and physical cores is not
universally recognized. This function will attempt to use operating
system definitions when available, which should usually match the CPU
itself.

## Value

An integer, or NA if no clear answer was obtained.

## Author

Vincent Dorie: <vdorie@gmail.com>.
