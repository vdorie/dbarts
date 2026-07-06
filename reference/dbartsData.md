# Discrete Bayesian Additive Regression Trees Sampler Data

Convenience function to create a data object for use with a
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) sampler.

## Usage

``` r
dbartsData(
    formula, data, test, subset, weights,
    offset, offset.test = offset,
    factors = c("categorical", "indicators"),
    missing = c("incorporate", "error"))
```

## Arguments

- formula, data, test, subset, weights, offset, offset.test, factors,
  missing:

  As in [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  Retains backwards compatibility with
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md), so that
  `formula`/`data` can be a
  [`formula`](https://rdrr.io/r/stats/formula.html)/[`data.frame`](https://rdrr.io/r/base/data.frame.html)
  pair, or a pair of `x.train`/`y.train` matrices/vector.

## Value

An object of class `dbartsData`.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
