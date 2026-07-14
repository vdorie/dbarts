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

## Details

### Data frame ingestion

Data frame columns - through the formula interface or a frame passed as
`x.train` - map to predictors as follows:

- numeric and logical columns enter as ordinal predictors;

- unordered factors each become a single categorical predictor under the
  default `factors = "categorical"`, with splits sending a subset of
  levels down each branch; `factors = "indicators"` instead expands each
  factor into binary indicator columns, as previous versions always did;

- ordered factors enter as ordinal predictors, coded numerically by
  level order;

- matrix columns (e.g. from [`poly`](https://rdrr.io/r/stats/poly.html))
  splice in as ordinal predictors, one per column;

- [`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
  and `dgCMatrix` columns (assigned into the frame; they do not survive
  `data.frame(...)` or [`I()`](https://rdrr.io/r/base/AsIs.html)) enter
  as sparse ordinal predictors;

- [`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
  columns enter as a single sparse categorical predictor, coded and
  binned identically to a dense factor of the same values. Only the x/y
  interface accepts them (a bare S4 column cannot survive
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html), so the
  formula path refuses one explicitly).

Factor level tables are retained so that test data are coded
identically; `test` accepts the same column types as `x.train`,
including `sparseFactor` and sparse ordinal columns. A `sparseFactor` or
sparse ordinal test column recodes over the training level table and
stays resident - through creation and `setTestPredictor` - rather than
densifying at ingestion, the same storage tier rule as training.
`predict` and `getTrees(newdata = )` code a test set the same way and
then materialize it to a dense matrix at the call boundary, since those
entry points evaluate already-built trees on resolved values.

Data frame input is stored columnar: no \\n \times p\\ double matrix is
retained. Code that previously reached into the data object's `x` slot
expecting a plain matrix should call
[`extract`](https://vdorie.github.io/dbarts/reference/bart.md)`(sampler, "predictors")` -
or `as.matrix` on the slot - to obtain the numeric predictor matrix,
factor columns as their integer codes.

## Value

An object of class `dbartsData`.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`extract.dbartsSampler`](https://vdorie.github.io/dbarts/reference/extract.dbartsSampler.md),
[`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
