# The Data of a dbartsData Object

Read the model's variables and the number of rows of a
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
object.

## Usage

``` r
# S3 method for class 'dbartsData'
model.frame(formula, ...)
# S3 method for class 'dbartsData'
nobs(object, ...)
```

## Arguments

- formula, object:

  A
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
  object.

- ...:

  Not used.

## Details

`model.frame` returns the model's variables over the rows in, after
`subset` and the `na.action`, as `model.frame` of an `lm` fit does: the
response first, factors as factors, the weights in a `(weights)` column
and the offsets in the [`offset()`](https://rdrr.io/r/stats/offset.html)
and `(offset)` columns, so that
[`model.response`](https://rdrr.io/r/stats/model.extract.html),
[`model.weights`](https://rdrr.io/r/stats/model.extract.html) and
[`model.offset`](https://rdrr.io/r/stats/model.extract.html) read them.
Its `"na.action"` attribute is the object's. On an object that keeps its
frame it is made from the frame; where the frame holds a sparse column,
which `model.frame` cannot read, it is the frame decoded from the block,
as for an object with no frame.

On an object with no frame (a sampler's data, `keep.data = FALSE`, the
matrix interface, an environment) it is decoded from the block where
that is exact: the response is returned in its recorded type and levels
(a factor, a logical, a character vector, a `Surv` response where the
object holds its time and status, the count matrix for a multinomial
one) and named as the formula wrote it (`log(price)`, and `y` where no
formula did), a factor predictor as a factor over its declared levels,
and a matrix term such as `poly(x, 2)` as one matrix column again, named
as the formula wrote it, with the column names the term gave. The block
holds no more than the values, so a rejoined matrix has none of the
attributes [`poly()`](https://rdrr.io/r/stats/poly.html) or
[`scale()`](https://rdrr.io/r/base/scale.html) give it (their class,
coefficients and centers); every other column is as held.

`nobs` counts every row in, rows of weight 0 included, which a BART fit
holds in its trees; `nobs` of an `lm` fit leaves them out.

## Value

A data frame, and a number.

## See also

[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
