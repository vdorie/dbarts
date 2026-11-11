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
and `(offset)` columns, with a terms record, so that
[`model.response`](https://rdrr.io/r/stats/model.extract.html),
[`model.weights`](https://rdrr.io/r/stats/model.extract.html) and
[`model.offset`](https://rdrr.io/r/stats/model.extract.html) read them.
Its `"na.action"` attribute is the object's. On an object that keeps its
frame it is made from the frame, a
[`poly()`](https://rdrr.io/r/stats/poly.html) term with its
coefficients; a sparse column (a `sparseVector`,
[`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
or `dgCMatrix`) comes back as stored, at the rows in, and is never made
dense, so the frame holds classes a frame for `lm` could not.

On an object with no frame (a sampler's data, `keep.data = FALSE`, the
matrix interface, an environment) it is decoded from the block where
that is exact: the response is returned in its recorded type and levels
(a factor, a logical, a character vector, a `Surv` response where the
object holds its time and status, the count matrix for a multinomial
one) and named as the formula wrote it (`log(price)`, and `y` where no
formula did), a factor predictor as a factor over its declared levels, a
matrix term such as `poly(x, 2)` or `scale(x)` as one matrix column
again, named as the formula wrote it, with the column names the term
gave, and a column the block holds sparse as a sparse vector, sparse
factor or sparse matrix. The columns come in the formula's order, those
the block holds sparse after the others, and a column no term of the
formula names (a hazard sampler's `period`) after them. The block holds
no more than the values, so a rejoined matrix has none of the attributes
[`poly()`](https://rdrr.io/r/stats/poly.html) or
[`scale()`](https://rdrr.io/r/base/scale.html) give it (their class,
coefficients and centers); the terms record carries them in its
`predvars`, as `lm`'s does. That record is the formula's, each sparse
column listed as a term after the others, so that `lm` on the frame
stops at it as it does on the frame a formula and data give; with no
formula (the matrix interface), or where a hazard sampler added its
`period`, the frame has none and `model.response` does not read it: its
first column is the response.

`nobs` counts the rows in with a positive weight, as the build counted
them, evaluating nothing again, as `nobs` of a fit and of an `lm` fit
do; rows of weight 0 stay in the object, and a BART fit holds them in
its trees.

## Value

A data frame, and a number.

## See also

[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
