# Mixed Dense/Sparse Predictor Container

The container a data frame's predictor columns assemble into whenever
any column is sparse-backed (a
[`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
column, or a
[`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)/`dgCMatrix`
column assigned into the frame): dense columns stay a per-column list
and sparse columns stay one shared `dgCMatrix`, rather than
materializing one \\n \times p\\ dense matrix.
[`makeModelMatrixFromDataFrame`](https://vdorie.github.io/dbarts/reference/makeind.md)
returns one whenever `x` carries a sparse column;
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)'s
x/y interface (and so
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s and
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s) returns
one under the same condition and stores it in the resulting `dbartsData`
object's `x`/`x.test` slot. A fit or sampler made from a formula or a
data frame holds one as its `data@x` whether or not a column is sparse.
The formula interface accepts a sparse column too, named like any other
(including through `.`): a bare S4 column cannot survive
[`model.frame`](https://rdrr.io/r/stats/model.frame.html), so it is
lifted out of `data` first and re-attached, row-subset under `subset`
and `na.action`, to the assembled predictor matrix afterward.

## Usage

``` r
# S3 method for class 'dbartsMixedMatrix'
dim(x)
# S3 method for class 'dbartsMixedMatrix'
dimnames(x)
# S3 method for class 'dbartsMixedMatrix'
x[i, j, drop = TRUE]
# S3 method for class 'dbartsMixedMatrix'
as.matrix(x, ...)
# S3 method for class 'dbartsMixedMatrix'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `dbartsMixedMatrix`.

- i, j:

  Row and column indices, as for a matrix. Supplying `j` densifies (see
  ‘Details’); `i` alone keeps the container.

- drop:

  Logical; as for a matrix, only consulted when `j` is supplied.

- row.names:

  `NULL`, the default, for a data frame's default row names, or a vector
  of row names.

- optional:

  Unused.

- ...:

  Unused.

## Details

`dim` and `dimnames` - and so the base `nrow`, `ncol`, and `colnames`,
which are defined in terms of them and need no method of their own -
read like a matrix's; row names are always `NULL`. `x[, j]` (or
`x[, "name"]`) returns the indexed/named column(s) as an ordinary
numeric vector or matrix. `x[i, ]` keeps the container, row-subset.
`x[i]` with `j` omitted also row-subsets, unlike a plain matrix's
single-index linear indexing. `as.matrix` materializes the whole object
to a dense \\n \times p\\ numeric matrix (a factor column's codes,
zero-based) and is the escape hatch for any operation the container
itself does not support.

`as.data.frame` gives the predictors as the model sees them, one column
per column of the container, named by `colnames(x)`: numeric columns as
numeric, and a factor column as a factor over its declared levels, which
keeps a level no row takes. For a fit made from a formula of plain
variables, `predict` and `setTestPredictor` take that data frame as new
data. A formula that transforms or expands a variable (`I(z^2)`,
`factor(i)`, a matrix term, `poly(z, 2)`) names columns that are not the
variables it reads, so new data for it is a data frame of those
variables, not this one.

What is **not** supported: matrix-style subassignment
(`x[i, j] <- value` errors). `is.matrix(x)` is `FALSE`, and a base
function that does not dispatch on the class, such as `cbind`, does not
see a matrix. Code that needs a matrix should call `as.matrix(x)` first.

## Value

`dim` returns an integer vector of length two. `dimnames` returns a
length-two list, the first element always `NULL`. `[` with `j` supplied
returns a numeric vector or matrix, matching `drop`; with `j` missing it
returns another `dbartsMixedMatrix`. `as.matrix` returns a plain numeric
matrix, and `as.data.frame` a data frame.

## See also

[`makeModelMatrixFromDataFrame`](https://vdorie.github.io/dbarts/reference/makeind.md),
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md),
[`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)

## Examples

``` r
df <- data.frame(a = 1:5)
df$b <- Matrix::sparseVector(c(1, 2), c(1, 3), length = 5)
mm <- makeModelMatrixFromDataFrame(df)

class(mm)
#> [1] "dbartsMixedMatrix"
dim(mm)
#> [1] 5 2
colnames(mm)
#> [1] "a" "b"
mm[, "a"]
#> [1] 1 2 3 4 5
as.matrix(mm)
#>      a b
#> [1,] 1 1
#> [2,] 2 0
#> [3,] 3 2
#> [4,] 4 0
#> [5,] 5 0
```
