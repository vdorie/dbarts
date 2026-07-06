# Make Model Matrix from Data Frame

Converts a data frame with numeric and factor contents into a matrix,
suitable for use with
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md). Unlike in
linear regression, factors containing more than two levels result in
dummy variables being created for each level.

## Usage

``` r
makeModelMatrixFromDataFrame(x, drop = TRUE)
makeind(x, all = TRUE)
makeTestModelMatrix(data, newdata)
```

## Arguments

- x:

  Data frame of explanatory variables.

- drop:

  Logical or list controling whether or not columns that are constants
  or factor levels with no instances are omitted from the result. When a
  list, must be of length equal to `x`. Elements correspond to `x`
  according to:

  - vector - single logical

  - matrix - vector of logicals, one per column

  - factor - table of factor levels to be referenced; levels with counts
    of 0 are to be dropped

- all:

  Not currently implemented.

- data:

  An existing
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
  object.

- newdata:

  Test data frame.

## Details

Character vectors are included as factors. If you have numeric data
coded as characters, convert it using `as.numeric` first.

Note that if you have train and test data frames, it may be best to
[`rbind`](https://rdrr.io/r/base/cbind.html) the two together, apply
`makeModelMatrixFromDataFrame` to the result, and then pull them back
apart. Alternatively, save the drop attribute used in creating the
training data and use it when creating a matrix from the test data, as
in the example given below.

Use of these functions is not required when using
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md), `bart2`, or
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md);
they exist to allow the user finer control and to assist with writing
packages that separate the creation of training from test data.

## Value

A matrix with columns corresponding to the elements of the data frame.
If `drop = TRUE` or is a list, the attribute `drop` on the result is set
to the list used when creating the matrix.

## Author

Vincent Dorie: <vdorie@gmail.com>.

## Examples

``` r
iv <- 1:10
rv <- runif(10)
f <- factor(rep(seq.int(3), c(4L, 4L, 2L)),
            labels = c("alice", "bob", "charlie"))
df <- data.frame(iv, rv, f)

mm <- makeModelMatrixFromDataFrame(df)

## create test and train matrices with disjoint factor levels
train.df <- df[1:8,]
test.df <- df[9:10,]
train.mm <- makeModelMatrixFromDataFrame(train.df)
test.mm <- makeModelMatrixFromDataFrame(test.df, attr(train.mm, "drop"))
```
