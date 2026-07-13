# Extract the Predictor Matrix From a Sampler

Returns the predictor matrix as the sampler currently sees it - one
numeric column per predictor, with factor columns as their integer
codes. This is the canonical source for tree replay: it is the same
matrix `getTrees` uses to count training observations through saved
trees.

## Usage

``` r
# S3 method for class 'dbartsSampler'
extract(object, type = "predictors", ...)
```

## Arguments

- object:

  A
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).

- type:

  A character string; only `"predictors"` is supported.

- ...:

  Unused; included for compatibility with the
  [`extract`](https://vdorie.github.io/dbarts/reference/bart.md)
  generic.

## Details

Data frame input is stored columnar in the sampler's data object, so its
`x` slot need not be a plain matrix. `extract` materializes the numeric
coding regardless of the storage form: numeric and logical columns as
doubles, and factor columns - ordered or not - as their integer codes,
`0` to one less than the number of levels, in level order. This is the
same coding a hand-built predictor matrix supplies, so the result can be
fed back to `predict` or used to construct an equivalent sampler.

After predictor mutation (`setPredictor` and relatives) the result
reflects the sampler's current predictor values, not the originally
ingested data. Each mutating call collects the accepted change into the
sampler's data object directly; nothing is maintained per iteration, so
the returned codes always reflect current values, materialized fresh on
every call.

## Value

A numeric matrix with one row per training observation and one column
per predictor, carrying the predictors' column names.

## See also

[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md),
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md),
[`extract`](https://vdorie.github.io/dbarts/reference/bart.md)

## Examples

``` r
df <- data.frame(
  x = runif(50L),
  g = factor(sample(c("a", "b", "c"), 50L, replace = TRUE))
)
df$y <- df$x + 0.5 * (df$g == "b") + rnorm(50L, 0, 0.1)

sampler <- dbarts(y ~ x + g, df, control = dbartsControl(
  n.chains = 1L, n.threads = 1L, updateState = FALSE))
predictors <- extract(sampler, "predictors")
str(predictors)
#>  num [1:50, 1:2] 0.923 0.302 0.36 0.786 0.309 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "x" "g"
stopifnot(all(predictors[, "g"] == as.integer(df$g) - 1L))
```
