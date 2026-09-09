# Drop Rows With A Missing Response, Keep Missing Predictors

The `na.action` that
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) and
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
default to. It drops the rows whose *response* is missing and keeps the
rows whose predictors are.

## Usage

``` r
na.keepPredictors(object, ...)
```

## Arguments

- object:

  A model frame, as
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html) hands its
  `na.action`; the response is the column the frame's `terms` attribute
  names, and a frame with no response loses no rows.

- ...:

  Ignored, for compatibility with the other `na.action` functions.

## Details

BART routes a missing predictor value down a learned side of each split
rule (“Missingness Incorporated in Attributes”), so an incomplete row is
data rather than a hole and dropping it discards information the model
can use. A missing response, on the other hand, is nothing to fit. This
function is that distinction:
[`na.omit`](https://rdrr.io/r/stats/na.fail.html) would drop both.
rpart's `na.rpart` is the behavioural precedent.

The dropped rows are recorded exactly as
[`na.exclude`](https://rdrr.io/r/stats/na.fail.html) records them - a
named integer vector of class `"exclude"` - so the training fits pad
back to the caller's own row count through
[`naresid`](https://rdrr.io/r/stats/nafns.html): `fitted` and
`residuals` on the fit are as long as the data, `NA` where a row was
dropped.

The base functions keep their usual meaning on these entry points:
[`na.omit`](https://rdrr.io/r/stats/na.fail.html) drops any row with a
missing value anywhere and pads nothing,
[`na.exclude`](https://rdrr.io/r/stats/na.fail.html) drops the same rows
and pads, [`na.fail`](https://rdrr.io/r/stats/na.fail.html) refuses any
missing value, and [`na.pass`](https://rdrr.io/r/stats/na.fail.html)
keeps every row - after which a missing response is an error, since
there is nothing to fit it against.

On the matrix (`x.train`, `y.train`) interface, where there is no model
frame, the `na.action` is applied to a frame standing for the (response,
predictors) pair, so every one of these functions means on that
interface what it means on the formula one.

## Value

`object` with the missing-response rows removed, carrying an
`"na.action"` attribute of class `"exclude"` naming them; or `object`
unchanged when no response value is missing.

## See also

[`na.omit`](https://rdrr.io/r/stats/na.fail.html),
[`na.exclude`](https://rdrr.io/r/stats/na.fail.html),
[`naresid`](https://rdrr.io/r/stats/nafns.html),
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)

## Examples

``` r
df <- data.frame(y = c(1, NA, 3, 4), x = c(1, 2, NA, 4))

# the response row goes, the predictor row stays
model.frame(y ~ x, df, na.action = na.keepPredictors)
#>   y  x
#> 1 1  1
#> 3 3 NA
#> 4 4  4

# and the fit's training values pad back to the data's own length
fit <- bart(
  y ~ x, df,
  n.samples = 20L, n.burn = 20L, n.chains = 1L, verbose = FALSE
)
length(fitted(fit)) == nrow(df)
#> [1] TRUE
```
