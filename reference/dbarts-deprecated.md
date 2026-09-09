# Deprecated and Removed Functions in dbarts

Names kept reachable for one release past their removal or rename, each
pointing at its successor. They are removed in dbarts 1.1-0.

`bart2` is the former name of the modern front door, which is now
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md); a `bart2`
call is forwarded unchanged, with a warning the first time in a session.

`rbart_vi` and its methods fit a grouped random-intercept BART model.
That model now lives in the stan4bart package (`stan4bart::stan4bart`),
whose prior on the group spread differs from the one dbarts drew from,
so a refit moves the results rather than reproducing them. These stubs
error.

## Usage

``` r
rbart_vi(...)

# S3 method for class 'rbart'
predict(object, ...)
# S3 method for class 'rbart'
extract(object, ...)
# S3 method for class 'rbart'
fitted(object, ...)
# S3 method for class 'rbart'
residuals(object, ...)
```

## Arguments

- object:

  An `rbart` fit object saved by dbarts 0.9-x.

- ...:

  Ignored; the stubs do not run. `bart2`'s own arguments are
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s, copied
  rather than restated - `formals(bart2)` is `formals(bart)` - so a
  caller reading the defaults programmatically still finds them.

## Value

`bart2` returns whatever
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) returns. The
`rbart_vi` stubs do not return.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md)
