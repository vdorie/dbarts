# Sparse Unordered Factors

Represents an unordered factor sparsely: entries at the given positions
carry the supplied levels, while every other position implicitly carries
a designated `reference` level. Analogous to
`Matrix::`[`sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html),
whose implicit entry is zero.

Accepted as a predictor by
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
through the x/y interface: a data frame passed as `x.train` may mix
dense numeric/factor columns with sparse ordinal columns
([`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
or `dgCMatrix`) and `sparseFactor` columns in any combination. A
`sparseFactor` column enters as one categorical predictor and bins
bitwise-identically to a dense factor of the same values. **The formula
interface does not accept it**: a bare S4 column cannot survive
[`model.frame`](https://rdrr.io/r/stats/model.frame.html), so the
formula path refuses a `sparseFactor` column explicitly, with a message
to supply it through the x/y interface instead.

## Usage

``` r
sparseFactor(x, levels, reference, i, length)
```

## Arguments

- x:

  The stored entries: a factor, a character vector, or integer level
  codes (which require an explicit `levels`). Missing values are not
  supported. When `i` is omitted, `x` is the complete (dense) vector and
  its non-`reference` entries become the stored ones.

- levels:

  Character vector of factor levels. Defaults to `levels(x)` for a
  factor and `sort(unique(x))` for a character vector, the
  [`factor`](https://rdrr.io/r/base/factor.html) convention.

- reference:

  The level every unstored position carries. Must be an element of
  `levels`; defaults to `levels[1]`, the baseline-contrast convention.
  Storage holds only the non-reference entries, so choosing the **most
  common level** as the reference maximizes the memory win; codes and
  every draw are identical under any choice of reference, so this is a
  storage tuning knob, not a modeling one.

- i:

  Optional integer vector of 1-based positions at which the entries of
  `x` sit, as in
  [`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html).
  When supplied, `x` and `i` must have equal length and `length` is
  required. Positions need not be sorted; duplicates are an error.

- length:

  Total number of observations. Required alongside `i`; defaults to
  `length(x)` for dense input.

## Details

The constructor canonicalizes its input: entries are re-ordered by
ascending position, entries whose level equals `reference` are dropped
(they are the implicit value), and positions are stored 0-based in the
`i` slot with 1-based level codes in the `values` slot. The remaining
slots are `levels`, `reference`, and `length`.

`show` prints the length, the number of stored entries, the level table,
and the reference level. `length` returns the observation count, which
lets a `sparseFactor` be a
[`data.frame`](https://rdrr.io/r/base/data.frame.html) column.

## Value

An object of S4 class `sparseFactor`.

## Author

Vincent Dorie: <vdorie@gmail.com>.

## See also

[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)

## Examples

``` r
f <- factor(c("a", "b", "a", "c", "a", "b"))
sf <- sparseFactor(f, reference = "a")
sf
#> sparseFactor of length 6, 3 stored entries
#>   levels: a, b, c
#>   reference (implicit): a

## equivalently, supply only the non-reference entries and their positions
sparseFactor(c("b", "c", "b"), levels = c("a", "b", "c"), reference = "a",
             i = c(2L, 4L, 6L), length = 6L)
#> sparseFactor of length 6, 3 stored entries
#>   levels: a, b, c
#>   reference (implicit): a
```
