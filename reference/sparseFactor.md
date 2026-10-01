# Sparse Unordered Factors

Represents an unordered factor sparsely: entries at the given positions
carry the supplied levels, while every other position implicitly carries
a designated `reference` level. Analogous to
`Matrix::`[`sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html),
whose implicit entry is zero.

Accepted as a predictor by
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
through either interface: a data frame passed as `x.train` may mix dense
numeric/factor columns with sparse ordinal columns
([`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
or `dgCMatrix`) and `sparseFactor` columns in any combination, and a
`formula`'s `data` accepts one too, named like any other column
(including through `.`): a bare S4 column cannot survive
[`model.frame`](https://rdrr.io/r/stats/model.frame.html), so
`dbartsData` lifts it out of `data` first and re-attaches it, row-subset
under `subset` and `na.action`, to the assembled predictor matrix
afterward. A `sparseFactor` column enters as one categorical predictor
and bins bitwise-identically to a dense factor of the same values.

A test frame (`test`, or `newdata` to `predict`) is read the way a
dense-factor frame is: the model's terms are replayed on its dense
columns (so a transformed term such as `log(z)` works and extra columns,
the response included, are dropped) and its `sparseFactor` columns are
used by name. The same holds for a test set (`test`/`x.test`): a
`sparseFactor` test column is recoded over the training level table and
stays resident - through creation and `setTestPredictor` - rather than
densifying at ingestion. `predict` and `getTrees(newdata = )` code it
the same way and then route its rows through the trees off that storage,
materializing no dense matrix of their own.

Mixing a `sparseFactor` (or any sparse ordinal) column with ordinary
dense columns changes how a default starting `sigma` is estimated. The
unmixed case fits an [`lm`](https://rdrr.io/r/stats/lm.html) on the
training design and uses its residual standard deviation; as soon as
*any* column is sparse-backed, the design is not densified for this
purpose and the estimate falls back to the unconditional `sd(y)`
instead - discarding every dense column's information along with the
sparse one's. The fallback warns (class
`dbartsSparseSigmaFallbackWarning`) rather than passing silently. This
is current behavior, not a documented guarantee; supply `sigest`
([`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)) explicitly
if the default matters to you.

## Usage

``` r
sparseFactor(x, levels, reference, i, length)
```

## Arguments

- x:

  The stored entries: a factor, a character vector, or integer level
  codes (which require an explicit `levels`). A missing value is stored
  as an explicit entry, never as the reference level, and a vector whose
  every value is missing needs no levels, as for a
  [`factor`](https://rdrr.io/r/base/factor.html). When `i` is omitted,
  `x` is the complete (dense) vector and its non-`reference` entries
  become the stored ones.

- levels:

  Character vector of factor levels. Defaults to `levels(x)` for a
  factor and `sort(unique(x))` for a character vector, the
  [`factor`](https://rdrr.io/r/base/factor.html) convention.

- reference:

  The level every unstored position carries. Must be an element of
  `levels` (`NA` when there are none); defaults to `levels[1]`, the
  baseline-contrast convention. Storage holds the non-reference entries
  and every missing value, which is always stored, so choosing the
  **most common non-missing level** as the reference maximizes the
  memory win; codes and every draw are identical under any choice of
  reference, so this is a storage tuning knob, not a modeling one.

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

A `sparseFactor` answers the operations a data frame reaches a factor
column with as a [`factor`](https://rdrr.io/r/base/factor.html) does.
`x[i]` subsets by position (positive, negative, zero, logical and
repeated indices) and returns a `sparseFactor` over the same levels,
mapping the stored positions without densifying; `drop = TRUE` drops the
levels no selected row takes, as does `droplevels`, so a vector of
missing values only is left with none. An `NA` or out-of-range index
selects a missing value. `x[i] <- value` and the double-bracket
assignment take labels (as a character vector, factor or `sparseFactor`)
and, as for a factor, store a missing value for `NA` and for a label
that is not a level (with a warning), drop an `NA` index for a value of
length one, and extend the vector past its end with missing values;
`length<-` truncates or pads with missing values; `levels(x) <- value`
renames the levels, the reference with them, merging repeated names as
it does for a factor and dropping a level named `NA`, whose entries
become missing. `c(x, ...)` combines with factors and `sparseFactor`s
over the union of their levels. `as.character`, `as.vector` and `format`
give the level labels, `as.integer` the level codes, and `levels` (read
from the class's slot), `is.na`, `anyNA`, `xtfrm` (so `order` and
`sort`), `unique`, `duplicated`, `rep`, `droplevels`, `summary`, `==`
and `!=` (which refuse two factors with different level sets, as for
factors) behave as for a factor, and `str` prints a factor's line.
`factor`, `as.factor` and `table` read it through those methods, over
the levels present. Together these let a data frame holding one be
subset by row, assigned into, printed and `str`-ed. `rbind` of such
frames returns the column as an ordinary factor; it needs R 4.6.0 or
later, since earlier versions of base R's `match` refuse an S4 object.
On any version, lengthening the first frame by row indexing, as in
`d[rep(seq_len(nrow(d)), 2), ]`, and assigning the other frames' rows
into it binds them and keeps the column a `sparseFactor`, provided their
labels are levels of the first frame's column.

A missing value is an explicit stored entry, and every method answers it
as it does for a factor; `show` counts the missing entries, and a fit
reads them as a missing predictor, bitwise as it reads the same rows of
a dense factor. The one shape without a level is a vector whose every
row is missing; its `reference` is `NA`. A level named `NA`, as
[`addNA`](https://rdrr.io/r/base/factor.html) makes, is not supported: a
missing value is not a level. A character index has no names to match
and is refused. Not supported: `complete.cases` and
[`na.fail`](https://rdrr.io/r/stats/na.fail.html) on a data frame
holding one, which fail in base R for any S4 column, and
[`na.omit`](https://rdrr.io/r/stats/na.fail.html) on such a frame, which
keeps its missing rows (use `d[!is.na(d$f), ]`); `complete.cases` of a
bare `sparseFactor` errors as well, while `na.omit` of one drops the
missing entries as it does for a factor; a fit's `na.action` does handle
them. Also not supported: `c(f, x)` with a factor `f` first, which
`c.factor` answers with a list (put the `sparseFactor` first), `rank`,
`relevel`, `as.numeric` (`as.integer` gives the level codes), `rep_len`
and `rep.int` (`rep` works), and the `exclude` argument of `droplevels`,
which is ignored; `droplevels` on a data frame leaves a `sparseFactor`
column untouched, as base R touches only factor columns, so call it on
the column. `table` drops the declared levels no row takes, where it
keeps them for a factor; use `levels(x)` for the declared table. A data
frame holding one cannot grow by assignment past its last row, since
base R strips the class before the method is reached.

A
[`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
or `dgCMatrix` column subsets in a data frame through Matrix's own
method, but printing such a frame shows a placeholder for the column and
warns of a corrupt data frame; that formatting belongs to Matrix. Use
`as.numeric` on the column to display it.

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
