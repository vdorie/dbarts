# sparse-formula-audit

status: SURVEY DELIVERED (2026-09-08); recommends one front-door idiom,
        no code written
agent: opus
rng: neutral (no engine change proposed)
budget: survey + recommendation

## Question

dbarts exports a `sparseFactor` class and accepts sparse columns in a
data frame through the x/y interface, but the formula path refuses
them. Before choosing an interface, this is what the established R
packages actually do when most predictors are ordinary data-frame
columns and one or more are sparse - either a sparse matrix block or a
categorical variable with thousands of levels.

## Base R's type rule

A data frame really can hold a sparse matrix. Only direct assignment
works - `d$S <- M` succeeds where `data.frame(S = M)` and `I(M)` both
refuse - and once assigned, `nrow`, row subsetting and `na.omit` all
handle it correctly, dropping the same rows from the sparse block as
from everything else. Printing the frame is the one thing that breaks,
with a "corrupt data frame" warning.

None of that helps, because `model.frame` refuses the column on type
before any of it runs: "invalid type (S4) for variable 'S'". A list
column fails the same way. The only multi-column object `model.frame`
carries is a plain base matrix, which is dense by construction - that
is the mechanism behind `poly()` and `ns()`, and it is why those work
and a sparse block cannot.

There is a partial exception worth knowing about: an S4 class whose
data part is a basic type (`contains = "integer"`) passes the type
check and reaches the model frame with its slots intact. It does not
survive row subsetting of the data frame, which silently strips it back
to a bare integer vector, so it is not a foundation to build on.

## Survey

**Matrix.** The user writes an ordinary formula over an ordinary data
frame and calls `sparse.model.matrix` instead of `model.matrix`.
Sparsity is a property of the output, never the input: a 1500-level
factor over 2000 rows comes out at 0.26 MB against 18 MB dense. It
calls `model.frame` first, so a sparse column in `data` dies with the
same S4 error. Test data are matched the base way, by replaying the
stored level tables.

**MatrixModels.** `model.Matrix(f, data, sparse = TRUE)` is the same
formula over the same dense data frame with a switch; `glm4` uses it.
Same refusal for a sparse column in `data`.

**lme4.** This is the canonical answer for a thousand-level factor. The
user writes

```r
lmer(y ~ x + (1 | subject), data = d)
```

with `subject` an ordinary dense factor column. `model.frame` carries
it densely; lme4's own term parser then builds the sparse random-effect
matrix from it, one nonzero per row. Nothing sparse is ever in the data
frame. Test data are matched against the stored level table, with
unseen levels an error unless explicitly allowed.

**glmnet.** Matrix-only: `x` may be a `dgCMatrix`, there is no formula.
glmnetUtils adds one and deliberately bypasses `model.frame` by
default, walking the formula's additive terms and building each term's
block itself - because the terms object is quadratic in the predictor
count, and because `model.matrix`'s dropped baseline level is wrong for
a penalized fit. Its `sparse = TRUE` sends each block through
`sparse.model.matrix` and column-binds the results, and it stores a
level table per term for test matching. It still does not accept a
sparse column in `data`.

**ranger.** The formula-and-data-frame path is dense throughout. Sparse
input is available only through `ranger(x = M, y = y)`, all-or-nothing,
and only for `dgCMatrix`. Handing either path a data frame with a
sparse column dies with a row-count mismatch inside the data frame's
own element-replacement method, because ranger flattens the frame.

**xgboost.** Version 3 dropped the formula interface outright. The
front door is `xgboost(x, y)`, where `x` is a matrix, a data frame of
numeric/integer/logical/factor columns with factors treated natively as
categorical, or a sparse matrix. The documentation states the either/or
plainly: categorical features are supported for data-frame input only,
not on sparse matrices. Test data are matched by feature name.

The consensus is uniform. Sparsity is either an output of the formula
path or an input to the matrix path, never both. No package carries a
sparse block through a formula, and the thousand-level factor is
handled by keeping the factor dense in the frame and letting a term
construct the sparse representation downstream.

## Recommendation

**A `sparse()` formula term whose argument is evaluated outside the
model frame.** The user writes `y ~ a + b + sparse(g)`, where `g` is
either an ordinary factor column - which dbarts turns into a
`sparseFactor` itself - or a name resolving to an already-sparse object
in the formula's environment. dbarts strips the term from the formula
before building the model frame and evaluates its argument separately,
which is exactly what the existing `forest()` term machinery already
does with its bases, so the sparse object never has to enter either the
data frame or the model frame and no fight with the type check arises.
This is lme4's idiom for the many-level factor and glmnet's for the
block, expressed once; it keeps the user's data frame printable and
ordinary; and test-data matching comes free, since a level table
replays exactly as `xlev` does today.

Runners-up. First, **pull out and re-attach**: detect S4 columns in
`data`, drop them before `model.frame`, then re-attach by name, aligned
by matching the model frame's row names back into the original frame's.
This is exact under `subset` and `na.action` together - data frames
guarantee unique row names - and it makes the existing class usable
from a formula unchanged. It loses on `y ~ .`, whose dot cannot see
columns that were removed before expansion, and it leaves the user with
a data frame that will not print. Second, **keep refusing**, sending
sparse input to the x/y interface. That is literally what every package
surveyed does, and it costs nothing; it is ranked last only because
dbarts already carries the mixed-column machinery, so the refusal is a
gap in its own surface rather than a missing feature.

## Ingestion work

Under the recommendation, ingestion does not need the pull-out step at
all. The formula walker gains one more recognized term, evaluates its
argument in the formula's environment, and appends the resulting
columns to the dense block by name after the model frame is built. The
row alignment machinery that a `forest()` basis already uses under
`subset` covers the sparse columns unchanged.

## Deferred checks (2026-09-09)

dec-B100 chose the pull-out/re-attach runner-up over the `sparse()` term:
"if there is a sparse matrix in a data frame, why do we need to identify it
as sparse in the formula?" This section runs the checks this survey left
open against a `d$S <- M` column (a 10-row frame, `M` a 3-column
`dgCMatrix`), to confirm nothing else in a data frame's ordinary machinery
trips on one before the pull-out ever reaches it.

**order, row subsetting.** `d[order(d$b), ]` reorders `S` right along with
every other column, rownames included; `d[i, ]` for an arbitrary row index
does the same. Both are exact - this is the same base mechanism `merge`
(next) and `model.frame`'s own row-dropping (below) rest on.

**merge.** `merge(d, d3, by = "a")` realigns `S` correctly to the merged row
order (verified by matching each output row's `S` back to its source row by
`a`), because merge's own row reordering is implemented as ordinary `[`
indexing over the whole frame, the same mechanism as plain subsetting.

**split.** `split(d, g)` produces one sub-frame per group, each with `S`
intact and correctly reduced to that group's rows (again plain `[`
indexing per group).

**rbind.** `rbind(d[1:5, ], d[6:10, ])` FAILS: "cannot set length of
non-(vector or list)". `rbind.data.frame` tries to reassign into the
combined `S` column directly rather than going through `[<-`, and an S4
matrix is neither. Irrelevant to ingestion (dbartsData never rbinds a
train and a test frame together) but worth knowing: a sparse column
survives every reshaping this package's formula path performs, and fails
the one it does not.

**head/tail.** `head(d, 3)` slices `S` to its first 3 rows correctly - `[`
again.

**str.** Prints the column as `Formal class 'dgCMatrix' ... with 6 slots`
and its raw `@i`/`@p`/`@x` - accurate, if not pretty; matches the survey's
"printing... breaks" finding, which is about `print`/`format.data.frame`
specifically, not `str`.

**saveRDS round trip.** `identical(readRDS(saveRDS(d)), d)` holds; ordinary
R serialization has no S4-column special case to get wrong.

**model.frame's row-dropping under subset.** `model.frame(a ~ b, data =
d5, subset = b > -100, na.action = na.omit)` where `d5` carries a sparse
`S` column (untouched by the formula) and NAs in `a`: the returned frame
drops exactly the NA rows, `attr(mf, "na.action")` records their original
positions, and `rownames(mf)` is the surviving subset of `rownames(d5)` -
confirms `pos <- match(rownames(modelFrame), rownames(data))` (step 12)
recovers the right original rows regardless of whether `subset` ran,
`na.action` ran, or both did, since both apply through the same row-name
mechanism a sparse column never has to enter.

Net: every check but `rbind` passes silently, and `rbind` never runs on
these frames. The pull-out/re-attach implementation (step 10) follows
directly - `data[!isSparse]`/`data[isSparse]` for the split, `match()`
against `rownames()` for the re-attach.
