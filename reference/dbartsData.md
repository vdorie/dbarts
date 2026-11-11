# Discrete Bayesian Additive Regression Trees Sampler Data

Convenience function to create a data object for use with a
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) sampler.

## Usage

``` r
dbartsData(
    formula, data, test, subset, weights,
    offset, offset.test = offset,
    factors = c("categorical", "indicators"),
    na.action = dbarts::na.keepPredictors,
    counts = NULL, n.cuts = NULL, useQuantiles = NULL,
    keep.data = TRUE)
```

## Arguments

- formula, data, test, subset, weights, offset, offset.test, factors,
  na.action:

  As in [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  Retains backwards compatibility with
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md), so that
  `formula`/`data` can be a
  [`formula`](https://rdrr.io/r/stats/formula.html)/[`data.frame`](https://rdrr.io/r/base/data.frame.html)
  pair, or a pair of `x.train`/`y.train` matrices/vector.
  `offset`/`offset.test` take the matrix interface's shape from whether
  `counts` is also supplied: with `counts`, an `n` x `K` (`nTest` x `K`)
  numeric matrix installs the per-category shift described under
  `counts` below, and a flat vector is refused; without `counts`, only a
  flat vector (or single number) is accepted, and a matrix is refused.
  `formula` may itself be a `dbartsData` object, in which case it is
  returned as-is; any other argument supplied alongside it is then
  ignored, with a warning naming each one. `weights` are held as given,
  validated here only for length, `NA` and non-negativity; the family's
  own rule applies when a sampler is built from the object, where a
  `"probit"`, `"ordinal"` or `"nbinom"` fit accepts only weights of 0
  and 1, which name the rows in the data set and compose with the
  sampler's active-row mask.

- n.cuts, useQuantiles:

  The maximum number of cut points per predictor, one number for every
  column, one per column, or by column or term name with `.default` for
  the rest, as
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s `n.cuts`
  reads, and whether to place them at the predictors' quantiles, which
  the data object then states; see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  `NULL`, the default, states none. A count is kept as written and read
  again by the sampler it goes to, over the predictors that sampler
  holds: handed to a hazard sampler's `$setData` or to
  `new("dbartsSampler", ...)` with a person-period matrix, `n.cuts = 7`
  gives the predictors 7 and the appended `period` column 100, as
  `dbarts(n.cuts = 7, family = "hazard")` does, and `period` takes a
  count only by its own name. A count one value short, or by name
  without `period`, is kept so a hazard sampler can read it, and any
  other sampler refuses it when it reads it; `dbartsData` itself refuses
  a count no sampler could read. Ignored, with the warning the other
  arguments draw, when `formula` is already a data object.

- keep.data:

  Whether an object built from a formula and a data frame or list keeps
  the frame. `TRUE`, the default, keeps the frame as given, unless the
  build draws random numbers (see ‘Data kept’), the positions of the
  rows kept and what it takes to read the model's variables from them
  again, and holds no predictors, response, weights or offset: see ‘Data
  kept’. `FALSE` keeps the model matrix, the response, the weights and
  the offset in the object, as it was in 0.9-34, and no frame. An object
  built from an environment, with no `data`, or from the matrix
  interface keeps its matrix either way. `formula` being a `dbartsData`
  object takes the argument quietly.

- counts:

  Optional `n` x `K` matrix of non-negative integer counts selecting the
  multinomial (softmax) model, column `k` holding category `k`'s
  successes. A row may carry no trial: it enters no likelihood and still
  receives fitted probabilities. A row's trials, its count total, may
  not exceed \\10^6\\, as each trial costs a Polya-Gamma draw per
  category per sweep. A row with any `NA` cell is a missing response,
  routed through `na.action` exactly as for every other family (see
  [`na.keepPredictors`](https://vdorie.github.io/dbarts/reference/na.keepPredictors.md)).
  It is restricted by `subset` exactly as `weights` is, and it IS the
  response: `y` is derived as its row sums - the trials \\n_i\\ - so no
  separate response need be supplied, and any that is supplied is
  replaced. At least two categories are required. Column names, when
  present, label the categories.
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)`(family = "multinomial")`
  writes this argument from a count-matrix or factor response rather
  than requiring it here. Under `family = binomial` or `"logistic"` an
  `n` x 2 matrix is instead the binomial response, its columns a row's
  successes and failures: `y` is again the row sums, the trials, and the
  offset is the flat shift on the log-odds a binary logistic fit takes;
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) and
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) write
  it from a `cbind(successes, failures)` response. Paired with `counts`,
  a matrix `offset`/`offset.test` is the per-category shift added to the
  raw per-category values BEFORE the softmax - which is the only place a
  shift means anything under one, since the blend is invariant to a
  common per-observation shift and a flat shift would be added after it.
  Only the row-centered part is identified: adding a constant to a whole
  row leaves every reported probability unchanged.

## Details

### Data kept

Built from a formula and a data frame or list, with `keep.data = TRUE`,
the object keeps `data` whole (in its `frame` slot), `rows`, the
positions in it of the rows kept after `subset` and the `na.action`,
read from the model frame as `lm` reads them, so an `na.action` that
records nothing of what it drops keeps the rows it returns (`NULL` when
every row is kept in order; a repeated row, and a `NA` in `subset` that
the `na.action` keeps, appear as they are), `rows.total`, the number of
rows handed over, which `rows` are positions within (`NULL` where `rows`
is), `frame.columns`, the names of its columns, `frame.source`,
`"frame"` (`"environment"` beside the names of an environment the data
were read from, and `NULL` where no data were given), and what reads the
formula, `weights` and `offset` against the frame again. A caller's
variable the formula, `weights` or `offset` names is kept as its value
at the build, a function of no package among them, with the names it
reads, so that a later change to the workspace, a function redefined or
removed included, leaves the object unchanged; a function local to a
caller's function is kept without the rest of that function's variables,
and a package's function by name. A later change to `data` leaves the
object unchanged. The slots `x`, `y`, and where they were given
`weights` and `offset` hold a stand-in that stops, naming
`model.frame(object)`, under `dim`, `length`, `[`, `as.matrix`,
`as.numeric`, arithmetic, `sum`, `mean`, `c`, `rep`, `unique`, `cbind`,
`rbind` and `sd`; `is.na`, `anyNA` and `is.finite` stop too, rather than
answer for values the slot does not hold; `var` and `paste` read it in
R's own code, which gives no way to stop them, and answer `NA` with a
warning or the slot's name; [`attr()`](https://rdrr.io/r/base/attr.html)
of one gives `NULL`, and a `weights` or `offset` that was not given is
`NULL`. Read the data with
[`model.frame`](https://rdrr.io/r/stats/model.frame.html)`(object)` and
[`nobs`](https://rdrr.io/r/stats/nobs.html)`(object)`. Every refusal and
warning of the build is given at the build.

A sampler, and
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)'s
`data`, are made from the block-holding object that `keep.data = FALSE`
gives, which the fitting functions make from a frame-keeping one by
evaluating the formula against the frame again, on the kept rows, with
the caller's variables as they were at the build; the sampler's data
hold no frame. A frame-keeping object loaded in another session does
this at its first fit, so the packages its formula calls must be
installed there. A build whose formula, `weights` or `offset` draw from
R's random-number generator, as `rnorm(n) ~ x` and `y ~ jitter(x)` do,
keeps no frame and holds the values drawn, as `keep.data = FALSE` does,
so that every fit from the object fits the same values; drawing `data`
or `subset` does not count, as the frame and the rows are kept. A term
drawing from a generator outside R's, which leaves R's stream where it
was, is not seen: the object keeps its frame and draws again at each
fit.

### Data frame ingestion

Data frame columns - through the formula interface or a frame passed as
`x.train` - map to predictors as follows:

- numeric and logical columns enter as ordinal predictors, as do classed
  numbers such as `Date`, `POSIXct` and `difftime`, by their numeric
  value (a `POSIXlt` column is refused; convert it with `as.POSIXct`);

- unordered factors each become a single categorical predictor under the
  default `factors = "categorical"`, with splits sending a subset of
  levels down each branch; `factors = "indicators"` instead expands each
  factor into binary indicator columns, as previous versions always
  did - a factor past 100 levels builds this expansion as a sparse
  (`dgCMatrix`-backed) block automatically, with no user-facing
  argument: the tree structure grown is the same as the dense
  expansion's, though leaf-mean draws can differ from it in the last
  bit;

- ordered factors enter as their own predictor kind, coded numerically
  by level order and split at a threshold on those codes; the split grid
  is the midpoints between consecutive *declared* levels, one candidate
  per level boundary, so `n.cuts` does not apply to such a column and
  adjacent levels are always separable;

- matrix columns (e.g. from [`poly`](https://rdrr.io/r/stats/poly.html))
  splice in as ordinal predictors, one per column;

- [`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
  and `dgCMatrix` columns (assigned into the frame; they do not survive
  `data.frame(...)` or [`I()`](https://rdrr.io/r/base/AsIs.html)) enter
  as sparse ordinal predictors. Both interfaces accept them for
  training: a bare S4 column does not survive
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html), so the
  formula path lifts any such column out of `data` before `model.frame`
  runs and re-attaches it, row-subset under `subset` and `na.action`
  alike, to the assembled predictor matrix afterward - named in
  `formula` like any other column, including through `.`, with no
  marker. A sparse name wrapped in
  [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()`,
  [`log()`](https://rdrr.io/r/base/Log.html),
  [`offset()`](https://rdrr.io/r/stats/offset.html), or a `':'`/`'*'`
  interaction is refused by name, as an interaction among dense terms
  already is. `test`/`x.test` accepts them under either interface, since
  a test set is never run through `model.frame`;

- [`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
  columns enter as a single sparse categorical predictor, coded and
  binned identically to a dense factor of the same values. Both
  interfaces accept them, the same pull-out/re-attach the other sparse
  column kinds get on the formula path.

On the formula path, a sparse column re-attaches after every dense term
regardless of where it was named, each in the order its name appears in
`formula` (after `.` expansion) rather than its position in `data`;
`test` is matched to `x.train` by name (or by position when neither
carries names), so the two stay aligned regardless of this reordering.
Factor level tables are retained so that test data are coded
identically, and a factor predictor's level count - ordered or
unordered - is the *declared* level table rather than the levels the
training rows happen to take. A level declared but unobserved in
training therefore keeps its own bin, and test or replacement data
carrying it are accepted at every entrance (creation,
`setTestPredictor`, `predict`, `setPredictor`, `setData`, which reads a
factor of the new data by label against the sampler's table whatever
levels it declares); conversely a value that is not an existing integer
level code is refused at those same entrances, naming the kind, for
ordered and unordered factors alike. An unordered factor is limited to
65535 levels and an ordered one to 65534 - one lower, since an ordered
factor spends a code on the upper bin of its \\K - 1\\ cut grid - and a
wider one is an error naming the column. This holds equally for dense
factors, `sparseFactor` columns, and dense factors inside a mixed
dense/sparse frame. Drop unused levels before fitting
([`droplevels`](https://rdrr.io/r/base/droplevels.html)) if you would
rather they not be modeled. `test` accepts the same column types as
`x.train`, including `sparseFactor` and sparse ordinal columns. A
`sparseFactor` or sparse ordinal test column recodes over the training
level table and stays resident - through creation and
`setTestPredictor` - rather than densifying at ingestion, the same
storage tier rule as training. `predict` and `getTrees(newdata = )` code
a test set the same way and then route its rows through the trees off
that storage, materializing no dense matrix of their own.

Data frame input is stored columnar: no \\n \times p\\ double matrix is
retained. Code that previously reached into the data object's `x` slot
expecting a plain matrix should call
[`extract`](https://vdorie.github.io/dbarts/reference/bartBT.md)`(sampler, "predictors")` -
or `as.matrix` on the slot - to obtain the numeric predictor matrix,
factor columns as their integer codes.

## Value

An object of class `dbartsData`. Its `frame`, `rows`, `frame.columns`
and `frame.source` slots are described under ‘Data kept’. Its `rowNames`
slot records the names of the kept training and test rows, as
`list(train, test)`, or `NULL` when neither is named; a fit carries them
on its outputs as the row names paragraph of
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s ‘Value’
section describes. Its `bases`, `basis.levels`, `basis.written` and
`basis.scale` slots are `NULL` as created: a forest's basis is stated on
the forest (see
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)), and
these slots are written on the data object a sampler of several forests
holds. There `bases` holds each forest's block as the sampler uses it, a
factor's level indicators among them; `basis.levels` records, for each
forest of levels, the levels in the order of the block's columns, and
`NULL` for a forest of numbers or with no basis; and `basis.written` the
text each basis written as code was written as, `NA` for a basis handed
over as a value and for a forest with none, by which a forest beside
these data keeps the held basis (`forest(basis = factor(z))` for one
written `factor(z)`). on the data object a sampler holds it records, for
each forest of numbers that states no `sd`, the standard deviation of
its column over the rows in the likelihood when the sampler was first
created, which that forest's default `sd` was taken from (see `sd` in
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)), and a
sampler created from such an object is handed it, so a copy, a reload
and a sampler made again from its own control, model and data have the
creation's prior whatever column is then in force. A forest that names a
held basis by its text keeps the entry, wherever it stands in the list,
and one given a new value takes the scale of that value; clearing an
entry by hand has the next creation take the scale of the column in
force. Its `missing.seen` slot is `NULL` as created; on the data object
a sampler holds it records, once a predictor has been changed, the
columns that sampler can hold missing values in, as a logical per
predictor column (see ‘Missing values in predictors’ in
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)),
and a sampler created from such an object takes them. Its `status` slot
is `NULL` as created; on the data object an `aft` sampler holds it
records the censoring status in force, 1 for an event and 0 for a
censored row, which `setResponse(status = )` rewrites and a copy or a
reload is created from. Its `useQuantiles` slot is `NA` as created; on
the data object a sampler holds it records whether the cut points are
placed at the predictors' quantiles (see `useQuantiles` in
[`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)),
the value the sampler's control then carries.

Four more slots record what a sampler derives from its data; where
`basis.scale` keeps the value of the sampler's first creation, they
follow what the sampler holds while the data move. Each is `NULL` as
created: `cutPoints`, a list with one entry per predictor column, `NULL`
for an unordered factor and otherwise the column's cut points, strictly
increasing, with their weights, where they are not equally likely, as
the entry's `"mass"` attribute (an increasing sequence whose successive
differences are the weights, as a state's `cutMass` holds them);
`response.range`, the `c(min, max)` of the response the leaf prior's
scale is taken from, `NULL` under a family that reads none (probit,
logistic, ordinal and multinomial), and under `nbinom`
`c(center, center + 1)` of the count fit's center, where a pair of
another width, as a range recorded under another family is, is read by
its first number, the sampler recording it as `c(center, center + 1)`,
and an equal pair is refused; `forest.calibration`,
`c(anchor = , unit = )` of a model of several forests (see
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)), `NULL`
elsewhere; and `standardization`, for a forest of `linear` or `gp`
leaves a list with one entry per forest holding
`list(center = , scale = )`, each named by the leaf covariates where
each one's column name is present, non-empty and the name of no other
column, and otherwise with `columns`, their column positions, beside
them, a `scale` of `NA` for a covariate without spread, and beside them,
for `gp` leaves whose model gives no `lengthscale`, the `lengthscale`
derived with them; `NULL` as a whole, and not an entry per forest, where
no forest's leaves take a covariate. The data object a sampler holds,
`sampler$data`, records them from creation on, and the sampler rewrites
them after every call that changes one. A sampler created from a data
object - by
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
or `new("dbartsSampler", control, model, data)`, and the copy and reload
of a sampler, which are created from its own - uses each value the
object records and derives the others, so `dbarts(sampler$data)` is made
on the sampler's cut points, response range and leaf standardization,
while a data object made afresh from the same rows derives its own. A
sampler's `setData` does the same unless its `updateCutPoints`,
`updateScale` or `updateStandardization` says to derive a value or to
keep the sampler's (see
[`dbartsSampler-class`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)),
so `sampler$setData(sampler$data)`, the object handed back whole, keeps
them. [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md)
takes `cutPoints` for every fold and derives the others in each fold,
from the fold's training rows, whatever the object records. A record
that fits no data - another number of columns, cut points for an
unordered factor, on an ordered factor other than those of its levels,
repeated or not finite, a range that is not two finite numbers, a
standardization whose centers are not finite, whose scales are neither
positive nor `NA`, whose lengthscales are not positive, whose names are
missing, empty or repeated, or whose `columns` are not one position per
center - is refused, naming the slot; setting it to `NULL` derives the
value instead. Cut points out of order are sorted. Which covariates the
leaves take is the model's, so a standardization is matched to them by
its keys, names to leaf columns whose names could key one, the others
deriving, `columns` where a leaf column's name could not, and one with
neither by position where the counts agree: an entry for a covariate the
leaves do not take is not read, a covariate without one is derived, and
a `lengthscale` the model supplies is used in place of the record's.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`model.frame.dbartsData`](https://vdorie.github.io/dbarts/reference/model.frame.dbartsData.md),
[`extract.dbartsSampler`](https://vdorie.github.io/dbarts/reference/extract.dbartsSampler.md),
[`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)

## Examples

``` r
set.seed(3)
n <- 60L
df <- data.frame(
  x = runif(n),
  g = factor(sample(c("a", "b", "c"), n, replace = TRUE))
)
df$y <- df$x + 0.5 * (df$g == "b") + rnorm(n, 0, 0.1)

## the formula interface; 'g' enters as one categorical predictor
data <- dbartsData(y ~ x + g, df)
data
#> An object of class "dbartsData"
#> Slot "y":
#> <'y' is read with model.frame(o)>
#> 
#> Slot "x":
#> <'x' is read with model.frame(o)>
#> 
#> Slot "varTypes":
#> [1] 0 1
#> 
#> Slot "x.test":
#> NULL
#> 
#> Slot "weights":
#> NULL
#> 
#> Slot "weights.test":
#> NULL
#> 
#> Slot "offset":
#> NULL
#> 
#> Slot "offset.test":
#> NULL
#> 
#> Slot "n.cuts":
#> [1] NA NA
#> 
#> Slot "sigma":
#> [1] NA
#> 
#> Slot "missing":
#> [1] "incorporate"
#> 
#> Slot "na.action":
#> NULL
#> 
#> Slot "rowNames":
#> NULL
#> 
#> Slot "missing.seen":
#> NULL
#> 
#> Slot "cutPoints":
#> NULL
#> 
#> Slot "response.range":
#> NULL
#> 
#> Slot "forest.calibration":
#> NULL
#> 
#> Slot "standardization":
#> NULL
#> 
#> Slot "response.type":
#> [1] "numeric"
#> 
#> Slot "response.n.levels":
#> [1] NA
#> 
#> Slot "response.levels":
#> NULL
#> 
#> Slot "bases":
#> NULL
#> 
#> Slot "basis.levels":
#> NULL
#> 
#> Slot "basis.scale":
#> NULL
#> 
#> Slot "basis.written":
#> NULL
#> 
#> Slot "counts":
#> NULL
#> 
#> Slot "counts.as.factor":
#> [1] NA
#> 
#> Slot "offset.category":
#> NULL
#> 
#> Slot "offset.category.test":
#> NULL
#> 
#> Slot "status":
#> NULL
#> 
#> Slot "useQuantiles":
#> [1] NA
#> 
#> Slot "testUsesRegularOffset":
#> [1] NA
#> 
#> Slot "frame":
#>             x g            y
#> 1  0.16804153 a  0.177577428
#> 2  0.80751640 c  0.760531131
#> 3  0.38494235 b  0.993684793
#> 4  0.32773432 a  0.265626355
#> 5  0.60210067 a  0.697894896
#> 6  0.60439405 a  0.519816750
#> 7  0.12463344 a  0.009659326
#> 8  0.29460092 a  0.244830643
#> 9  0.57760992 b  1.014011813
#> 10 0.63097927 b  1.005901046
#> 11 0.51201590 a  0.335204293
#> 12 0.50502391 c  0.506964701
#> 13 0.53403535 c  0.651273625
#> 14 0.55724944 b  1.127721920
#> 15 0.86791949 b  1.276516351
#> 16 0.82970869 a  0.841611818
#> 17 0.11144915 c -0.011854548
#> 18 0.70368836 b  1.330332559
#> 19 0.89748826 a  0.905008414
#> 20 0.27973255 b  0.857943974
#> 21 0.22820188 a  0.244923711
#> 22 0.01532989 a -0.034990768
#> 23 0.12898156 b  0.499862880
#> 24 0.09338193 a  0.052575161
#> 25 0.23688501 b  0.621328776
#> 26 0.79114741 c  0.745552247
#> 27 0.59973157 c  0.703889559
#> 28 0.91014771 b  1.424777615
#> 29 0.56042455 b  1.032684790
#> 30 0.75570477 b  1.387674186
#> 31 0.37917189 c  0.320629695
#> 32 0.37328098 c  0.481533368
#> 33 0.17029064 a  0.168516403
#> 34 0.45330732 b  0.924658735
#> 35 0.25841396 b  0.806339550
#> 36 0.33626595 a  0.152153052
#> 37 0.88958303 c  0.883719720
#> 38 0.20194630 c  0.120276483
#> 39 0.57918604 b  1.272290346
#> 40 0.20763203 b  0.589972076
#> 41 0.28146879 a  0.344054654
#> 42 0.78628120 a  0.874243472
#> 43 0.17301935 b  0.697341157
#> 44 0.57074752 c  0.628454551
#> 45 0.41928296 c  0.279732169
#> 46 0.26762217 c  0.232456352
#> 47 0.04780944 c  0.026730705
#> 48 0.10349305 a  0.133342322
#> 49 0.31403146 a  0.441749469
#> 50 0.80064106 b  1.251978519
#> 51 0.22932470 c  0.253527651
#> 52 0.21299844 a  0.275760383
#> 53 0.87710091 c  0.727522408
#> 54 0.99322196 b  1.454490179
#> 55 0.84424702 a  1.054661032
#> 56 0.91043655 c  0.964278766
#> 57 0.47126973 c  0.405456058
#> 58 0.22441841 a  0.254249106
#> 59 0.12781466 c  0.132526675
#> 60 0.27968351 a  0.272371030
#> 
#> Slot "rows":
#> NULL
#> 
#> Slot "rows.total":
#> NULL
#> 
#> Slot "frame.columns":
#> [1] "x" "g" "y"
#> 
#> Slot "frame.source":
#> [1] "frame"
#> 
#> Slot "build":
#> $formula
#> y ~ x + g
#> <environment: 0x5fd7495c25c8>
#> 
#> $weights
#> NULL
#> 
#> $offset
#> NULL
#> 
#> $factors
#> [1] "categorical"
#> 
#> $coding
#> $coding$x
#> $dense
#> $dense[[1]]
#> numeric(0)
#> 
#> $dense[[2]]
#> factor()
#> Levels: a b c
#> 
#> 
#> $sparse
#> NULL
#> 
#> $map
#> [1] 1 2
#> 
#> $numObservations
#> [1] 0
#> 
#> $columnNames
#> [1] "x" "g"
#> 
#> attr(,"class")
#> [1] "dbartsMixedMatrix"
#> attr(,"term.labels")
#> [1] "x" "g"
#> attr(,"varTypes")
#> [1] 0 1
#> attr(,"factor.levels")
#> attr(,"factor.levels")[[1]]
#> NULL
#> 
#> attr(,"factor.levels")[[2]]
#> [1] "a" "b" "c"
#> 
#> attr(,"terms")
#> ~x + g
#> attr(,"variables")
#> list(x, g)
#> attr(,"factors")
#>   x g
#> x 1 0
#> g 0 1
#> attr(,"term.labels")
#> [1] "x" "g"
#> attr(,"order")
#> [1] 1 1
#> attr(,"intercept")
#> [1] 1
#> attr(,"response")
#> [1] 0
#> attr(,".Environment")
#> <environment: base>
#> attr(,"predvars")
#> list(x, g)
#> attr(,"dataClasses")
#>         y         x         g 
#> "numeric" "numeric"  "factor" 
#> 
#> $coding$missing.seen
#> [1] FALSE FALSE
#> 
#> 
#> $nobs
#> [1] 60
#> 
#> 

## the same data through the backwards-compatible x/y interface, with a
## test set that is coded against the training level table
train <- df[1:40, ]
test  <- df[41:60, c("x", "g")]
data <- dbartsData(train[, c("x", "g")], train$y, test)

## build a sampler on it, rather than letting dbarts() ingest again
sampler <- dbarts(data, control = dbartsControl(
  n.chains = 1L, n.threads = 1L, n.samples = 10L, n.burn = 10L,
  updateState = FALSE), n.trees = 25L)
samples <- sampler$run()
str(samples$test)
#>  num [1:20, 1:10] 0.283 0.78 0.663 0.536 0.478 ...
```
