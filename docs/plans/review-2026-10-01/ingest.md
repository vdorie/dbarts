# Review 3, wave 2 - lens: ingest

Pinned tree .claude/worktrees/review3 (01dee4b4), library r3-lib; 0.9-34 from r3-rfit-cranlib.
Probes: scratchpad/r3-ingest-p*.R (p1-p26).

Covered: column kinds through makeCategoricalModelMatrix and the C indicator builder (numeric,
integer, logical, character, factor, ordered, Date/POSIXct/difftime, complex, matrix columns), the
three doors (formula, x/y matrix and frame, dbartsData) on the same data with NA in every kind and
under na.keepPredictors/na.omit/na.exclude and subset, test-level mapping (unused declared levels,
reordered levels, character for factor, numeric for factor), categorical NA position at K = 63, 64,
65, 128 dense and sparse, all-NA and one-level columns, dense vs dgCMatrix vs sparseVector vs
sparseFactor equivalence (uniform and quantile, n.cuts 100 and 5), other sparse Matrix classes,
the cut grid (ties, constant columns, few distinct values, huge/tiny/offset scales, Inf, NaN),
n.cuts validation, setCutPoints, setPredictor introducing NA, save/load/copy/setState round trips
of the grid.

Not covered: hazard/hurdle/multinomial ingestion specifics (forest and family lenses), GP/linear
leaf covariate ingestion, the flat C API's translateSource (bridge lens), Windows.

Excluded as already found: offset() terms, poly()/ns()/cut()-style data-dependent terms (cut()
is the same root as rfit-02), indicator-route test levels (rfit-04, rgen-05; a character column
hits the same bare column-count error), sampler categorical column updates (rgen-03),
offset.test length, vartypes-attribute (TODO).

## ingest-01 - BLOCKER - a constant predictor column makes every stored state unrestorable

Location: src/bartcore/data.hpp DataStore::fillCutsOverRange / fillCutsUniformly (the uniform grid
of a degenerate column repeats one value n.cuts times) against src/bartcore/sampler.hpp
Sampler::setState (`if (state.cutPoints[j][k] <= state.cutPoints[j][k - 1]) return false;`).

Claim: under the default uniform grid, a column that is constant in the training rows (an all-zero
dummy, an all-TRUE logical, a column made constant by subset) gets 100 identical cut points, and
setState refuses any state carrying that grid, so a saved fit cannot predict after loading, and
copy() and setState fail on the live sampler; 0.9-34 restores the same fit.

Probe (p25, p25c, p25o):
```
df <- data.frame(y = rnorm(100), a = runif(100), b = 0)
f <- bart(y ~ ., df, n.samples = 20, n.burn = 20, n.chains = 2, keepTrees = TRUE, verbose = FALSE)
predict(f, df[1:3, ])                      # 40 x 3
f$fit$storeState(); saveRDS(f, fn); f2 <- readRDS(fn)
predict(f2, df[1:3, ])  # Error: state is not consistent with this sampler
# same with bart(y ~ a + c, subset = 1:50) where c is constant on rows 1:50
# dbarts(cbind(a = runif(n), b = 0), y): s$copy() and s$setState(s$state) give the same error
# 0.9-34, bart2 on the same frame: after save/load predict returns 40 x 3
```
The same refusal follows from any uniform grid whose increment is below the column's ulp
(x = 1e15 + 4 * runif(n): 100 cuts, 33 unique, p23) and from a quantile grid whose midpoints
overflow (|x| near 1.7e308, p15).

Why gates missed: test-sampler-degenerate-cuts.R covers the constant column only under
useQuantiles = TRUE, and its comment says the uniform grid places "the single degenerate cut" -
it places n.cuts copies. The engine lens checked constant predictors for NaN and crashes, not for a
state round trip; save/load tests use non-degenerate columns.

Fix: let setState (and installForests' cross-grid check) accept a non-decreasing grid, or build a
degenerate uniform column as the quantile path does (one cut), which moves draws for such fits.

## ingest-02 - BLOCKER - formula fits with a sparse column turn repeated subset rows into missing values

Location: R/data.R dbartsData, formula branch: `pos <- match(rownames(modelFrame), rownames(data))`
before subsetSparseColumn.

Claim: model.frame makes repeated rows' names unique ("1", "1.1"), so the match returns NA for every
repeat and the sparse column (sparseVector, dgCMatrix or sparseFactor) is read as missing on those
rows, while the dense columns and the x/y door keep their values; the fit silently treats them as
NA-routed rows.

Probe (p21, p21b): bootstrap subset `boot <- sample(n, n, replace = TRUE)` (78 repeats of 200),
sparse column s with y = 5 where s > 0:
```
dbartsData(y ~ ., df,  subset = boot)   # sparse column: NA cells 78 (source has none)
dbartsData(y ~ ., dfDense, subset = boot) and dbartsData(df[-1], y, subset = boot): 0 NA
sparseFactor column, same subset: NA cells 78
bart(y ~ ., df, subset = boot) vs the dense-column frame: max |yhat.train diff| 5.24
```

Why gates missed: the sparse formula tests subset with unique indices only.

Fix: carry a row index through model.frame as the sparse-missing marker already rides
(an extra variable such as dbartsRowIndex = seq_len(nrow(data))) and index the sparse columns by it.

## ingest-03 - MAJOR - Date, POSIXct and difftime predictors are refused on the default route

Location: R/utility.R makeCategoricalModelMatrix (`is.numeric(column) || is.logical(column)`,
else "cannot be converted to a predictor").

Claim: a Date, POSIXct or difftime column fails bart(), dbarts(), dbartsData() and xbart() under
the default factors = "categorical" on both the formula and x/y doors, while factors =
"indicators", bartBT and 0.9-34 (bart and bart2) fit it as its numeric value; NEWS does not list
it.

Probe (p1):
```
1.0-0 bart(df, y) / bart(y ~ ., df), z Date|POSIXct|difftime: ERROR: column 'z' cannot be converted to a predictor
1.0-0 bart(df, y, factors = "indicators"): ok, varcount x1, z
0.9-34 bart(df, y) and bart2(y ~ ., df): ok, z used as numeric, predict works
```
(POSIXlt and complex are refused too; 0.9-34 silently dropped complex.)

Why gates missed: no test fits a date-time column; the two builders were never compared on
column kinds.

Fix: in makeCategoricalModelMatrix treat any atomic double/integer column that is not a factor
(Date, POSIXct, difftime, other classed numerics) as ordinal via as.double(unclass(x)), as the C
builder does.

## ingest-04 - MAJOR - a numeric newdata column for a factor predictor is read as 0-based codes

Location: R/utility.R mapFactorColumnsToTrainingLevels (`if (!is.factor(column) &&
!is.character(column)) next`), reached from R/data.R validateXTest for predict, bart(test =),
dbartsData(test =) and $setTestPredictor.

Claim: when training had a factor (or ordered factor) and the new data frame carries that column
as numbers, the numbers pass straight through as level codes 0..K-1; the natural user value,
as.integer(factor), is 1-based, so every row is silently predicted at the next level (only the top
code errors, and only if present). The reverse case (factor test, numeric training) is refused;
0.9-34 errored here.

Probe (p18; levels a, b, c, d with effects 0, 10, 20, 30):
```
predict(f, data.frame(x1 = .5, g = factor(c("a","b","c"), levels(g))))  # -0.1 9.6 19.9
predict(f, data.frame(x1 = .5, g = 1:3))                                  #  9.6 19.9 29.9
0.9-34 bart2 same call: ERROR: LOGICAL() can only be applied to a 'logical', not a 'integer'
```
Related to rgen-03 (the setter column branch), but a different entrance and the opposite offset.

Why gates missed: tests feed factor or character test columns; numeric-for-factor is tested only
in the opposite direction.

Fix: refuse a non-factor, non-character data-frame column whose training column has a level table,
as the opposite direction already is (a raw-code matrix remains the documented code channel).

## ingest-05 - MAJOR - setCutPoints accepts 65535 cuts, whose top bin is the missing code

Location: src/R_interface_bartcore.cpp bartcore_setCutPoints (`if (numCuts > 65535)`), and the
same bound in src/bartcore/sampler.hpp Sampler::setState; data.hpp maxNumCutsRepresentable is
65533 and naCode is 0xFFFF.

Claim: with 65535 cut points a value above the last cut quantizes to code 65535 = naCode, so on a
column with missing values those rows follow the learned missing route and are pooled with the NA
rows; the strictly-increasing check also passes NaN cut points (comparisons with NaN are false),
giving a column that silently cannot split or splits on an unordered grid. Silent, but reachable
only through an explicit setCutPoints of exactly 65535 cuts (or any NaN), hence MAJOR rather than
BLOCKER.

Probe (p5, p4; y = 5 for x < 0.1 or x NA, else 0; cuts seq(0, 0.9, length = nc)):
```
65533 cuts: fit x<0.1 4.98  mid 0     NA 5.01  x>0.9 0.01
65534 cuts: fit x<0.1 4.98  mid 0     NA 5.01  x>0.9 0.02
65535 cuts: fit x<0.1 4.76  mid 0.03  NA 2.19  x>0.9 2.19
s$setCutPoints(c(0.2, NaN, 0.6), 1)  # accepted; s$setCutPoints(NaN, 1) accepted, cor(fit, y) -0.27
```
NaN acceptance is pre-existing in 0.9-34; the 65535 collision is in the new store.

Why gates missed: test-bartcore.R pins 65535 for category codes, not for cut counts; no test sets a
NaN cut.

Fix: bound numCuts by bartcore::maxNumCutsRepresentable in bartcore_setCutPoints and setState, and
refuse NaN cuts (test `!(cuts[i] > cuts[i - 1])`, and isnan on a single cut).

## ingest-06 - MINOR - an infinite predictor value silently voids its column (pre-existing)

Location: src/bartcore/data.hpp DataStore::fillCutsUniformly (range over all non-NaN values).

Claim: one Inf or -Inf makes the uniform range infinite, every cut Inf or NaN, and the column
unsplittable, with no message; bart() then fails with "unable to obtain a starting estimate of
sigma", which names the wrong input. Same in 0.9-34. A range above DBL_MAX (|x| near 1.7e308)
does the same.

Probe (p16, p15): `z <- runif(300); z[1] <- Inf`, y a step in z: uniform cor(fit, y) 0.561,
sigma 1.99 (quantile grid: 0.994, 0.22).

Why gates missed: no test passes non-finite predictors.

Fix: take the range over finite values (an Inf then codes past the last cut) or refuse non-finite
predictors at ingestion, naming the column.

## ingest-07 - MINOR - training refuses integer and logical matrices and non-dgC sparse classes with wrong messages

Location: src/R_interface_bartcore.cpp (`if (!Rf_isReal(slotExpr)) Rf_error("'x' must be numeric")`);
R/data.R dbartsData final `stop("unrecognized 'formula' type; ...")`.

Claim: an integer matrix or vector x is refused as "'x' must be numeric" (it is numeric; 0.9-34
said "x must be of type real"); a logical matrix and a dgTMatrix/dgRMatrix/lgCMatrix x are refused
as "unrecognized 'formula' type", although the test path and setPredictor coerce exactly these
(storage.mode<- and asDgCMatrix) and data-frame integer/logical columns are accepted.

Probe (p12, sparse probe): int matrix: bart/dbarts/bartBT "'x' must be numeric"; logical matrix,
dgTMatrix, dgRMatrix, lgCMatrix: "unrecognized 'formula' type; must be coercible to numeric or a
valid formula object"; dgCMatrix ok.

Why gates missed: tests build x with runif/rnorm and dgCMatrix only.

Fix: coerce integer/logical matrices with storage.mode and other sparse classes with asDgCMatrix in
dbartsData's x/y branches, as validateXTest does.

## ingest-08 - MINOR - factors is silently ignored when dbarts gets a dbartsData object

Location: R/data.R dbartsData (the inherits(formula, "dbartsData") short-circuit warns for data,
test, offset, offset.test, bases, counts but not factors).

Probe (p20): `dbarts(dbartsData(y ~ ., df, factors = "indicators"), factors = "categorical")` runs
with the 4-column indicator design and no warning.

Fix: add factors (when supplied) to the ignored-argument warning.

## ingest-09 - MINOR - bart.Rd's n.cuts length rule is false (pre-existing)

Location: man/bart.Rd \item{n.cuts} ("otherwise must be a vector of length equal to the number of
predictor columns"); R/dbarts.R `data@n.cuts <- rep_len(control@n.cuts, ncol(data@x))`.

Probe (p19): `bart(X2col, y, n.cuts = 5:7)` runs with cuts 5, 6, the 7 dropped silently.
dbartsControl.Rd states the recycling correctly.

Fix: refuse a length other than 1 or ncol in bart, or align bart.Rd with dbartsControl.Rd.

## Checked and found correct

- Formula vs x/y vs dbartsData doors on one frame with NA in numeric, factor, ordered, logical,
  character and sparseVector columns: bitwise identical draws under na.keepPredictors, na.omit and
  na.exclude, with and without a unique subset; omitted rows agree (p6, p17, p20).
- Categorical missing position at K = 63, 64, 65, 128, dense and sparseFactor: NA rows, first and
  last level never conflated, train and predict agree (p11).
- Declared-unused levels, reordered test levels, character test for factor, ordered test given as
  character or re-leveled: mapped by label (p6, p8).
- dgCMatrix vs dense (NA entries, -0, all-nonzero column; uniform and quantile; n.cuts 100 and 5):
  equal to 4e-15 with sigest fixed; dense fit predicts a sparse test bitwise; sparseFactor vs
  factor, both doors and all four train/test pairings, including a test reference level differing
  from training's: identical (p7, p8).
- All-NA columns of every kind refused by name; one-level factor never split; "" level name and an
  NA level (addNA) fit and predict (p14).
- Quantile grid over ties and few distinct values; constant column under quantiles (1 cut,
  restorable); 0.9-34's bus error on that case is gone (p9).
- setPredictor introducing NA into a training-complete column: the route is learned (p13).
- Wide factor (150 levels) on the indicators route: sparse dummy block, NA rows fit, predict equals
  training fit to 2e-14 (p26).
- n.cuts validation: 0, negative, NA, fractional, 65534 and 3e9 refused by name (p19).
