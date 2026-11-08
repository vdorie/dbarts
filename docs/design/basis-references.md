# A forest's basis names a column of the data

Status: PROPOSED, 2026-11-07, revised 2026-11-08 to the rulings dec-B674 to dec-B678.
Follows dec-B663 (a forest's basis is a reference to a column of the data), dec-B666 (the frame kept on the
data object), dec-B670 (the whole frame kept, and beside it as little derived from it as can be), dec-B671
(what a fit of several forests keeps under `keepFits = FALSE`), dec-B674 (`keep.data = FALSE`, the frame kept
by default), dec-B675 (the response, weights and offsets made from the frame), dec-B676 (`o@x` a stand-in
that stops when read), dec-B677 (`model.frame(o)` and `nobs(o)` the accessors; `extract()` keeps to samples)
and dec-B678 (a sampler's data holds the model frame and what its forests need, the setters writing into it).

## In short

Two kinds of data object, one class and one set of slots:

- **A user's object**, `dbartsData(formula, data)`, keeps the caller's frame whole, the positions of the rows
  it kept and the record of how its variables are coded. It stores no predictors, response, weights or
  offsets: `o@x`, `o@y`, `o@weights` and `o@offset` are stand-ins that stop, naming `model.frame(o)`.
  `model.frame(o)` makes the model's variables from the frame; `nobs(o)` counts its rows.
  `keep.data = FALSE` keeps no frame and keeps the block, response, weights and offsets in their slots, as
  today; bartCause and stan4bart pass it.
- **A sampler's data**, `s$data`, made where a sampler is handed its data (`dbarts(o)`, `dbartsSpec(o)`,
  `bart(o)`, `$setData(o)`, `new("dbartsSampler", data = o)`), holds the model frame: one column for each
  predictor column the engine reads, a factor kept a factor, the response as the sampler reads it, the
  weights, the summed offset and each forest's basis. The caller's frame is dropped. The setters write into
  it; `model.frame(s$data)` reads it; `s$data@x` and its kin are stand-ins as on a user's object.

The engine's input is unchanged: before each call that hands the bridge a data object, R makes today's
layout from the frame, sharing its columns, so nothing under `src/` reads a new shape. A forest's basis
given in a forests list beside a user's object is read against its frame as `dbarts(f, d, forests = )` reads
it; `dbartsData(bases = )` goes.

Run on a shipped build of main of 2026-11-07, at 20,000 rows, 14 numeric columns and a factor of 26 levels, one
thread; today's main differs on these paths in the setters' argument checks and the bridge's standardization
code, so a rerun on today's main is owed before the build:

- **Setters**: `$setResponse` 138 us today, about 8 us more; `$setPredictor` of one column 2.65 ms today,
  the same within the timer; one iteration of 75 trees is 3.0 ms. On the matrix interface a column's
  setter gets 0.45 ms faster (it no longer copies the whole matrix). A whole-matrix `$setPredictor` would
  pay 0.22 ms more (3 percent of the call) if R split the matrix into columns; a split in the bridge avoids
  it.
- **Saved sizes** stay within 0.001 MB of today: a sampler 4.55 MB, a `bart` fit kept with its sampler
  11.82 MB, bartCause's fits 2.18 and 2.34 MB, stan4bart's fit 0.46 MB, as long as stan4bart's bundle
  keeps what it keeps today (see Costs dec-B678 did not weigh). A saved sampler on main writes its data
  twice, once through its engine pointer; 0.9-34 wrote it once, 2.32 MB.
- **Cost**: about 480 new lines of R and 360 lines changed where R reads a sampler's slots, about 15 in the
  bridge, up to about 740 test lines moved, seven help pages, one NEWS entry, about 8 lines in bartCause, and
  in stan4bart 4 lines plus its bundle (open).

## 1. What exists

Run on main, 2026-11-07 and 2026-11-08.

- **A formula's `forest()` term** (`bart()`, `dbarts()`), `forest(x2, basis = factor(z))`: code, its names
  found in `data` and then where the formula was written; `predict(fit, newdata)` reads `z` from `newdata` by
  name.
- **A forests list beside a frame** (`dbarts(f, data = d, forests = )`): a column of `data` is that column and
  hides a caller's variable of the same name; a value is full length and cut by `subset` and the na.action.
- **The matrix interface**: a value, full length, cut the same way.
- **`dbartsData(bases = list(NULL, v))`**: values by forest position. bartCause's `bcf()` uses it, then
  `dbarts(o, forests = )` with knobs only.
- **`dbartsSpec(o, forests = )`**: no frame, so a basis is the caller's value at the object's kept length.
  Handed a sampler's own data, a declaring list replaces its bases whole, with no message.
- **`$setForestBasis`** by label or position; **`$setData`** refused on a sampler of several forests.
- **A data object** stores the predictors in `x` (a columnar container on the formula path, a factor one
  column; a plain matrix or a `dgCMatrix` on the matrix interface), `y`, `weights`, `offset`, `x.test` and
  their kin, and no frame. In a session its columns are the caller's frame's columns, shared; on the matrix
  interface `x` is a copy of the caller's matrix.
- **A sampler's data** is that object, rewritten by the setters: `$setResponse` writes `y`; `$setPredictor`
  of a column replaces that column of the container, and of a plain matrix copies the whole matrix first; a
  whole-matrix `$setPredictor` replaces `x` with a plain matrix; `$setOffset`, `$setWeights`,
  `$setTestPredictor` and `$setForestBasis` write their slots. The engine keeps no predictor block: it reads
  the data object when a sampler is created, at `$setData`, at a model change and when a state is put back.
- **A saved sampler** writes its data twice, 4.55 MB at this size: once as its `data` field, once as the
  protected value of its engine pointer. 0.9-34's saved sampler was 2.32 MB, its pointer empty.
- **0.9-34** had no forests, no bases and no `dbartsSpec()`; its data object held the model matrix,
  factors spread into indicator columns, in `x`.

## 2. The design

### 2.1 Rows: what is numbered and how

A user's object records the number of rows handed over and the positions of the rows kept. It is the index
`dbarts()` already builds to cut a basis written beside a frame, so a value and a code basis are cut by one
index (12 rows; rows 3 and 9 missing a predictor):

- no `subset`, default na.action: 1 to 12 (`na.keepPredictors` keeps rows with a missing predictor)
- `subset = 2:10`: 2 to 10; `subset = -1`: 2 to 12; `subset = z == 1`: 5, 9, 10, 11, 12
- `subset = c(10, 2, 5, 7)`: in that order; `c(2, 2, 5, 7, 8)` keeps the repeat
- a character `subset` by row name: the positions of those names
- `na.action = na.omit`: rows 3 and 9 go; after `subset = 2:10` too: 2, 4, 5, 6, 7, 8, 10

An integer vector, 4 bytes a row; empty when every row is kept in order. A sampler's data needs none: its
frame holds the rows in, in order. The internal paths that cut or expand rows (the hurdle's two parts and
its merge, the hazard's person-period expansion, `xbart()`'s folds) act on a sampler's data and cut its frame.

### 2.2 Columns: what a user's object keeps

The whole frame, as given, when `data` is a data frame or a list (dec-B670). In a session it costs nothing:
R shares the columns, and a later `d$z[1] <- 9` copies `z` alone, the object keeping the old one. Saved, at
20,000 rows and 14 columns, the frame is 1.19 MB. A column added to `d` later is not in the object. A name
that is no column is looked up where `forest()` was called. Set aside: the frame less the formula's columns
(a basis naming a predictor then falls through to the caller), numbers into the caller's frame (a reorder or
a dropped row misreads them in silence) and named columns only (dec-B666).

`keep.data = FALSE` keeps no frame (dec-B674); the object keeps the block, `y`, `weights` and `offset` as
today, and the frame's column names, so a basis named by one of them is refused rather than read from the
caller (call 6).

### 2.3 A user's object, under dec-B670 and dec-B674 to dec-B677

Kept, beside the frame and the rows:

- **The coding record**: the block with no rows, that is its column names, `varTypes`, term labels, factor
  levels (or indicator levels and the dropped columns), the frozen terms (dec-B185: `poly()` coefficients,
  `scale()` centers, a caller's value), and for each column whether it held a missing value. 750 bytes where
  the formula names only columns; a caller's vector named in the formula rides at full length, 8 bytes a row.
- **The response's coding**: `response.type`, `response.n.levels`, `response.levels`.
- **The frozen expressions** of the response, `weights` and `offset`, so a column stays a name and a
  caller's value is held as its value (dec-B675).
- **`na.action` and `rowNames`**: the rows the na.action dropped, numbered within the subset, and the names
  model.frame gave the kept rows ("2", "2.1", "5" under a repeated subset), which the positions alone do not
  give (call 19).
- **From the caller's other arguments**, kept as built: `x.test`, `weights.test`, `offset.test`, `n.cuts`,
  `sigma`, `offset.category` and its test twin, `useQuantiles`, `testUsesRegularOffset`, `missing`.

Not stored: the predictors, the response, the weights and the offsets. `o@x`, `o@y`, `o@weights` and
`o@offset` hold a small classed stand-in whose `dim()`, `nrow()`, `colnames()`, `[` and `as.matrix()` stop
with a message naming `model.frame(o)`; `attr(o@x, ...)` gives `NULL` without a word, R reading attributes
without asking the object, and the help says so (dec-B676).

`model.frame(o)` returns, as `model.frame()` of an `lm` fit does, the model's variables over the rows in,
after `subset` and the na.action, factors kept as factors, the weights in a `(weights)` column and the offsets
in `offset()` and `(offset)` columns, so `model.response()`, `model.weights()` and `model.offset()` read them;
`nobs(o)` is its row count (dec-B677). It is made in the build's order: the frozen terms evaluated on the whole
frame, the variables cut to the kept positions, then checked against the recorded levels. Run, the predictor
columns this makes, coded, equal the built block in every case tried: plain columns, `y ~ .`, logical and
matrix columns, `offset()` terms, `I(x2^2)`, `log(x2) + scale(k)` and `poly()` under `subset` (one with a
repeated row), factors and ordered factors with a `subset` that empties a level, character columns in either
coding, `na.omit`, `na.exclude`, `na.pass`, the default `na.keepPredictors` with a missing factor value, a
value or a missing value that only a dropped row holds, a caller's scalar, a caller's vector under `subset`
and under `na.omit`, and a sparse factor. Replaying `predict()`'s path on the kept rows instead fails on a
caller's vector once a row is dropped (call 20). Time: 0.0014 s at 20,000 rows and 15 columns under a
subset of 15,000, against 0.004 s for the build.

`makeTestModelMatrix(data, newdata)` (help `makeind`) reads the coding record and its missing-value flags,
so on a default object it needs no block and gives what it gives today.

Saved, a user's object is the frame and the record: 2.22 MB at 20,000 rows with a 26-level factor, where
frame and block together would be 4.5 MB.

**An object built on an environment** (call 7) holds the model frame made at the build, the shape a
sampler's data holds (2.4), and the environment's variable names.

### 2.4 The sampler's data, under dec-B678

**Where it is made.** At the step `dbarts(o)`, `dbartsSpec(o)`, `bart(o)`, `rbart_vi(o)`, `xbart(o)`, the
hurdle's entry, `new("dbartsSampler", data = o)` and `$setData(o)` share, the resolution. It is the object
`dbartsSpec()` returns as its `data`, so stan4bart's saved bundle holds it. Handed a sampler's data again,
the step takes it as it is.

**What it holds**, in one slot, `model.frame`, a data frame carrying `terms` (call 24 for the shape):

| column | holds | made from |
|---|---|---|
| each predictor column the engine reads | a numeric column as doubles; a factor as a factor of the recorded levels, a character column coded to one; a matrix-valued term (`poly()`, `scale()`, a matrix variable) split into its columns; under `factors = "indicators"` each indicator column; a sparse column as a `sparseVector` or `sparseFactor` | a user's object: `model.frame(o)`, coded; a block-holding object: the block's columns, shared |
| the response, first | the response as the sampler reads it: the coded doubles; for `aft` a two-column (time, status) matrix; for `multinomial` the count matrix | `model.response()`, coded (call 23) |
| `(weights)` | the weights, absent when there are none | `model.weights()` or the slot |
| `(offset)` | the one offset the sampler adds, the formula's `offset()` terms summed with the argument | `model.offset()` or the slot (call 27) |
| `(offset.category)` | the multinomial category offset | the slot |
| `(basis.k)` | forest k's basis, a number or a two-level factor, for each forest that has one | the resolved basis, cut by the rows (call 28) |

Its row names are the kept rows' names, and its `na.action` attribute what the na.action dropped, as
`model.frame()` gives them, so those two slots go on a sampler's data (call 32). Kept beside it as records,
none a copy of a row: the coding record (column names, `varTypes`, term labels, levels), the response's
coding, `basis.levels` and `basis.scale`, `n.cuts`, `sigma`, `useQuantiles`, `missing`, and what a sampler
writes (`missing.seen`, `cutPoints`, `response.range`, `forest.calibration`, `standardization`). The caller's
frame and the row positions are dropped; the frame's column names are kept for call 6.

The test rows stay where they are: `x.test`, `weights.test`, `offset.test`, `offset.category.test` (call 29).
They come from `test`, not from the data, and no model frame holds them.

**Slots that go on a sampler's data**: `y`, `x`, `weights`, `offset`, `bases`, `counts`, `status`,
`offset.category`, `na.action` and the training part of `rowNames`. `x`, `y`, `weights` and `offset` hold the
stand-ins of 2.3, naming `model.frame(s$data)`.

**The engine's input.** Each call that hands the bridge a data object (creation, `$setData`, a model change,
`xbart()`'s data handle, `$setState`'s raw predictors) is handed today's layout, made by R from the frame
just before: the predictor columns as the columnar container (a double or factor column shared, not copied;
an integer or logical one turned to doubles), the sparse columns assembled into the container's `dgCMatrix`,
the record's attributes set, and the response, weights, offset and bases as the slots the bridge reads. Run:
0.006 ms to make the container at 20,000 rows and 15 columns. The cut points and codes are what the engine
derives from it, as today. Nothing under `src/` changes for this (call 30).

**At each setter and door**, what it writes and what it costs per call. Times at 20,000 rows, 14 numeric
columns and a 26-level factor, 75 trees, one thread, where one iteration is 3.0 ms:

| setter or door | writes into the frame | today | under dec-B678 |
|---|---|---|---|
| `$setResponse(y)` | the response column | 138 us | about 146 us: the column write costs 33 us against 25 us for the slot's |
| `$setResponse(y, status)` (aft) | both columns of the response | as above | as above |
| `$setPredictor(v, j)` | column j | 2.65 ms | the same: 31 us to write the column against 33 us today |
| the same, a factor column | the factor as given, matched to the recorded levels | as today | as today; a factor with the recorded levels is stored as is (2 us) |
| the same, on the matrix interface | column j | 3.10 ms | about 2.65 ms: today's 0.45 ms copy of the whole matrix goes |
| `$setPredictor(x)`, whole | every predictor column | 6.8 ms | about 6.8 ms with the split in the bridge (one copy of the block, about 0.08 ms, in place of today's `as.double` copy, 0.21 ms); 7.0 ms if R split it (0.43 ms, call 31) |
| `$setPredictor(x, forceUpdate = "partial")` | the installed rows of column j | as today | as today |
| `updatePredictorPerObservationJointly` | each sampler's column | as today | as today |
| `$setOffset(o)` | `(offset)`, removed for `NULL` | | one column write, as `$setResponse` |
| `$setWeights(w)` | `(weights)`, removed for `NULL` | | one column write |
| `$setCounts` (multinomial) | the response's count matrix | | the trials are made from it where read, not stored |
| `$setForestBasis` | `(basis.k)` | | one column write |
| `$setTestPredictor`, `$setTestOffset` | `x.test`, `offset.test`, as today | | unchanged |
| `$setActiveRows` | nothing (the mask is the sampler's) | | unchanged |
| `$setData(o)` | a new frame, made at the resolution | 7.7 ms | about 7.7 ms from a block-holding object; from a user's object 1.4 ms more to make the frame; from a matrix-interface object 0.08 ms to split it (0.43 ms in R) |
| `$setState` | nothing; the raw predictors are read from the frame | | unchanged: the container is made in 0.006 ms |
| `$copy()` | the copy shares the frame in a session | 12 ms | unchanged |
| save and reload | the frame is written; a reload creates the engine from it | | the same size (2.5) |
| `xbart()`'s folds | the folds cut the container made once at its start | | unchanged |

The R-side reads of the old slots (`data@x` and its kin) move to internal accessors that read the frame:
the predictors' container, the response, the weights, the offset, the bases, the row count.

**A factor column given to `$setPredictor`**: a factor or a character vector of the column's labels, matched
by label, a level the column does not declare refused by name, as today; a factor whose levels are the
recorded ones in order is stored as given, another is recoded to them. Numbers and codes are refused, as
today: a code is no label. Under `factors = "indicators"` the frame holds the indicator columns, and a 0/1
vector given to one of them is taken, as in 0.9-34 and today (run: both take it, and both refuse the factor
by name); call 25.

**The matrix interface.** `dbarts(x, y)` makes the frame of the matrix's columns, one column each under the
matrix's column names, the response first, and the weights and offset as above. An integer matrix is turned
to doubles, as today. A `dgCMatrix` gives one `sparseVector` column each (call 33).

**An unnamed matrix** (dec-B653). The frame's columns need names, so they take placeholder names, and the
coding record says the predictors are unnamed: `colnames()` of the predictors stays `NULL`, an unnamed matrix
given later is read in column order with no message, a named one given to an unnamed sampler still warns
"'x' and 'test' are not both named", and `model.frame(s$data)` shows the placeholder names. `y ~ X` on an
unnamed `X` keeps the names model.frame gives, `X.1` and `X.2`, as today.

**The test rows, a frozen state, and the memory rule.** A sampler's state holds the trees, the variances and
the chains' scalars, about 0.016 MB here, and no data (state never carries model). It is unchanged: a saved
sampler is its control, its model, its data with the frame, and its state.

### 2.5 What a block-holding object and its sampler hold, and what the consumers read

Under dec-B674 an object built with `keep.data = FALSE`, and the matrix interface's object, hold the block,
`y`, `weights` and `offset` in their slots, as today, and the frame's column names. A sampler made from one
holds the model frame made from those slots: the container's columns shared in a session, a matrix's split
into columns. A user who keeps a matrix-interface object beside its sampler then holds the predictors twice
in a session, n times p times 8 bytes (2.4 MB at 20,000 by 15); `dbarts(x, y)` with no object kept does not.

What bartCause and stan4bart read, against the examination of what each should do:

- **bartCause** passes `keep.data = FALSE` at its two data call sites (three objects), so its saved fits keep
  `data.rsp` as today, and its reads of `nrow()`, `colnames()`, `x[, z]`, `as.matrix(x)`, `@y`, `@weights` and
  its write of `@x.test` stay. Its `bcf()` gives the treatment forest `basis =` with the full-length
  `factor(z)` it already builds, needing calls 1, 2 and 6 as recommended (run: draws identical to the data
  door, with and without `subset`). New under dec-B678: it reads `attr(fit$fit$data@x, "term.labels")` off
  the samplers its fits keep (`fit.trt$fit`, `fit.rsp$fit`) in two places; on a stand-in that gives `NULL`
  without a word, so both move to the term labels of `model.frame(fit$fit$data)`'s terms. About 8 lines.
- **stan4bart** passes `keep.data = FALSE` in `glFormula`, so its fit's `bartData` reads (`@x`, `@varTypes`,
  `factor.levels`, `makeTestModelMatrix()`) stay. New under dec-B678: it reads `sampler$data@y` and
  `sampler$data@weights` when it installs a response range, and `data@offset` on the specification's data at
  a restore (3 lines, to `model.response()`, `model.weights()` and `model.offset()` of the model frame); and
  its saved bundle strips the block off the specification's data (`@x` and `@x.test` cut to no rows) and
  splices `bartData@x` back at a restore, which cannot be written so once the data holds a model frame. The
  ways out are in Costs dec-B678 did not weigh.
- The examination's other recommendations (bartCause's test frame in place of a hand-built `x.test`, its
  `plot_support` on a fit with missing responses, stan4bart's hand-written slots moved to arguments) are
  untouched by this design.

### 2.6 Saved sizes

Run on a shipped build of main; the dec-B678 column projects each sampler's data with the frame of 2.4 in
place of `x`, `y`, `weights`, `offset` and `bases` (the frame's saved size, from the same values).

| object | today | under dec-B678 |
|---|---|---|
| a `dbarts()` sampler, formula, 20,000 rows, 14 numeric and a 26-level factor | 4.55 MB (data 2.23, the engine pointer 2.22) | 4.55 MB (data: 2.220 to 2.219) |
| a `bart` fit with `keepTrees = TRUE`, its sampler kept, 50 draws | 11.82 MB | 11.82 MB |
| the same with `keepTrees = FALSE` | 7.46 MB | 7.46 MB |
| a matrix-interface sampler, 20,000 by 15 | 4.81 MB (data 2.36) | 4.81 MB (the frame of columns 2.34) |
| stan4bart, 750 rows, 25 BART columns, one chain | 0.46 MB; bundle 0.050, its data 0.007 | 0.46 MB if the bundle strips the frame's predictors; 0.50 MB, bundle 0.088, if it keeps the specification's data whole |
| bartCause `bartc()`, 750 rows, 25 covariates, bart response | 2.18 MB | 2.18 MB (each kept sampler's data within 0.001 MB) |
| the same, bcf response | 2.34 MB | 2.34 MB |

0.9-34 saved the same sampler at 2.32 MB, writing its data once.

### 2.7 The doors

- **`dbartsData()`**: keeps the frame, the row record and the coding record (2.3); `keep.data = FALSE` the
  block; `bases =` gone; its formula still refuses a `forest()` term, pointing to a forests list.
- **`dbartsSpec(o, forests = )`** and **`dbarts(o, forests = )`**: a basis is read against `o`'s frame, then
  where `forest()` was called, built on `o`'s rows by the function `dbarts()` uses, so
  `dbartsSpec(dbartsData(f, d), forests = L)` and `dbarts(f, d, forests = L)` draw the same. A value has the
  length of the rows handed over (call 1). `dbarts(o)` takes a declaring list (call 2).
- **A sampler's data beside a declaring list**: refused where it carries basis columns (call 3); a name of a
  column of the caller's frame refused where it carries none (call 6).
- **`$setForestBasis`** and **`predict`**: as built.
- **The matrix interface**: a value, full length, cut by the record.
- **Fit and reload**: `dbartsSpec()` and `dbarts(o)` write the model record `dbarts(f, d)` writes.

### 2.8 Data objects that carry bases, and setData

Once `bases =` goes, a data object carries basis columns only where a sampler of several forests made them,
`s$data` and what is made from it; a sampler of one forest makes none.

- **`bases = list(NULL)`**, one forest with no basis, has no path: a record of fewer than two forests stays
  refused at creation (call 17).
- **`$setData` on a sampler of one forest**, given another sampler's data that carry basis columns: refused
  before anything is stored, the text saying they are the data of a sampler of several forests and pointing
  to the data object they were made from or to `dbarts()`. Its own data handed back is taken.
- **`$setData` on a sampler of several forests**: refused whatever it is given, as today.
- **`dbarts(s$data)`, `dbartsSpec(s$data)` and a reload**: the basis columns and `basis.levels` re-create the
  forests.

Three texts that advise the data door (the `$setData` refusal, `dbarts()`'s refusal of a declaring list, the
refusal of a model of one forest with a basis) are reworded.

### 2.9 What else this needs

- **The slots**: `frame` (the caller's data, a user's object only), `rows` (its positions, `NULL` when every
  row is kept in order), `model.frame` (a sampler's data and an environment's object), `frame.columns` (the
  frame's or the environment's names) and the coding record. Where `frame` or `model.frame` is set, `x`, `y`,
  `weights` and `offset` hold stand-ins; where neither is, they hold the block and the values. The class's
  validity checks each kind. The help describes the slots as read-only (dec-B554) and shows no slot write.
- **The record's check**: the positions lie within the rows handed over; a model frame's row count, its
  response's length and each basis column's agree. Sampler creation and the validity call it.
- **`dbartsData(o)`** handed a data object returns it as it is.
- **A data object saved by 0.9-34** already fails on main ("no slot of name "bases""); nothing changes there
  (call 22).
- **A comment to fix**: the one above the row numbers `dbarts()` hands to the data object says a code basis is
  built only on the rows kept; it is read over every row and then cut.

## 3. The landscape

Run where installed (mgcv 1.9.4, nlme 3.1.171, lme4 2.1.0, brms 2.23.0); stochtree and bcf read from source.

- `lm`, `glm`: a varying slope is a term; names found in `data`, then the formula's environment; the fit keeps
  its model frame (`model = TRUE`) and `model.frame(fit)` returns it, weights in `(weights)` and the offset
  argument in `(offset)`.
- mgcv `s(x2, by = z)`, lme4 `(1 + z | g)`: the pre-built object (`gam(fit = FALSE)`, `lFormula()`) holds the
  model's columns; `predict` needs them.
- nlme: the fit keeps the whole data frame unless `keep.data = FALSE`.
- brms: by-variables must be in `data` or `data2`.
- stochtree, bcf (Hahn, Murray and Carvalho): values only. BART and bartMachine: no basis.

## 4. What a user meets

- A formula user, a forests-list user beside a frame, a matrix-interface user of `dbarts()`: no draw moves.
- A reader of `s$data@x`, `s$data@y`, `s$data@weights` or `s$data@offset`: an error naming
  `model.frame(s$data)`, which holds them. In 0.9-34 `s$data@x` held the model matrix and `s$data@y` the
  response; code that runs the BART sampler inside a larger sampler and reads them meets this.
- A reader of `o@x` and its kin on an object built from a frame: the same, naming `model.frame(o)`.
- A caller of `$setPredictor` on the matrix interface: a column's setter is faster.
- treatSens uses the matrix interface and reads no slot; nothing moves.
- bartCause and stan4bart: 2.5.
- A `dbartsSpec()` caller with a basis: a name that is a column now finds it; where it is also a caller's
  variable, the column wins, as at `dbarts()`: the one call whose draws move. A kept-length value on an object
  that dropped rows is refused.
- Anyone who saves a data object built from a frame: it holds the frame in place of the block, about the same
  size where the frame is mostly the model's columns, larger where it holds others.

## 5. The cost

- R, about 480 new lines and 360 changed, nothing in `inst/include/dbarts/dbarts.h`:
  - the slots, the stand-ins on both kinds and the validity, about 40
  - `dbartsData()` keeping the frame, the row and coding records, `keep.data`, about 45
  - `bases =` leaving the exported function, about 20
  - the resolution making the sampler's frame: from a user's object in the build's order, from a
    block-holding object or an environment's, the response coded, the offsets summed, matrix terms split,
    the basis columns, about 110
  - the engine view made from the frame before each bridge call, about 40
  - the internal accessors, about 40, and the reads moved to them: about 360 lines that read or write `@x`,
    `@y`, `@weights` or `@offset` today (182, 73, 64 and 41)
  - the setters writing the frame (`$setResponse`, the three `$setPredictor` paths and the joint update,
    `$setOffset`, `$setWeights`, `$setCounts`, `$setForestBasis`, the hazard expansion, the hurdle parts),
    about 70
  - `model.frame()` and `nobs()` methods, about 23
  - `dbartsSpec()` and `dbarts(o)` resolving against the frame, calls 1, 2, 3, 6 and 7, about 85
  - `makeTestModelMatrix()` reading the record, about 5
- `src/`: about 15 lines in the bridge, an entry that splits a double matrix into column vectors for the
  whole-matrix setter and `$setData` (call 31). No engine change, no draw moves.
- Tests: the 83 calls of `dbartsData(bases = )` in 24 files become forests lists drawing what they drew (45
  fixtures), pins of the forest door (36) or oracles on the list's value (2). Reads of the old slots, about
  740 lines in up to 90 files (472 `@x`, 135 `@y`, 59 `@weights`, 74 `@offset`, by name, an upper bound), move
  to `model.frame()`, the internal accessors or `keep.data = FALSE`. New, about 250 lines: the sibling
  identity, the row positions, the frame's contents per family and interface, each setter writing its column
  with draws identical to today's, the stand-ins' errors, a factor column's setter, indicator coding, an
  unnamed and a sparse matrix, the saved size, `makeTestModelMatrix()` on a default object, the refusals of
  calls 3, 6 and 7.
- Help: `dbartsData` (the slots, `keep.data`), `dbartsSampler-class` (the data field), `dbartsSpec`,
  `dbarts`, `makeind`, a page for the `model.frame()` and `nobs()` methods, and the forests help where it
  names the data door. Four design notes that name the data door (forest-kind-by-class, public-surface,
  written-surface, multiplier-combiner) change at landing.
- NEWS, against 0.9-34, one entry: "A data object made by dbartsData() from a formula and a data frame keeps
  the frame and the positions of the rows it keeps, and a sampler's data holds the model frame, the setters
  writing into it. In both, the x, y, weights and offset slots, which held the model matrix and the response
  in 0.9-34, no longer hold them: model.frame(object) returns them and nobs(object) counts the rows. A forest
  given to dbarts() or dbartsSpec() beside a data object finds its basis in the frame. keep.data = FALSE keeps
  the matrix and the response in the object instead, as before."
- bartCause: about 8 lines and its suite. stan4bart: 4 lines, its bundle (open), about 5 test lines and its
  suite.

## Calls made

Ruled: call 4 by dec-B674, call 8 by dec-B675, call 9 by dec-B676, call 10 by dec-B677 (`model.frame(o)` and
`nobs(o)`, no `extract()` method for a data object), call 5 by dec-B678 (the caller's frame dropped where a
sampler is handed its data, the sampler's data holding the model frame and what its forests need), the
former call 11 by dec-B671.

Taken as the design's calls, for the maintainer to mark, each following the rulings' direction (one shape
for a user's object and a sampler's data, values read from the model frame, no copy of the caller's data in
a saved object):

1. **The length of a basis given as a value**: the length of the rows handed over, at every creation door; a
   value already cut to the kept rows refused where rows were dropped. One rule for one argument, as at
   `dbarts(f, d, forests = )`. Cost: within `dbartsSpec()`'s lines. Alternative: either length at
   `dbartsSpec()`, two rules for one argument.
2. **`dbarts(o, forests = )` declaring a basis**: taken, as `dbartsSpec()` takes it, so both read one data
   object alike; bartCause's `bcf()` needs it. Cost: about 25 lines. Alternative: `dbarts()` refusing,
   bartCause moving to `dbartsSpec()` and `new("dbartsSampler")`.
3. **A sampler's data carrying basis columns, beside a list declaring a basis**: refused at both doors, naming
   the data object to start from; a list of knobs only keeps the record. A sampler's model frame holds the
   model's variables and `(basis.k)` columns, not the caller's columns a basis expression names, so a new
   basis has nothing to be read against. Cost: about 10 lines. Alternatives: the list replacing the bases
   whole, where a plain forest in it loses its basis; or any list refused, knobs too.
6. **A basis naming a column beside an object with no caller's frame** (a sampler's data of one forest, a
   `keep.data = FALSE` object): a name that is a column of the model frame it holds is read from there, the
   rows in, as a column of `data` is read at `dbarts()`; a name of another column of the caller's frame (kept
   in `frame.columns`) is refused, pointing to the object built with its frame; any other name is the
   caller's value. A transformed predictor's source (`x2` under `log(x2)`) is no column of the model frame and
   is refused, not taken from the caller. Cost: about 10 lines. Alternatives: every frame name refused, the
   model frame not read; refusing any declaring list beside a sampler's data, which leaves a
   `keep.data = FALSE` object reading from the caller; every name read from the caller, a column silently
   replaced by a variable of its name.
7. **An environment given as `data`**: the object holds the model frame made at the build, as a sampler's
   data does, and the environment's variable names; a basis naming a model variable is read from that frame,
   another of the environment's names is refused, pointing to a data frame. Every value is read at one time,
   the build. Cost: about 10 lines. Alternatives: the block kept as today with every environment name refused;
   or the environment kept and read at the resolution, a basis read later than the predictors.

Made here, in building dec-B678:

23. **The response column holds the response as the sampler reads it**, coded doubles (0/1 for a binary
    response, codes for an ordinal one), the coding record keeping its type and levels, because
    `$setResponse` takes numbers and the column must hold what it was given. Alternative: the type the build
    read (a factor or a logical), each `$setResponse` coding back, about 70 us at 20,000 rows.
24. **One frame column for each predictor column the engine reads**: a matrix-valued term split into its
    columns, named as the block names them, so every setter writes one vector and none copies a matrix.
    `model.frame(s$data)` then shows `poly(x3, 2)1` and `poly(x3, 2)2` where `model.frame(o)` shows one matrix
    column. Alternative: matrix columns kept as `model.frame()` keeps them, a column's setter copying that
    matrix (0.45 ms for a 15-column matrix at 20,000 rows).
25. **Under `factors = "indicators"` the frame holds the indicator columns**, so a 0/1 value given to one is
    taken, as in 0.9-34 and today. Alternative: the factor held, the indicators made at creation, a write to
    one refused and the factor taken by name.
26. **A character predictor is held as a factor of its recorded levels**, coded once at the resolution.
    Alternative: held as characters and coded at each creation.
27. **The offsets are summed into one `(offset)` column**, as the sampler adds one; `$setOffset` then writes
    the total, which separate `offset()` columns could not take. Alternative: none that keeps
    `model.offset()` right after a setter.
28. **A forest's basis column is named `(basis.k)`** by the forest's position, as `bases` is indexed and as
    model.frame names `(weights)`. Alternative: by the forest's label, or by the basis's expression, where two
    forests naming one expression would share a column that `$setForestBasis` must split.
29. **The test rows stay in `x.test` and their kin**: they come from `test`, and no model frame holds them.
    Alternative: a test frame, about 40 lines more.
30. **The bridge is handed today's layout, made by R from the frame before each call**, sharing its columns,
    so the engine and the bridge's parsing are untouched. Alternative: the bridge parses the frame, about 150
    lines under `src/` and every parse site retested.
31. **The whole-matrix setter splits the matrix in the bridge**, about 15 lines, so it costs what it does
    today. Alternative: the split in R, 0.22 ms more per call at 20,000 by 15 (3 percent of the call, 7
    percent of an iteration).
32. **On a sampler's data the training row names and `na.action` are the frame's**, as `model.frame()` gives
    them; the slots keep the test row names only. Alternative: the slots kept beside the frame.
33. **A sparse column is held as a `sparseVector` or `sparseFactor` column**, a `dgCMatrix` split into one
    each, the container's `dgCMatrix` assembled at each bridge call. A column's setter then rewrites one
    column's entries, not the whole matrix's. Alternative: the `dgCMatrix` kept whole in a slot beside the
    frame.
34. **An `aft` response is a two-column (time, status) matrix**, as `model.frame()` keeps a `Surv` response;
    the `status` slot goes. Alternative: a `(status)` column.
35. **A multinomial response is its count matrix**; the trials are made from it where read. Alternative:
    the trials stored beside it, a second copy.

Also made here:

12. Rows stored as positions, nothing stored when every row is kept in order. Alternative: always stored.
13. The internal row-cutting paths act on a sampler's data and cut its frame. Alternative: cutting a user's
    object's row record in step, which nothing would read.
14. A forest's record names columns, never numbers. Alternative: numbers, which reordered newdata breaks.
15. A forest is chosen by label or position; the matrix interface takes values only.
16. `bases =` gets R's own "unused argument" error, the argument never having shipped.
17. A record of fewer than two forests stays refused at creation, and the `$setData` refusal is reworded
    (2.8). Alternative: a record of one plain forest read as none.
18. A user's object keeps `x.test` as built: it comes from `test`, not the frame.
19. A user's object keeps `na.action` and `rowNames`: the first is not given by the positions alone, the second
    is the names model.frame gave the rows. Alternative: replaying model.frame's naming.
20. `model.frame(o)` is made in the build's order: evaluate on the whole frame, cut, code. Alternative:
    predict's path on the kept rows, which fails on a caller's vector once a row is dropped.
21. A fit keeps its own basis vector wherever it keeps per-forest fits, its sampler kept or not: it cannot be
    made again (dec-B185), and a kept sampler's basis is its current one. This answers the question dec-B671
    sent here. Alternative: read from the kept sampler, wrong after `fit$fit$setForestBasis()`.
22. No reader for a data object saved by 0.9-34, which already fails on main.

### Costs dec-B678 did not weigh

- **The reads of a sampler's slots.** dec-B678 removes `x`, `y`, `weights` and `offset` from a sampler's data,
  so every internal read of them moves to an accessor (about 360 lines of R) and so do the tests that read
  them (up to about 740 lines). No draw moves and no setter slows; it is the bulk of the build.
- **stan4bart's saved bundle.** It keeps the specification's data with the block cut to no rows and splices
  the fit's `bartData@x` back at a restore. Under dec-B678 that data holds the model frame, and neither the cut
  nor the splice can be written as a slot write. Shown: (a) the bundle keeps the specification's data whole,
  4 lines removed from stan4bart, the bundle growing by the predictor columns (0.050 to 0.088 MB at 750 rows
  and 25 columns, about n times p times 8 bytes in general, so 2.2 MB at 20,000 by 14); (b) the bundle keeps
  no data and a restore makes the specification again from `bartData` by the call the fit made, about 10 to 15
  lines in stan4bart, the bundle as today less the response and weights; (c) dbarts offers a way to make a
  sampler's data from a block-holding object and a response, stan4bart calling it at a restore, about 15
  lines in dbarts and 5 in stan4bart. Not decided here; (b) is recommended, keeping a saved fit free of a
  second copy of its predictors.
- **bartCause's reads of the samplers' term labels**, `attr(fit$fit$data@x, "term.labels")` in two places,
  which a stand-in answers with `NULL` in silence. Taken: the two lines move to the model frame's terms.
  Alternative: the stand-in carrying the coding record's attributes, so `attr()` reads answer rightly, which
  changes what dec-B676 says `attr()` gives.
- **A matrix-interface object kept beside its sampler** holds the predictors twice in a session, the
  sampler's frame splitting the matrix (2.4 MB at 20,000 by 15). Taken, for call 24's faster setters.
  Alternative: the matrix kept as one matrix column on such a sampler, today's sharing and today's copy at each
  column's setter.
- **A saved sampler writes its data twice** on main, as its `data` and as its engine pointer's protected
  value, 4.55 MB against 0.9-34's 2.32 MB; under call 30 the pointer protects the made layout, so this stays.
  It predates this design; whether the bridge must keep the data object behind the pointer is an engine
  question, and two comments in the bridge disagree on whether anything borrows from it.
