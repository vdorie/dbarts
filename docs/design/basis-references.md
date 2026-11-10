# A forest's basis names a column of the data

Status: PROPOSED, 2026-11-07, revised 2026-11-08 to the rulings dec-B674 to dec-B688, dec-B692 and dec-B704 to dec-B707; part 1 LANDED 2026-11-10
([basis-references-part-1.md](../plans/basis-references-part-1.md)), part 2 not built.
Follows dec-B663 (a forest's basis is a reference to a column of the data), dec-B666 (the frame kept on the
data object), dec-B670 (the whole frame kept, and beside it as little derived from it as can be), dec-B671
(what a fit of several forests keeps under `keepFits = FALSE`), dec-B674 (`keep.data = FALSE`, the frame kept
by default), dec-B675 (the response, weights and offsets made from the frame), dec-B676 (`o@x` a stand-in
that stops when read), dec-B677 (`model.frame(o)` and `nobs(o)` the accessors; `extract()` keeps to samples),
dec-B678 (a sampler's data holds the model frame and what its forests need, the setters writing into it),
dec-B679 (two landings), dec-B680 (a basis vector has one entry per row handed over), dec-B681 (`dbarts(o)`
reads a declared basis as `dbartsSpec(o)` does), dec-B682 (beside an object with no caller frame, a name it
holds is read from it; it revises dec-B674's words), dec-B683 (an environment's object holds what it read at
the build), dec-B684 (a forests list states the forests in full), dec-B685 (`model.frame(s$data)` decodes
where exact), dec-B686 (the `s$data@x` stand-in keeps the coding record's attributes), dec-B687 (stan4bart's
bundle rebuilds from `fit$bartData`), dec-B688 (the row mask is an `(active)` column) and dec-B692
(`extract()` on a sampler returns only draws).

## Two landings

The design lands in two parts (dec-B679):

- **Part 1, the data object and the forest doors**: sections 2.1 to 2.3, 2.7 to 2.9, the user's object's
  stand-ins, `model.frame()` and `nobs()`, `keep.data`, and the rulings on calls 1, 2, 3, 6 and 7. A sampler's
  data is today's object with real slots: the caller's frame dropped where a sampler is handed its data, its
  column names kept. About 495 lines of R added, above dec-B679's 305 by about three fifths (the plan says
  why).
- **Part 2, the sampler's model frame** (dec-B678): sections 2.4 to 2.6, the setters, the stand-ins on a
  sampler's data, the lint rule and the bridge entry. Calls 39, 36 and 41 and stan4bart's bundle are ruled
  (dec-B685, dec-B686, dec-B688, dec-B687), and `extract(s, "predictors")` goes with it (dec-B692).

Between the two, main holds a user's object of the new shape beside a sampler's data of today's.

## In short

One class, one set of slots, two kinds of data object:

- **A user's object**, `dbartsData(formula, data)`, keeps the caller's frame whole, the positions of the rows
  it kept and the record of how its variables are coded. It stores no predictors, response, weights or
  offsets: `o@x`, `o@y`, `o@weights` and `o@offset` are stand-ins that stop, naming `model.frame(o)`.
  `model.frame(o)` makes the model's variables from the frame; `nobs(o)` counts its rows.
  `keep.data = FALSE` keeps no frame and keeps the block, response, weights and offsets in their slots, as
  today; bartCause and stan4bart pass it.
- **A sampler's data**, `s$data`, made where a sampler is handed its data, holds the model frame as the
  sampler runs on it: one column for each predictor column the engine reads, a factor kept a factor, the
  response as the sampler reads it, the weights, the summed offset, the active-row mask, and each forest's
  basis and weights. The caller's frame is dropped. The R setters write into it. `model.frame(s$data)`
  returns it decoded to the model's own shape, so that it matches `model.frame(o)` but for the differences of
  call 39. `s$data@x` and its kin are stand-ins carrying the coding record's attributes.

The engine's input is unchanged: before each call that hands the bridge a data object, R makes today's
layout from the frame, sharing its columns. A forest's basis given in a forests list beside a user's object
is read against its frame as `dbarts(f, d, forests = )` reads it; `dbartsData(bases = )` goes.

These figures are from today's main, a shipped build, at 20,000 rows, 14 numeric columns and a factor of 26
levels, 75 trees, one thread. Each is the median of 7 timed loops after 3 warm-up loops, with the range in
brackets:

- **Setters**: one iteration takes 2.94 ms (2.90 to 2.98). `$setResponse` is 120 us (120 to 130) and gains
  about 7 us, the cost of a frame column's write against a slot's. A one-column `$setPredictor` is 2.08 ms
  (2.08 to 2.10) and does not change. On the matrix interface a column's setter takes 3.84 ms (2.98 to 3.86)
  today and drops to about 2.1 ms, because it no longer copies the whole matrix. A whole-matrix
  `$setPredictor` keeps today's 6.3 ms (6.2 to 6.5), its split made in the bridge (call 31).
- **Saved sizes** stay within 0.001 MB of today's: a sampler is 4.55 MB, a `bart` fit kept with its sampler
  11.82 MB, bartCause's fits 2.18 and 2.34 MB, and stan4bart's fit 0.46 MB if its bundle keeps what it keeps
  today (Costs dec-B678 did not weigh). A saved sampler on main writes its data twice; 0.9-34 wrote it once,
  2.32 MB (call 30).
- **Cost**: about 600 new lines of R, 360 lines changed where R reads a sampler's slots, a custom lint rule,
  about 15 lines in the bridge, up to about 740 test lines moved and 300 new, seven help pages, one NEWS
  entry, 6 lines in bartCause, and in stan4bart 4 lines plus its bundle (open).

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
- **A data object** stores the predictors in `x`: a columnar container on the formula path, with a factor
  as one column, or a plain matrix or a `dgCMatrix` on the matrix interface. It also stores `y`,
  `weights`, `offset`, `x.test` and their kin, and no frame. In a session its columns are the caller's
  frame's columns, shared; on the matrix interface `x` is a copy of the caller's matrix.
- **A sampler's data** is that object, rewritten by the R setters:
  - `$setResponse` writes `y`.
  - `$setPredictor` of a column replaces that column of the container; on a plain matrix it copies the whole
    matrix first.
  - A whole-matrix `$setPredictor` replaces `x` with a plain matrix.
  - `$setOffset`, `$setWeights`, `$setTestPredictor` and `$setForestBasis` write their slots.
  - Two per-row quantities live outside it, in fields of the sampler: `$setForestWeights`'s per-forest weights
    and `$setActiveRows`'s mask. Getting the engine back, `$setState` and `$copy()` each apply them again.
  - A probit or ordinal fit's 0/1 weights ride `s$data@weights` and survive `$setData(s$data)` (run on main
    at 7a54b2b1d); only a mask set by `$setActiveRows` is cleared. The two channels compose.
  - The flat C API's `dbarts_sampler_setResponse` and `dbarts_sampler_setOffset` change the engine and write no
    slot. stan4bart calls the second every iteration.
  - The engine keeps no predictor block. It reads the data object at creation, at `$setData`, at a model
    change, when a state is put back and at a warm start (`$installTrees`).
- **0.9-34** had no forests, no bases and no `dbartsSpec()`. Its data object held the model matrix in `x`,
  with factors spread into indicator columns; `attr(s$data@x, "term.labels")` gave the term labels (run).

## 2. The design

### 2.1 Rows: what is numbered and how

A user's object records the number of rows handed over and the positions of the rows kept. It is the index
`dbarts()` already builds to cut a basis written beside a frame, so a value and a code basis are cut by one
index. With 12 rows, rows 3 and 9 missing a predictor:

- no `subset`, default na.action: 1 to 12 (`na.keepPredictors` keeps rows with a missing predictor)
- `subset = 2:10`: 2 to 10; `subset = -1`: 2 to 12; `subset = z == 1`: 5, 9, 10, 11, 12
- `subset = c(10, 2, 5, 7)`: in that order; `c(2, 2, 5, 7, 8)` keeps the repeat
- a character `subset` by row name: the positions of those names
- `na.action = na.omit`: rows 3 and 9 go; after `subset = 2:10` too: 2, 4, 5, 6, 7, 8, 10

An integer vector, 4 bytes a row; empty when every row is kept in order. A sampler's data needs none: its
frame holds the rows in, in order. The internal paths that cut or expand rows (the hurdle's two parts and
its merge, the hazard's person-period expansion, `xbart()`'s folds) act on a sampler's data and cut its
frame. Each cut clears the frame's `na.action` attribute, which `[` keeps (run) and which would otherwise
number rows of another frame for `naresid()` (call 43).

### 2.2 Columns: what a user's object keeps

The whole frame, as given, when `data` is a data frame or a list (dec-B670). In a session it costs nothing:
R shares the columns, and a later `d$z[1] <- 9` copies `z` alone, the object keeping the old one. Saved, at
20,000 rows and 14 columns, the frame is 1.19 MB. A column added to `d` later is not in the object. A name
that is no column is looked up where `forest()` was called. Set aside: the frame less the formula's columns
(a basis naming a predictor then falls through to the caller), numbers into the caller's frame (a reorder or
a dropped row misreads them in silence) and named columns only (dec-B666).

`keep.data = FALSE` keeps no frame (dec-B674). The object keeps the block, `y`, `weights` and `offset` as
today, the positions of its rows (dec-B680 cuts a full-length basis by them) and the frame's column names.
Beside it a basis naming a column the object holds is read from it, and another of the frame's names is
refused, naming the object built with its frame (dec-B682, revising dec-B674's words).

### 2.3 A user's object, under dec-B670 and dec-B674 to dec-B677

Kept, beside the frame and the rows:

- **The coding record**. This is the block with no rows: its column names, `varTypes`, term labels, factor
  levels (or indicator levels and the dropped columns) and frozen terms (dec-B185: `poly()` coefficients,
  `scale()` centers, a caller's value), with a flag for each column saying whether it held a missing value.
  It is 750 bytes where the formula names only columns; a caller's vector named in the formula rides at full
  length, 8 bytes a row.
- **The response's coding**: `response.type`, `response.n.levels`, `response.levels`, `counts.as.factor`.
- **The frozen expressions** of the response, `weights` and `offset`. A column stays a name, and a caller's
  value is held as its value (dec-B675).
- **`na.action`**: the rows the na.action dropped, numbered within the subset, which the positions alone do
  not give.
- **Row names only where the frame cannot give them** (call 19). Under a repeated subset, model.frame names
  the kept rows "2", "2.1", "5"; `make.unique()` of the frame's names at the positions gives them back (run).
  It does not when a repeated row is one the na.action then drops, and only then are the names stored. Today
  the object holds a character vector of every kept row's name: 938 KB at 14,999 rows.
- **From the caller's other arguments**, kept as built: `x.test`, `weights.test`, `offset.test`, `n.cuts`,
  `sigma`, `offset.category` and its test twin, `useQuantiles`, `testUsesRegularOffset`, `missing`.

Not stored: the predictors, the response, the weights and the offsets. `o@x`, `o@y`, `o@weights` and
`o@offset` each hold a small classed stand-in. Its `dim()`, `nrow()`, `colnames()`, `[` and `as.matrix()`
stop with a message naming `model.frame(o)`. `attr(o@x, ...)` gives `NULL` without a word, as dec-B676
rules; the help says so (call 36 asks whether to change this).

`model.frame(o)` returns the model's variables over the rows in, after `subset` and the na.action, as
`model.frame()` of an `lm` fit does. Factors stay factors, the weights sit in a `(weights)` column and the
offsets in `offset()` and `(offset)` columns, so `model.response()`, `model.weights()` and `model.offset()`
read them.

`nobs(o)` counts every row in. Unlike `nobs()` of an `lm` fit it keeps the rows of weight 0, which a BART fit
holds in its trees; the help says so (dec-B677).

`model.frame(o)` is made in the build's order: the frozen terms are evaluated on the whole frame, the
variables cut to the kept positions, then checked against the recorded levels. Run, the coded predictor
columns this makes equal the built block in every case tried:

- plain columns, `y ~ .`, logical and matrix columns
- `offset()` terms, `I(x2^2)`, and `log(x2) + scale(k)` and `poly()` under `subset`, once with a repeated row
- factors and ordered factors under a `subset` that empties a level
- character columns in either coding
- `na.omit`, `na.exclude`, `na.pass`, and the default `na.keepPredictors` with a missing factor value
- a value, or a missing value, that only a dropped row holds
- a caller's scalar, and a caller's vector under `subset` and under `na.omit`
- a sparse factor

Replaying `predict()`'s path on the kept rows instead fails on a caller's vector once a row is dropped
(call 20). It takes 0.0014 s at 20,000 rows and 15 columns under a subset of 15,000, against 0.004 s for the
build.

`makeTestModelMatrix(data, newdata)` (help `makeind`) reads the coding record and its missing-value flags. On
a default object it needs no block and gives what it gives today.

Saved, a user's object is the frame and the record: 2.22 MB at 20,000 rows with a 26-level factor, where frame
and block together would be 4.5 MB.

**An object built on an environment** (dec-B683, call 7) holds the model frame made at the build, the shape a
sampler's data holds (2.4; in part 1, the block, as a sampler's data holds it then), and the environment's
variable names. A basis naming a model variable is read from it; another of the environment's names is
refused, naming a data frame or list. Every value is read at the build.

### 2.4 The sampler's data, under dec-B678

**Where it is made.** The sampler's data is made in the resolution, the step these doors share: `dbarts(o)`,
`dbartsSpec(o)`, `bart(o)`, `rbart_vi(o)`, `xbart(o)`, the hurdle's entry, `new("dbartsSampler", data = o)`
and `$setData(o)`. It is the object `dbartsSpec()` returns as its `data`, which stan4bart keeps in its saved
bundle. A sampler's data handed in again is taken as it is.

**What it holds.** Each row below is a column of one data frame in the `model.frame` slot. The frame carries
terms built from its own columns, the response at position 1, its term labels its predictor columns and no
`offset` attribute (call 38):

| column | holds |
|---|---|
| the response | as the sampler reads it: coded doubles (0/1, ordinal codes); for `aft` the time, its status in `(status)` (call 34); for `multinomial` the count matrix (call 35) |
| each predictor column the engine reads | a numeric column as doubles; a factor as a factor of the recorded levels, a character column coded to one (call 26); a matrix-valued term (`poly()`, `scale()`, a matrix variable) as one column each (call 24); under `factors = "indicators"`, each indicator column (call 25) |
| `(weights)` | the weights, absent where there are none; a probit or ordinal fit's 0/1 weights among them, which ride the data today |
| `(active)` | the mask `$setActiveRows` installs, in every family, composing with the weights; absent where every row is in (call 41) |
| `(offset)` | the one offset the sampler adds, the formula's `offset()` terms summed with the argument (call 27) |
| `(offset.category)` | the multinomial category offset |
| `(basis.<label>)` | a forest's basis, a number or a two-level factor, for each forest that has one (call 28) |
| `(forest.weights.<label>)` | a forest's per-row weights, for each forest given some (call 40) |

The frame's row names are the names of the kept rows and its `na.action` attribute records what the
na.action dropped, both as `model.frame()` gives them, so those two slots go (call 32).

Sparse predictor columns are not columns of this frame. A sparse-bearing design keeps them as one
`dgCMatrix` in the slot `x.sparse`, as the container keeps them today (call 33). No dense column is in it,
so it is no second copy. A data frame of `sparseVector` columns prints as corrupt, and `model.frame.default`
refuses it ("invalid type (S4)").

Kept as records, none a copy of a row:

- the coding record (with the formula's terms) and the response's coding
- `basis.levels` and `basis.scale`
- `n.cuts`, `sigma`, `useQuantiles`, `missing`
- what a sampler writes: `missing.seen`, `cutPoints`, `response.range`, `forest.calibration` and
  `standardization`

The caller's frame and the row positions are dropped. The frame's column names are kept for dec-B682.

The test rows stay where they are, in `x.test`, `weights.test`, `offset.test` and `offset.category.test`
(call 29). They come from `test`.

**Slots that go on a sampler's data**: `bases`, `counts`, `status`, `offset.category`, `na.action` and the
training part of `rowNames`. `x`, `y`, `weights` and `offset` hold stand-ins naming
`model.frame(s$data)`. The stand-in in `x` carries the coding record's attributes (`term.labels`,
`factor.levels`, `varTypes`, `terms`, `drop`) but no `dim`, so `nrow()` and `[` still stop and `attr()` reads
answer rightly (call 36). Two sampler fields go as well, `forestWeights` and `activeRows`.

**What `model.frame(s$data)` returns.** The frame decoded on read to the model's own shape (call 39):

- the response in the type and levels the build read: a factor, a logical, a `Surv` response rebuilt from
  the time and `(status)`
- a character predictor back to characters
- a split matrix term joined back into its matrix column, with its attributes (`poly()`'s coefficients)
- the formula's terms, with its `offset()` terms removed and their sum in `(offset)`
- `(weights)` as the build read them

Run, the decoded frame gives `model.offset()`, `model.response()`, `model.weights()` and `model.matrix()`
equal to those of `model.frame(o)`. The stored frame gives `model.offset()`, `model.response()` and
`model.weights()` equal to what the sampler runs on, and `model.matrix()` runs on it. Under the formula's own
terms a stored frame would add a `poly()` column into the offset with no message. Decoding costs time on
read only, 0.07 ms for a factor response at 20,000 rows, and no setter pays it.

**The engine's input.** Each call that hands the bridge a data object is handed today's layout, which R makes
from the frame just before. That covers creation, `$setData`, a model change, `xbart()`'s data handle, and
`$setState`'s and `$installTrees`' raw predictors.

- The predictor columns go in as the columnar container. A double or factor column is shared, not copied; an
  integer or logical one is turned to doubles.
- `x.sparse` is passed as today's sparse part.
- The record's attributes are set.
- The response, weights, mask, offset and bases go in as the values the bridge reads.

Run: 0.0035 ms (0.0035 to 0.0045) to make the container at 20,000 rows and 15 columns. The cut points and
codes are what the engine derives from it, as today (call 30).

**The R setters and doors**, what each writes and what it costs per call. The times are today's main at
20,000 rows, as in In short, where one iteration takes 2.94 ms. A frame column's write costs 31.5 us
(31.5 to 32.0), against 24.5 us (24.0 to 24.5) for today's slot write:

| setter or door | writes into the frame | today | under dec-B678 |
|---|---|---|---|
| `$setResponse(y)` | the response column | 120 us (120 to 130) | about 127 us |
| `$setResponse(y, status)`, aft | the time and `(status)` | 280 us (270 to 320) | about 287 us; a two-column matrix rebuilt at each call would add 35 us |
| `$setPredictor(v, j)` | column j | 2.08 ms (2.08 to 2.10) | the same: one column replaced, as today's container replaces one |
| the same, a factor column | the factor, matched to the recorded levels | as today | as today; a factor with the recorded levels in order is stored as given |
| the same, on the matrix interface | column j | 3.84 ms (2.98 to 3.86) | about 2.1 ms: the whole matrix is no longer copied |
| `$setPredictor(x)`, whole | every predictor column | 6.3 ms (6.2 to 6.5) | about 6.3 ms, the bridge splitting the matrix in one copy (0.08 ms) in place of today's `as.double` copy (call 31) |
| `$setPredictor(x, forceUpdate = "partial")`, `updatePredictorPerObservationJointly` | the installed rows of the column | as today | as today |
| `$setOffset(o)` | `(offset)`, removed for `NULL` | 140 us (140 to 145) | about 147 us |
| `$setWeights(w)` | `(weights)`, removed for `NULL` | 160 us (160 to 170) | about 167 us |
| `$setActiveRows(a)` | `(active)` | 240 us (240 to 260) | about 247 us |
| `$setForestWeights` | `(forest.weights.<label>)` | | one column write |
| `$setForestBasis` | `(basis.<label>)` | | one column write |
| `$setCounts`, multinomial | the response's count matrix | | one column write; the trials are made from it where read |
| `$setTestPredictor`, `$setTestOffset` | `x.test`, `offset.test` | | unchanged |
| `$setData(o)` | a new frame, made at the resolution | 7.7 ms | about 7.7 ms from a block-holding object; 1.4 ms more from a user's object, to make the frame |
| `$setState` | nothing; the raw predictors are read from the frame | 4.05 ms on a sparse design (p = 15), 3.85 ms (p = 500) | the same: the container is made in 0.0035 ms and `x.sparse` is passed as it is |
| `$installTrees` | nothing; the raw predictors and the donor's levels are read from the frame and the record | | unchanged |
| `$copy()` | the copy shares the frame in a session | 12 ms | unchanged |
| save and reload | the frame is written; a reload creates the engine from it | | the same size (2.6) |
| `xbart()`'s folds | the folds cut the container made once at its start | | unchanged |

**Samplers driven through the flat C API.** `dbarts_sampler_setResponse` and `dbarts_sampler_setOffset`
change the engine and write no frame, as they write no slot today (call 42). Only the R setters write the
frame, so `model.frame(s$data)` shows the last values given through R. A copy or a reload rebuilds from
those.

**Reads of the old slots.** The R code reads a sampler's slots in about 360 lines, 16 of them `attr()` reads
of `data@x`: 8 of `term.labels`, 5 of `terms` (one in `partialDependence.R` chooses its path by it) and 3 of
`factor.levels` (in `data.R`, a missing one falls back to counting levels from the codes). All of them move
to internal accessors that read the frame and the record. A custom lint rule, run with the others, refuses
`@x`, `@y`, `@weights` and `@offset` outside the accessors' file, so no read is left behind (call 37). The
stand-in's attributes (call 36) keep a read through a local variable, which the lint rule cannot see,
answering rightly.

**A factor column given to `$setPredictor`** takes a factor or a character vector of the column's labels,
matched by label, as today. A level the column does not declare is refused by name, and so are numbers and
codes, as today: a code is no label. A factor whose levels are the recorded ones in order is stored as
given; another is recoded to them. Under `factors = "indicators"` the frame holds the indicator columns, and
a 0/1 vector given to one of them is taken, as 0.9-34 and today take it; both refuse the factor by name
(run, call 25).

**The matrix interface.** `dbarts(x, y)` makes the frame of the matrix's columns, one column each under the
matrix's column names, the response first, and the weights and offset as above. An integer matrix is turned
to doubles, as today. A `dgCMatrix` is kept as `x.sparse`.

**An unnamed matrix** (dec-B653). The frame's columns need names, so they take placeholder names, and the
coding record says the predictors are unnamed:

- `colnames()` of the predictors stays `NULL`.
- An unnamed matrix given later is read in column order with no message.
- A named one given to an unnamed sampler still warns "'x' and 'test' are not both named".
- `model.frame(s$data)` shows the placeholder names.

`y ~ X` on an unnamed `X` keeps the names model.frame gives, `X.1` and `X.2`, as today.

**The state.** A sampler's state holds the trees, the variances and the chains' scalars, about 0.016 MB
here, and no data (state never carries model). It is unchanged. A saved sampler is its control, its model,
its data with the frame, and its state.

### 2.5 What a block-holding object and its sampler hold, and what the consumers read

Under dec-B674 an object built with `keep.data = FALSE`, and the matrix interface's object, hold the block,
`y`, `weights` and `offset` in their slots, as today, and the frame's column names. A sampler made from one
holds the model frame made from those slots: a container's columns are shared in a session, and a matrix's
are split into columns. A user who keeps a matrix-interface object beside its sampler holds the predictors
twice in a session, n times p times 8 bytes (2.4 MB at 20,000 by 15); `dbarts(x, y)` with no object kept does
not.

What bartCause and stan4bart read, set against the examination of what each should do:

- **bartCause** passes `keep.data = FALSE` at its two data call sites (three objects). Its saved fits keep
  `data.rsp` as today, and its reads of `nrow()`, `colnames()`, `x[, z]`, `as.matrix(x)`, `@y` and `@weights`
  and its write of `@x.test` stay.
  - Its `bcf()` gives the treatment forest `basis =` with the full-length `factor(z)` it already builds. This
    needs calls 1 and 2 (run: draws identical to the data door, with and without `subset`).
  - It reads `attr(fit$fit$data@x, "term.labels")` off the samplers its fits keep, in two places, and under
    call 36 that read still answers rightly.
  - About 6 lines.
- **stan4bart** passes `keep.data = FALSE` in `glFormula`, so its fit's `bartData` reads (`@x`, `@varTypes`,
  `factor.levels`, `makeTestModelMatrix()`) stay.
  - New under dec-B678: it reads `sampler$data@y` and `sampler$data@weights` when it installs a response
    range, and `data@offset` on the specification's data at a restore. These 3 lines move to
    `model.response()`, `model.weights()` and `model.offset()` of the model frame.
  - Its saved bundle strips the block off the specification's data (`@x` and `@x.test` cut to no rows) and
    splices `bartData@x` back at a restore. Once the data holds a model frame, neither can be written as a
    slot write. The ways out are in Costs dec-B678 did not weigh.
  - Its per-iteration `dbarts_sampler_setOffset` is call 42.
- The examination's other recommendations are untouched by this design: bartCause's test frame in place of a
  hand-built `x.test`, its `plot_support` on a fit with missing responses, and stan4bart's hand-written slots
  moved to arguments.

### 2.6 Saved sizes

Run on a shipped build of main. The dec-B678 column takes each sampler's data with the frame of 2.4 in place
of `x`, `y`, `weights`, `offset` and `bases`, sized as the frame saves from the same values.

| object | today | under dec-B678 |
|---|---|---|
| a `dbarts()` sampler, formula, 20,000 rows, 14 numeric and a 26-level factor | 4.55 MB (data 2.23, the engine pointer 2.22) | 4.55 MB (data: 2.220 to 2.219) |
| a `bart` fit with `keepTrees = TRUE`, its sampler kept, 50 draws | 11.82 MB | 11.82 MB |
| the same with `keepTrees = FALSE` | 7.46 MB | 7.46 MB |
| a matrix-interface sampler, 20,000 by 15 | 4.81 MB (data 2.36) | 4.81 MB (the frame of columns 2.34) |
| stan4bart, 750 rows, 25 BART columns, one chain | 0.46 MB; bundle 0.050, its data 0.007 | 0.46 MB if the bundle keeps no predictors; 0.50 MB, bundle 0.088, if it keeps the specification's data whole |
| bartCause `bartc()`, 750 rows, 25 covariates, bart response | 2.18 MB | 2.18 MB (each kept sampler's data within 0.001 MB) |
| the same, bcf response | 2.34 MB | 2.34 MB |

0.9-34 saved the same sampler at 2.32 MB, writing its data once.

### 2.7 The doors

- **`dbartsData()`** keeps the frame, the row record and the coding record (2.3); `keep.data = FALSE` keeps
  the block. `bases =` goes. Its formula still refuses a `forest()` term, pointing to a forests list.
- **`dbartsSpec(o, forests = )`** and **`dbarts(o, forests = )`** read a basis against `o`'s frame, then where
  `forest()` was called, built on `o`'s rows by the function `dbarts()` uses. So
  `dbartsSpec(dbartsData(f, d), forests = L)` and `dbarts(f, d, forests = L)` draw the same.
  - A value has one entry per row handed over, before `subset` and the na.action, and is cut with them; a
    vector of the cut length on an object that dropped rows is refused (dec-B680, call 1).
  - `dbarts(o)` reads a declaring list exactly as `dbartsSpec(o)` does (dec-B681, call 2).
- **Beside an object with no caller frame** (a sampler's data, a `keep.data = FALSE` object), a basis naming
  a column the object holds, a model variable or a held basis, is read from it; another of the caller's
  frame's names is refused, naming the data object built with its frame; any other name is the caller's
  value (dec-B682, call 6). A held basis is named by the expression it was written as, matched deparsed and
  read from the held columns; a source column of it is refused naming that expression (dec-B706).
- **A basis name is read from the data first, then where `forest()` was called, at every door**, a matrix's
  columns included: `dbarts(x, y)`, `dbartsData(x, y)` and its sampler's data read `x`'s column; a formula
  with no data reads where it was written (dec-B705). A single string names its column (dec-B700).
- **Beside any data, a sampler's data included, a forests list states the forests in full** (dec-B684, call
  3): a forest declaring no basis has none, keeping a held basis means naming it, a list of tree counts alone
  gives plain forests and no list gives one forest. dec-B663's interim refusal is not built. A forest
  naming a held basis keeps its recorded scale, matched by name; a new value takes its own (dec-B704).
- **A hazard fit's basis** has one value per subject, cut by `subset` and the na.action, repeated over its
  periods as weights and offsets are; a period-length vector is refused (dec-B707).
- **`$setForestBasis`** and **`predict`**: as built.
- **The matrix interface**: a value, full length, cut by the record, or a name of one of its columns
  (dec-B705).
- **Fit and reload**: `dbartsSpec()` and `dbarts(o)` write the model record `dbarts(f, d)` writes.

### 2.8 Data objects that carry bases, and setData

Once `bases =` goes, a data object carries basis columns only where a sampler of several forests made them:
`s$data` and what is made from it. A sampler of one forest makes none.

- **`bases = list(NULL)`**, one forest with no basis, has no path: a record of fewer than two forests stays
  refused at creation (call 17).
- **`$setData` on a sampler of one forest**, given another sampler's data that carry basis columns, is
  refused before anything is stored. The text says they are the data of a sampler of several forests and
  points to the data object they were made from or to `dbarts()`. Its own data handed back is taken, and now
  keeps its mask and forest weights, which live in the frame.
- **`$setData` on a sampler of several forests**: refused whatever it is given, as today.
- **`dbarts(s$data)` and `dbartsSpec(s$data)`**: the list states the forests (dec-B684), a held basis kept
  where a forest names it. **A reload** re-creates the forests from the model record and the data's basis
  columns.

Three texts that advise the data door are reworded: the `$setData` refusal, `dbarts()`'s refusal of a
declaring list, and the refusal of a model of one forest with a basis.

### 2.9 What else this needs

- **The slots**:
  - `frame`: the caller's data, on a user's object only.
  - `rows`: its positions, `NULL` when every row is kept in order.
  - `model.frame`: on a sampler's data and an environment's object.
  - `x.sparse`: a sparse-bearing sampler's sparse columns.
  - `frame.columns`: the frame's or the environment's names.
  - the coding record.

  Where `frame` or `model.frame` is set, `x`, `y`, `weights` and `offset` hold stand-ins; where neither is,
  they hold the block and the values. The class's validity checks each kind. The help describes the slots as
  read-only (dec-B554) and shows no slot write.
- **The record's check**: the positions lie within the rows handed over, and a model frame's row count, its
  response's length and each per-row column's length agree. Sampler creation and the validity call it.
- **The lint rule** (call 37) in the package's lint configuration, with the accessors' file exempt.
- **`dbartsData(o)`** handed a data object returns it as it is.
- **A data object saved by 0.9-34** already fails on main ("no slot of name "bases""), and nothing changes
  there (call 22).
- **A comment to fix**: the one above the row numbers `dbarts()` hands to the data object says a code basis is
  built only on the rows kept, but it is read over every row and then cut.

## 3. The landscape

Run where installed (mgcv 1.9.4, nlme 3.1.171, lme4 2.1.0, brms 2.23.0); stochtree and bcf read from source.

- `lm`, `glm`: a varying slope is a term, and names are found in `data`, then the formula's environment. The
  fit keeps its model frame (`model = TRUE`) and `model.frame(fit)` returns it, the weights in `(weights)` and
  the offset argument in `(offset)`. `nobs()` drops rows of weight 0.
- mgcv `s(x2, by = z)` and lme4 `(1 + z | g)`: the pre-built object (`gam(fit = FALSE)`, `lFormula()`) holds
  the model's columns, and `predict` needs them.
- nlme: the fit keeps the whole data frame unless `keep.data = FALSE`.
- brms: by-variables must be in `data` or `data2`.
- stochtree, bcf (Hahn, Murray and Carvalho): values only. BART and bartMachine: no basis.

## 4. What a user meets

- A formula user, a forests-list user beside a frame, a matrix-interface user of `dbarts()`: no draw moves.
- A reader of `s$data@x`, `s$data@y`, `s$data@weights` or `s$data@offset` meets an error naming
  `model.frame(s$data)`, which holds them. `attr(s$data@x, "term.labels")` still answers.
  - In 0.9-34, `s$data@x` held the model matrix and `s$data@y` the response. Code that runs the BART sampler
    inside a larger sampler and reads them meets this.
- A reader of `o@x` and its kin on an object built from a frame meets the same error, naming `model.frame(o)`.
- A caller of `$setPredictor` on the matrix interface: a column's setter is faster.
- `$setData(s$data)` keeps the sampler's mask and forest weights, where today it drops the mask.
- treatSens uses the matrix interface and reads no slot, so nothing moves.
- bartCause and stan4bart: 2.5.
- A `dbartsSpec()` or `dbarts(o)` caller with a basis (part 1):
  - A name that is a column now finds it. Where it is also a caller's variable, the column wins, as at
    `dbarts()`, and the draws move.
  - Beside a sampler's data or a `keep.data = FALSE` object, a name the object holds is read from it where
    today the caller's variable is (dec-B682), and the draws move.
  - Beside a sampler's data of several forests, no list gives one forest and a list of tree counts plain
    forests, where today the held bases are kept (dec-B684), and the draws move.
  - A kept-length value on an object that dropped rows is refused; `dbarts(o)` takes a declaring list.
- Anyone who saves a data object built from a frame: it holds the frame in place of the block, about the same
  size where the frame is mostly the model's columns, larger where it holds others.

## 5. The cost

- R, about 600 new lines and 360 changed, with nothing in `inst/include/dbarts/dbarts.h`:
  - the slots, the stand-ins on both kinds with the attributes on a sampler's data, and the validity, about 45
  - `dbartsData()` keeping the frame, the row and coding records, the narrowed row names, `keep.data`,
    about 50
  - `bases =` leaving the exported function, about 20
  - the resolution making the sampler's frame: from a user's object in the build's order, or from a
    block-holding object or an environment's, with its own terms, the response coded, the offsets summed,
    matrix terms split, the basis, mask and forest-weight columns, and `x.sparse`; about 130
  - `model.frame()` and `nobs()` methods, decoding a sampler's frame on read, about 60
  - the engine view made from the frame before each bridge call, about 40
  - the internal accessors, about 40, and the reads moved to them: about 360 lines that read or write `@x`,
    `@y`, `@weights` or `@offset` today (182, 73, 64 and 41)
  - the setters writing the frame, about 85: `$setResponse`, the three `$setPredictor` paths and the joint
    update, `$setOffset`, `$setWeights`, `$setActiveRows`, `$setForestWeights`, `$setCounts`,
    `$setForestBasis`, the hazard expansion, the hurdle parts, the `na.action` cleared at each cut. The two
    sampler fields go.
  - `dbartsSpec()` and `dbarts(o)` resolving against the frame, calls 1, 2, 3, 6 and 7, about 85
  - `makeTestModelMatrix()` reading the record, about 5
- The lint rule: about 25 lines in the package's lint configuration.
- `src/`: about 15 lines in the bridge, an entry that splits a double matrix into column vectors for the
  whole-matrix setter and `$setData` (call 31). No engine change, and no draw moves.
- Tests:
  - The 83 calls of `dbartsData(bases = )` in 24 files become forests lists that draw what they drew: 45
    fixtures, 36 pins of the forest door and 2 oracles on the list's value.
  - Reads of the old slots, about 740 lines in up to 90 files (472 `@x`, 135 `@y`, 59 `@weights`, 74
    `@offset`, counted by name, an upper bound), move to `model.frame()`, the internal accessors or
    `keep.data = FALSE`.
  - New tests, about 300 lines:
    - the sibling identity and the row positions
    - the frame's contents per family and interface
    - `model.offset()`, `model.response()` and `model.weights()` of `model.frame(s$data)` and of the stored
      frame, each equal to what the sampler runs on, with `offset()` terms, a `poly()` term and a factor
      response
    - the decoded frame against `model.frame(o)`
    - each setter writing its column, with draws identical to today's
    - `$setData(s$data)` keeping the mask and forest weights
    - the stand-ins' errors and attributes
    - a factor column's setter, indicator coding, an unnamed and a sparse matrix
    - the `na.action` at a cut
    - the saved size
    - `makeTestModelMatrix()` on a default object
    - the refusals of calls 3, 6 and 7
- Help: `dbartsData` (the slots, `keep.data`), `dbartsSampler-class` (the data field, the mask and forest
  weights), `dbartsSpec`, `dbarts`, `makeind`, a page for the `model.frame()` and `nobs()` methods, and the
  forests help where it names the data door. Four design notes that name the data door (forest-kind-by-class,
  public-surface, written-surface, multiplier-combiner) change at landing.
- NEWS, against 0.9-34, one entry: "A data object made by dbartsData() from a formula and a data frame keeps
  the frame and the positions of the rows it keeps, and a sampler's data holds the model frame, the setters
  writing into it. In both, the x, y, weights and offset slots, which held the model matrix and the response
  in 0.9-34, no longer hold them: model.frame(object) returns them and nobs(object) counts the rows. A forest
  given to dbarts() or dbartsSpec() beside a data object finds its basis in the frame. keep.data = FALSE keeps
  the matrix and the response in the object instead, as before."
- bartCause: about 6 lines and its suite. stan4bart: 4 lines, its bundle (open), about 5 test lines and its
  suite.

## Calls made

Ruled:

- call 4, by dec-B674
- call 8, by dec-B675
- call 9, by dec-B676
- call 10, by dec-B677: `model.frame(o)` and `nobs(o)`, and no `extract()` method for a data object
- call 5, by dec-B678: the caller's frame is dropped where a sampler is handed its data, and the sampler's
  data holds the model frame and what its forests need
- the former call 11, by dec-B671
- call 1, by dec-B680: a basis vector has one entry per row handed over, before `subset` and the na.action, at
  every door that builds a model; a vector of the cut length on an object that dropped rows is refused, the
  message saying the length expected
- call 2, by dec-B681: `dbarts(o, forests = )` reads a declared basis exactly as `dbartsSpec(o, forests = )`
  does; a data object holds the frame and no forest settings
- call 6, by dec-B682, which takes the alternative and revises dec-B674's words: beside a sampler's data or a
  `keep.data = FALSE` object, a basis naming a column the object holds (a model variable or a basis column) is
  read from it; a name that was another column of the caller's frame is refused, naming the data object built
  with its frame; any other name is the caller's value
- call 7, by dec-B683: an environment's object holds the model frame made at the build and the environment's
  names; a basis named later is read from that frame or refused, naming a data frame or list
- call 3, by dec-B684, which takes neither the design's refusal nor its alternatives: beside any data, a
  sampler's data included, a forests list states the forests in full; a forest declaring no basis has none,
  keeping a held basis means naming it, no list gives one forest; dec-B663's interim refusal is replaced
- the split into two landings, by dec-B679

Ruled for part 2: call 39 by dec-B685, call 36 by dec-B686, stan4bart's bundle by dec-B687 ((b)), call 41 by
dec-B688. Open, part 1: the calls its plan lists.

The calls below are agent calls: the design takes each so that it can be stated whole, and each is marked at the
session's close-out (dec-B629); a landing does not wait for it. They follow the rulings' direction: one shape for a user's object and a
sampler's data, base R's verbs, one source of truth with no derived copies, and matching by name.

Made in building dec-B678, agent calls as above:

19. **A user's object stores row names only where its frame cannot give them**: under a repeated subset that
    the na.action then cuts. Otherwise they are made by `make.unique()` of the frame's names at the positions,
    so there is no derived copy. Cost: about 10 lines. Alternative: the names always stored, 938 KB at 14,999
    rows.
23. **The response column holds the codes the sampler reads**, and `model.frame(s$data)` decodes them to the
    type and levels the build read. `$setResponse` takes numbers, so the column holds what it is given, and a
    reader gets the response as written in the formula. Cost: about 10 lines; 0.07 ms on read, none at the
    setter. Alternatives: the coded doubles shown as they are; or the factor stored, each `$setResponse`
    coding back, 0.07 ms a call.
24. **One frame column for each predictor column the engine reads**: each is a variable of the sampler that
    the setters address one at a time, so a matrix-valued term is split into its columns, named as the block
    names them, and `model.frame(s$data)` joins them back. Cost: within the resolution's lines. Alternative:
    matrix columns stored as `model.frame()` keeps them, a column's setter then rewriting that matrix.
25. **Under `factors = "indicators"` the frame holds the indicator columns**: the user asked for indicators as
    the predictors, and a 0/1 value given to one is taken, as in 0.9-34 and today. This departs from "factors
    as factors" for this coding only. Alternative: the factor held and the indicators made at creation, a
    write to one refused and the factor taken by name.
26. **A character predictor is held as a factor of its recorded levels**, coded once at the resolution and
    decoded on read. Alternative: held as characters and coded at each creation.
27. **The offsets are summed into one `(offset)` column**, since the sampler adds one; `$setOffset` then
    writes the total, which separate `offset()` columns could not take. The terms carry no `offset`
    attribute (call 38). Alternative: none that keeps `model.offset()` right after a setter.
28. **A forest's per-row columns are named by its label**, `(basis.<label>)` and
    `(forest.weights.<label>)`, matched by name as forests are. Alternatives: by position; or by the basis's
    expression, where two forests naming one expression would share a column that `$setForestBasis` must
    split.
29. **The test rows stay in `x.test` and their kin**: they come from `test`, and no model frame holds them.
    Alternative: a test frame, about 40 lines more.
30. **The bridge is handed today's layout, which R makes from the frame before each call**, sharing its
    columns, so the engine and the bridge's parsing are untouched. The engine pointer's protected value is
    then that layout: the bridge writes it at creation, at creation from a handle and at `$setData`, and
    reads it nowhere, the protection vector's only access being that write. It keeps the layout alive while
    the pointer lives, and `save()` writes it a second time. Whether the bridge must keep it at all is examined
    apart, because a comment in the R setter says the engine borrows `data@x` and one in the bridge says it
    retains no pointer. Until then the design keeps today's protection, and with it the second copy, 4.55 MB
    against 0.9-34's 2.32 MB. Alternative: the bridge parses the frame, about 150 lines under `src/`, with
    every parse site retested.
31. **The whole-matrix setter splits the matrix in the bridge**, about 15 lines, so no setter slows. In R the
    split takes 0.39 ms (0.37 to 0.84), against 0.08 ms for one copy, a call 5 percent slower and 11 percent
    of an iteration. Alternative: the split in R, one layer and nothing under `src/`, at that cost.
32. **On a sampler's data the training row names and `na.action` are the frame's**, as `model.frame()` gives
    them; the slots keep the test row names only. Alternative: the slots kept beside the frame.
33. **A sparse-bearing design keeps its sparse columns as one `dgCMatrix` in `x.sparse`**, as today's
    container does, and the frame holds the dense ones. Measured: `$setState` on a sparse design takes
    4.05 ms (p = 15) and 3.85 ms (p = 500); assembling the `dgCMatrix` from per-column `sparseVector`s at
    each bridge call would add 0.40 and 3.4 ms (10 and 88 percent), and base R cannot hold such a frame.
    `model.frame(s$data)` returns the dense columns and, after them, the sparse ones as the `sparseVector`
    and `sparseFactor` columns a data frame given to `dbartsData()` may hold; the help says base R neither
    prints nor models such a frame. The dense path is untouched. Alternatives: the per-column frame, at the
    measured cost; or the sparse columns made dense on read.
34. **An `aft` response is the time column and a `(status)` column**, and `model.frame(s$data)` rebuilds the
    `Surv` response from them, so `$setResponse` writes two vectors and builds no matrix. Alternative: a
    two-column matrix, 35 us more per `$setResponse`.
35. **A multinomial response is its count matrix**, and the trials are made from it where read. Alternative:
    the trials stored beside it, a second copy.
36. **A sampler's data's stand-in in `x` carries the coding record's attributes** (`term.labels`,
    `factor.levels`, `varTypes`, `terms`, `drop`), but not `dim`. So `attr()` reads answer as in 0.9-34, among
    them bartCause's term labels and the 16 internal ones, while `nrow()` and `[` still stop. dec-B676 rules a
    user's object's stand-in, not a sampler's data. Cost: about 5 lines. Alternatives: `NULL` in silence, as
    on a user's object; or the same attributes on a user's object too, one stand-in for both kinds, which
    revises dec-B676 and is the maintainer's. Ruled (dec-B686).
37. **A custom lint rule refuses `@x`, `@y`, `@weights` and `@offset` outside the accessors' file**, so every
    read is shown to have moved, now and after. Cost: about 25 lines. Alternative: a grep at review, which
    misses nothing written later only if someone runs it.
38. **The stored frame carries terms built from its own columns**: the response at position 1, its term
    labels its predictor columns, and no `offset` attribute. `model.frame(s$data)` carries the formula's terms
    with the `offset()` terms removed. Run, each reads back what the sampler runs on, where the formula's
    terms on the stored frame added a `poly()` column into the offset. Alternative: none that keeps
    `model.offset()` right.
39. **What still differs between `model.frame(o)` and `model.frame(s$data)`**, after the decoding of calls
    23, 24, 26 and 34:
    - the offsets, summed in `(offset)` where `model.frame(o)` has `offset()` columns too, since after a
      `$setOffset` only the total exists
    - under `factors = "indicators"`, the indicator columns in place of the factor (call 25)
    - a multinomial response, its count matrix where the user may have given a factor
    - integer and logical predictors, shown as doubles, since the setters write numbers
    - sparse columns, made on read (call 33)
    - the sampler's own columns: `(active)`, `(basis.<label>)`, `(forest.weights.<label>)`

    Each is what the sampler runs on and has no lossless way back. Alternative: each stored in the build's
    type and coded at each setter, a cost on every call. Ruled (dec-B685).
40. **A forest's per-row weights live in the frame**, `(forest.weights.<label>)`: they are per-row values a
    forest needs, so `$setData(s$data)`, a copy and a reload carry them from one place, and the sampler field
    goes. Cost: about 15 lines; a `$setForestWeights` costs one column write more. Alternative: the field as
    today, applied again on each re-creation.
41. **The active-row mask lives in the frame**, in `(active)` in every family. A probit or ordinal fit's 0/1
    weights ride the data already and survive `$setData(s$data)` (corrected 2026-11-08; the earlier premise
    said they became the mask); what `$setData(s$data)` drops today is a mask set by `$setActiveRows`, which
    it would then keep, and the sampler field goes. Cost: about 15 lines; `$setActiveRows` is 7 us more.
    Alternative: the field as today, `$setData(s$data)` dropping that mask. Ruled (dec-B688).
42. **The flat C API's setters do not write the frame**: they change the engine only, as they write no slot
    today, and the header does not change. `model.frame(s$data)` shows the last values given through R, and a
    copy or a reload rebuilds from them; the help says so. Alternative: an R-side mirror kept by the C entries,
    an ABI change and the maintainer's.
43. **A row cut clears the frame's `na.action`**, since the cut's own rows drop no missing value and the old
    record numbers another frame's rows. Alternative: rebuilding it through the cut, which nothing reads.

Also made here, agent calls as above:

12. Rows are stored as positions, with nothing stored when every row is kept in order. Alternative: always
    stored.
13. The internal row-cutting paths act on a sampler's data and cut its frame. Alternative: cutting a user's
    object's row record in step, which nothing would read.
14. A forest's record names columns, never numbers. Alternative: numbers, which reordered newdata breaks.
15. A forest is chosen by label or position; the matrix interface takes values only.
16. `bases =` gets R's own "unused argument" error, the argument never having shipped.
17. A record of fewer than two forests stays refused at creation, and the `$setData` refusal is reworded
    (2.8). Alternative: a record of one plain forest read as none.
18. A user's object keeps `x.test` as built: it comes from `test`, not the frame.
20. `model.frame(o)` is made in the build's order: evaluate on the whole frame, cut, code. Alternative:
    predict's path on the kept rows, which fails on a caller's vector once a row is dropped.
21. A fit keeps its own basis vector wherever it keeps per-forest fits, its sampler kept or not: it cannot be
    made again (dec-B185), and a kept sampler's basis is its current one. This answers the question dec-B671
    sent here. Alternative: read from the kept sampler, wrong after `fit$fit$setForestBasis()`.
22. No reader for a data object saved by 0.9-34, which already fails on main.

### Costs dec-B678 did not weigh

- **The size of the build.** The design was first sized at about 305 lines of R. dec-B678 removes `x`, `y`,
  `weights` and `offset` from a sampler's data, so every internal read of them moves to an accessor: about
  360 changed lines of R beside about 600 new, a lint rule, and up to about 740 test lines. No draw moves and
  no setter slows.
- **stan4bart's saved bundle.** The bundle keeps the specification's data with the block cut to no rows, and
  splices the fit's `bartData@x` back at a restore. Under dec-B678 that data holds the model frame, and
  neither the cut nor the splice can be written as a slot write. Shown:
  - (a) The bundle keeps the specification's data whole. 4 lines go from stan4bart; the bundle grows by the
    predictor columns, 0.050 to 0.088 MB at 750 rows and 25 columns, about n times p times 8 bytes in
    general, so 2.2 MB at 20,000 by 14.
  - (b) The bundle keeps no data, and a restore makes the specification again from `bartData` by the call the
    fit made. About 10 to 15 lines in stan4bart; the bundle is today's less the response and weights. It
    depends on `dbartsSpec()` deriving at a restore what it derived at the fit, as today's bundle already does
    (the cut points are derived at creation, and `setState(forceUpdate = TRUE)` puts the trees onto them).
  - (c) dbarts offers a way to make a sampler's data from a block-holding object and a response, and
    stan4bart calls it at a restore: about 15 lines in dbarts and 5 in stan4bart.

  Ruled (b) (dec-B687), in a reviewed stan4bart commit with the second landing.
- **A matrix-interface object kept beside its sampler** holds the predictors twice in a session, because the
  sampler's frame splits the matrix (2.4 MB at 20,000 by 15). The design takes this cost for call 24's
  setters. Alternative: on such a sampler, the matrix kept as one matrix column, with today's sharing and
  today's copy at each column's setter.
- **The second copy of a saved sampler's data** (call 30). This is a regression against 0.9-34 that predates
  the design, and call 30 keeps it until the engine question is settled.
