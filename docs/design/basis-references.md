# A forest's basis names a column of the data

Status: PROPOSED, 2026-11-07, revised three times 2026-11-07 for the maintainer's rulings. Follows dec-B663
(a forest's basis is a reference to a column of the data, and the per-forest bases a data object holds are the
package's record of that resolution), dec-B666 (the frame kept on the data object is the direction; the
numbering and the drop at the sampler to be examined) and dec-B670 (the whole frame is kept, and beside it as
little derived from it as can be: the predictor block is made from the frame when a sampler needs it).

## In short

A data object built from a formula and a data frame keeps the frame, whole and as given, the positions of the
rows it kept, and the record of how the predictors are coded, not the predictor block: that is made from the
frame when a sampler needs it, by the path `predict(newdata)` already takes. It keeps the response, the
weights, the offsets and what the caller gave beside the frame, as read. `keep.data = FALSE` keeps the block in
place of the frame, today's object. `dbartsSpec(o, forests = )` and `dbarts(o, forests = )` read a forest's
basis against the frame exactly as `dbarts(formula, data, forests = )` reads it. The resolution that makes a
sampler's specification makes the block and drops the frame and the rows, so a sampler, its data, a copy, a
saved sampler and a fit carry the block and no frame, as today. `dbartsData(bases = )` goes; the slots
`bases`, `basis.levels` and `basis.scale` stay, written by a sampler.

Examined, run on main:

- Rows are stored as their positions among the rows handed over, in the fit's order. That is the index the
  package already computes to cut a basis. The object records nothing when every row is kept in order.
- Columns: keeping the whole frame costs nothing in a session, since R shares its columns until one side
  changes one, and then copies only that column. Keeping fewer columns, or numbers pointing into the caller's
  frame, either breaks a basis that names a predictor or reads the wrong column after the caller edits the
  frame, and saves nothing in a session.
- The block made again from the frame's kept rows equals today's block in every case run (poly() and scale()
  under `subset`, na.omit, a sparse factor, a caller's variable changed after the build), in 0.0024 s at 20,000
  rows and 15 columns. Saved, a data object is then the frame and 0.25 MB, where frame and block would double it.
- Dropping at the sampler keeps a saved fit free of the frame, as dec-B185 has fits. It also keeps stan4bart's
  saved sampler bundle from growing by the frame.
- A data object with no frame to read names in (a sampler's own data, an opt-out object) keeps the frame's
  column names, and a basis naming one of them is refused, so a column is never read from the caller instead.
- The opt-out still serves bartCause, whose saved fits would otherwise carry the caller's whole frame, and lets
  stan4bart and bartCause keep reading the block off their data objects.
- A fit's own copy of each basis cannot be made again, the fit keeping no frame; it stays beside the
  per-forest fits it multiplies. A `keepFits = FALSE` fit stores none and prints no forest block today; the
  few numbers that block reads, stored per forest, would print it on every fit.

About 250 lines of R, nothing under `src/`, 83 test calls in 24 files and about 100 test lines reading the
block, five help pages, a NEWS sentence, a few lines in bartCause and one in stan4bart. Eleven calls for the
maintainer, listed at the end.

## 1. What exists

Run on main, 2026-11-07; the R code is unchanged between the library run and the branch point.

- **A formula's `forest()` term** (`bart()`, `dbarts()`), `forest(x2, basis = factor(z))`: code, its names
  found in `data` and then where the formula was written; `predict(fit, newdata)` reads `z` from `newdata` by
  name, in any column order, and refuses newdata without it.
- **A forests list beside a frame** (`dbarts(f, data = d, forests = )`): read when `forest()` is called; a
  column of `data` is that column and hides a caller's variable of the same name; a value is full length and
  cut by `subset` and the na.action.
- **The matrix interface**: a value, full length, cut the same way. A column of `x` is no name.
- **`dbartsData(bases = list(NULL, v))`**: values by forest position; the list fixes the forest count.
- **A data object handed to `dbarts()`**: carrying bases, with a list of knobs only, taken; a list declaring a
  basis, refused.
- **`dbartsSpec(o, forests = )`**: no frame, so a basis is the caller's value; it must be cut to the object's
  rows already, where `dbarts()` wants full length. Handed a sampler's own data, a declaring list replaces the
  object's bases whole, with no message.
- **`$setForestBasis`** by label or position; a value or a formula read where it was written.
- **`$setData`** is refused on a sampler of several forests, so it is no door for a basis.
- **Copy and reload** carry the blocks, the levels and the model's record of each code basis (its terms, its
  levels, its label and the columns it names), bit for bit.
- **0.9-34** had no forests, no bases and no `dbartsSpec()`.

bartCause's `bcf()` calls `dbartsData(bases = list(NULL, factor(z)))`, then `dbarts(o, forests = )` with knobs
only.

## 2. The design

### 2.1 Rows: what is numbered and how

The record is the number of rows handed over and the positions of the rows kept. It is the index `dbarts()`
already builds to cut a basis written beside a frame, so a value and a code basis are cut by one index. Read
off a basis that is the row number (12 rows; rows 3 and 9 missing a predictor):

- no `subset`, default na.action: 1 to 12, every row (`na.keepPredictors` keeps rows with a missing
  predictor)
- `subset = 2:10`: 2 to 10; `subset = -1`: 2 to 12; `subset = z == 1`: 5, 9, 10, 11, 12
- `subset = c(10, 2, 5, 7)`: 10, 2, 5, 7, in that order; `c(2, 2, 5, 7, 8)` keeps the repeat
- a character `subset` by row name: the positions of those names
- `na.action = na.omit`: rows 3 and 9 go; after `subset = 2:10` too: 2, 4, 5, 6, 7, 8, 10
- the matrix interface gives the same positions, and a full-length value is cut to them

Today the object records its row names and the na.action's drops counted within the subset, but not the subset,
so the positions cannot be rebuilt afterwards. They are recorded when the object is made, `subset` being read
once.

Storage: an integer vector, 4 bytes a row in a session. Saved, 1,000,000 positions take 3.0 MB shuffled and 2.0
MB in order with one gap. When every row is kept in order the record is empty: the row count alone, nothing
saved per row.

Keeping the frame cut to the kept rows would need no record, but it copies every column in a session, and a
code basis reads every row (run: `scale(x2)` under `subset = 1:6` is `scale(d$x2)[1:6]`, as in `lm`).

A later drop of rows. No user-facing door cuts a data object's rows: the sampler's setters keep the row count,
and `setActiveRows` masks rows without dropping them. Four internal paths do cut or expand them: the hurdle's
two parts, its zero part's merge, the hazard's person-period expansion and `xbart()`'s folds. None of them
resolves a basis afterwards, so each drops the frame and the record with the rows. The alternative cuts the
record in step at each, which nothing would read.

How a value is cut: a value has the length of the rows handed over and is cut by the positions. A value already
at the kept length on an object that dropped rows is refused, as `dbarts()` refuses one (2.5).

### 2.2 Columns: what is kept

**K1. The whole frame, as given.** `data` itself, every row and column, when it is a data frame or a list. An
environment is call 7. No `data` keeps nothing, and names are found where `forest()` was called, as today. A
matrix as `data` is refused by `dbartsData()` today (run: "data must be a data frame, list, or environment"),
so there is nothing to keep.

    o <- dbartsData(y ~ x1 + x2, data = d)
    s <- dbarts(o, forests = list(forest(), forest(x2, basis = factor(z))))

- In a session: nothing (run: the slot holds the same vectors as `d`). A later `d$z[1] <- 9` copies `z` alone;
  the other columns stay shared, and the object keeps the old `z`.
- Saved, 20,000 rows: with 3 predictors of 14 columns, the frame 1.19 MB, the data object today 0.45 MB, the
  two together 1.64 MB; with 13 predictors of 14 columns, 1.47, 1.47 and 2.94 MB. With the block made from
  the frame instead of kept (dec-B670, 2.3), the object is the frame and about 0.25 MB.
- A name that is no column is looked up where `forest()` was called, a value of full length; a name found
  nowhere is refused by name, as today.
- A column added to `d` after the object was made is not in it (run): the object is the data as it was handed
  over, as its predictors are.
- A basis naming a predictor (`forest(x1, basis = x2)`, which runs today) finds the predictor's column on
  every row, so `scale()` and `poly()` read every row, as at `dbarts()`.

**K2. The frame less the columns the formula reads.** Saved, nothing is written twice: stan4bart's frame
holds formula columns only. But a predictor is then no column of the kept frame. Its values on dropped rows
exist nowhere, and a basis naming it falls through to the caller, where a variable of the same name would be
taken in silence. Rebuilding it from `x` works only when no row was dropped and the column is plain numbers.
Not recommended.

**K3. Numbers pointing into the caller's frame**, a column number or the binding `d` read later. Nothing kept
or saved. Run: after the caller reorders `d`, column 11 is `x10` and no longer `z`; after `d <- d[-1, ]`, row 2
holds what was row 3. Both mean a silent misreading. A reload has nothing to point to, and a column added
later changes what an object made earlier fits. Not recommended.

**K4. Only columns named on the object** (`dbartsData(..., columns = ~ z)`, option B of the first note). Saved
small; the column is named twice, and the object states which columns are bases. Set aside by dec-B666.

**K1 with an opt-out.** As K1, and `dbartsData(..., keep.data = FALSE)` keeps no frame, as nlme's
`lme(keep.data = FALSE)` keeps no data and `lm(model = FALSE)` no model frame. The object then keeps the
predictor block, as today, and the frame's column names, and a basis naming one of them beside it is refused
(call 6). bartCause passes it:
its `bartc()` returns the response fit's data object as `data.rsp`, and its `bcf()` returns its own as `data`,
so under K1 alone every saved bartCause fit carries the caller's whole frame, columns the model never names
included, which dec-B185 keeps out of fits; slots being read-only (dec-B554), bartCause cannot strip it. Cost:
an argument and about 10 lines.

Ruled: K1 (dec-B670). Recommended: the opt-out, the frame kept by default (call 4; 2.3 on what it still
serves). In R, keeping the frame already is keeping a
reference: nothing is copied until one side changes it. Numbers into the caller's object would save nothing in
a session, and saved they have nothing to point to. The default keeps `dbartsSpec(dbartsData(f, d), L)` the
same as `dbarts(f, d, L)` for a user who never heard of the argument.

The forest's record of its basis names columns, not numbers. `predict(newdata)` with the columns in another
order gives the same prediction (run), so a name is what survives new data.

### 2.3 Beside the frame: what else is kept

dec-B670: beside the frame a data object keeps as little derived from it as it can. Saved sizes below are at
20,000 rows, 14 numeric columns and a factor of 26 levels: the frame 2.22 MB, today's data object 2.32 MB.

**The predictor block.** Today `x` is the columnar container (one entry a column, a factor kept a factor).
The bridge assembles from it the coded block the engine quantizes, inside sampler creation and `$setData`,
and frees it after. In a session the container shares the frame's double and factor columns (run); an integer
column, a cut to `subset` and indicator columns are copies. Saved, it writes the predictors again: 2.07 MB.

Made from the frame, it is what `predict(newdata)` already makes from new rows: the predictor terms, with
their data-dependent parts frozen (dec-B185: `poly()` coefficients, `scale()` centers, a caller's variable),
and the factor levels the block records. Run on the kept rows, `frame[rows, ]`, it gives today's block exactly:

- plain columns, and `y ~ .`
- `poly()` and `scale()` under `subset = 2:150`
- `na.omit`, and a missing response dropped
- a `subset` that leaves a factor level empty (the level kept)
- a sparse factor column, with and without `subset`
- a caller's vector named in the formula and changed after the build (the stored terms hold its value)
- indicator coding: the same values

Two adjustments. Predict refuses a missing value in a column that had none in training; on the training rows,
which are what defines that, the check is off (run: the default `na.keepPredictors` keeps such rows). And an
indicator matrix loses its attributes when cut to no rows, so the coding is stored apart.

Time: 0.0024 s at 20,000 rows and 15 columns, 0.062 s at 2,000 rows and 1,000 columns, beside 0.026 s to
create a sampler of 10 trees on the first.

Kept in its place, the coding record: the block with no rows, that is its column names, `varTypes`, term
labels, factor levels (or indicator levels and the dropped columns) and the frozen terms; 750 bytes saved.
Where it sits is call 9.

The block is made once, at the resolution, the step that drops the frame (2.4): `dbarts(o)`, `dbartsSpec(o)`,
`bart(o)`, `rbart_vi(o)`, `xbart(o)` (once at its start; its folds cut rows inside the engine), the hurdle's
entry, `new("dbartsSampler", data = o)` and `$setData(o)`. What comes out is today's object, the block and no
frame, so a sampler's data, its copy and reload, a fit's predict and stan4bart's per-iteration setters, which
act on the sampler, see no change. A saved `o` writes the frame and the record; the block is made at use.

**The other slots**, by source:

- The response's coding, kept with the record: `response.type`, `response.n.levels`, `response.levels`.
- `y`, the coded response, 0.15 MB. It can be made from the frame by the same replay, the response's
  expression frozen as the predictors' are (the stored terms drop it today) and coded again. Call 8.
- `weights` and `offset`, 0.10 and up to 0.15 MB: each a column of the frame or a value of the caller's
  (`weights = w`); the offset adds the formula's `offset()` terms. A value the frame does not hold is kept
  in any case. Call 8.
- `na.action`, the rows the na.action dropped, numbered within the subset: the positions alone do not give
  it, the subset not being kept apart. Kept (call 19).
- `rowNames`: the frame's row names at the kept rows, and the test rows'. Saved compactly for automatic row
  names; read through one accessor. Kept (call 19).
- From the caller's other arguments, kept as built: `x.test`, `weights.test`, `offset.test`, `n.cuts`, `sigma`,
  `counts` (the matrix interface only), `offset.category` and its test twin, `useQuantiles`,
  `testUsesRegularOffset`, `missing`. `x.test` comes from `test`, not `data` (call 18).
- Written by a sampler: `missing.seen`, `cutPoints`, `response.range`, `forest.calibration`,
  `standardization`, `status` on an aft sampler, `bases`, `basis.levels`, `basis.scale`. `NULL` on an object
  built from a frame, which no sampler holds (2.4); unchanged.

Saved, an object built from the frame is then the frame and 0.25 MB (2.47 MB), where frame and block together
would be 4.5 MB.

**An object with no frame** keeps the block, as today: the matrix interface (a data frame as `x` included),
`keep.data = FALSE`, a sampler's own data, and an environment as `data` (call 7). It also keeps the frame's
column names (call 6).

**Does the opt-out survive?** With the block gone, a saved object costs the frame alone, and
`keep.data = FALSE` keeps the block in place of the frame: today's object. It still serves two consumers.
bartCause's saved fits keep its data objects, and `bartc()` builds them from the caller's whole frame with a
column added, so they would carry columns the model never names, which dec-B185 keeps out of fits. And
stan4bart (about 6 places) and bartCause (about 12) read the block's values or row count off their data
objects; under the opt-out those keep working. Neither needs the frame: bartCause's `bcf()` gives its basis as
a value, stan4bart has none. Without the opt-out each moves its reads to an accessor (call 10), and bartCause
cuts the frame it hands over to the model's columns, which `subset` makes awkward (it may name other columns).
Call 4.

**What a user reads.** `o@x` on an object built from a frame no longer holds the rows. In 0.9-34 it held the
model matrix, factors as indicator columns (run); on main it is the columnar container, and the help sends a
reader to `as.matrix()` on it or to `extract(sampler, "predictors")`. Call 10.

**A fit's own copies.** A `bart` fit of several forests stores `bases` (each forest's basis at the training
rows) and `basis.terms` (the code of each basis written as code) only where it reports per forest, that is
under `keepFits = TRUE`, whatever `keepTrees` says (run). Their readers:

- `predict(bases = )` at the training rows and `extract()`'s per-forest contribution read the vector, with
  the per-forest fits it multiplies;
- `predict(newdata)` replays `basis.terms` on newdata and reads from `bases` only each forest's kind, levels
  and width;
- `print()` and `summary()` read, over the rows in the likelihood, a numeric basis's mean, sd, range and
  whether it is constant, and a factor's levels.

Can the vector be made again in place of kept?

- `keepTrees = FALSE`: no. The fit keeps no frame (dec-B185); `basis.terms` can be replayed only on rows
  handed to it.
- `keepTrees = TRUE`: the sampler's data holds the same block (run: equal), but as the sampler's current
  basis: after `fit$fit$setForestBasis()` it differs (run), and the fit's draws were taken on the first.
  Under D1 the sampler holds no frame either.
- A basis written as code against the frame, or given as values (bartCause's `bcf()`): alike; a value has no
  code to replay, and a code basis has no rows to replay on.

So the fit keeps it where it reports per forest. It is small beside what it serves: 16 KB against 6.4 MB of
per-forest fits at 2,000 rows and 200 draws.

`keepFits = FALSE` stores neither, nor the forest labels and tree counts, and prints no forest block at all
(run). The numbers the block reads, stored per forest at the fit (a few numbers, or a factor's levels, with
the labels and tree counts), serve `print()` and `summary()` on every fit, this one included, and
`predict(newdata)`'s kind, levels and width. The vector stays only beside the per-forest fits. Call 11.

### 2.4 Dropping at the sampler

**D1. Dropped where the specification is resolved**, in the step `dbarts()`, `dbartsSpec()` and `bart()`
share, and again by `new("dbartsSampler", ...)` for an object handed straight to it. The same step makes the
predictor block from the frame (2.3):

- `s$data` holds the predictors, the response and the record (blocks, levels, scale), with no frame and no row
  record, and the frame's column names. Its rows are its own. `dbartsSpec()`'s returned `data` is likewise not
  the object handed in: a caller holding `spec$data` starts a second model from the object it built.
- `dbarts(s$data)` and `dbartsSpec(s$data)` with no list, or a list of knobs only, re-create the forests from
  the record, as today.
- `dbartsSpec(s$data, forests = )` and `dbarts(s$data, forests = )` with a list declaring a basis, where
  `s$data` carries blocks: refused at both, naming both places (call 3). The original `o` still holds its
  frame, and a new model starts from it.
- The same beside a one-forest sampler's data, which carries no blocks: run today, `dbartsSpec(s1$data,
  forests = list(forest(), forest(basis = factor(z))))` takes the caller's `z` where `dbarts(f, d, L)` takes the
  column. This is call 6: a basis naming one of the frame's columns is refused, pointing to `o`; a name that
  is no column of the frame is the caller's value, as at `dbarts()`. The alternative widens call 3 to every
  sampler's data, refusing any declaring list beside it, which leaves an opt-out object uncovered.
- `$setData(o)` with an object that keeps its frame, on a sampler of one forest: taken, the frame dropped.
  Several forests: refused, as today.
- `$copy()`, `saveRDS` and `readRDS` of a sampler, and a `bart` fit with `keepSampler = TRUE`: no frame. Run:
  the fit saves at 1.49 MB.
- stan4bart: `dbartsSpec()`'s data, which its saved sampler bundle keeps with the design matrix stripped,
  has no frame, so the bundle does not grow by it (it still holds the response, the weights and the row
  names). Its fit keeps the data object made before the specification (`bartData`), which carries the BART
  frame: one more copy of columns its fit already saves in `frame`, unless it passes `keep.data = FALSE`.
- Memory: nothing on any sampler.

**D2. Kept on the sampler:**

- `s$data` carries the frame and the rows. A declaring list beside it can be resolved against the frame, the
  record replaced whole.
- `$setData(o)` brings `o`'s frame.
- A copy shares the frame. A saved sampler and a fit with `keepSampler` write it: a projection, main having no
  frame slot, from the fit and the frame saved together, 2.86 MB against 1.49. A fit then carries the caller's
  whole frame, which dec-B185 removed from fits for this reason.
- stan4bart's bundle gains the frame.
- What D2 adds is that `s$data` can start a new model, and the original `o` already can.

Recommended: D1. The rows then follow one rule. Creation takes a value of the length handed over;
`$setForestBasis` and `predict` with no newdata take the sampler's own rows, as today.

### 2.5 The doors, under K1 and D1

- **`dbartsData()`**: keeps the frame, the row record and the coding record, not the block (2.3);
  `bases =` gone; its formula still refuses a
  `forest()` term, pointing to a forests list.
- **`dbartsSpec(o, forests = )`** and **`dbarts(o, forests = )`**: a basis is read against `o`'s frame, then
  where `forest()` was called, built on `o`'s rows by the function `dbarts()` uses, so
  `dbartsSpec(dbartsData(f, d), forests = L)` and `dbarts(f, d, forests = L)` draw the same.
- **Full-length values at every creation door** (call 1). One rule: a value has the length of the rows handed
  over. The cost is that a kept-length value on an object that dropped rows, taken at `dbartsSpec()` today, is
  refused. Where no row was dropped the two lengths are one and nothing changes. The alternative takes either
  length at `dbartsSpec()`, two rules for one argument.
- **`dbarts(o)` taking a declaring list** (call 2). Recommended: taken, as at `dbartsSpec()`, so the two read
  one data object alike. bartCause's route then needs no data door. The alternative keeps `dbarts()` refusing,
  `dbartsSpec()` the only door that resolves beside a pre-built object.
- **A sampler's own data beside a declaring list** (call 3). Recommended: where it carries blocks, refused at
  both doors, naming both places; knobs only keep the record. Where it carries none, call 6. This is the
  interim refusal of dec-B663 made permanent and extended to `dbarts()`, whose broader refusal narrows to
  it. The alternatives are the list replacing the record whole, where a plain forest in it loses the basis and
  nothing is left to resolve a name against, or any forests list refused beside a record, knobs too.
- **`$setForestBasis`** and **`predict`**: as built. A forest is chosen by label or position; `predict(newdata)`
  reads each code basis by name; `bases =` stays for a basis that is no column of newdata and for the training
  rows.
- **The matrix interface**: no frame; a value, full length, cut by the record, at `dbarts(x, y)` and at
  `dbartsSpec(dbartsData(x, y))` alike.
- **Fit and reload**: nothing to add. `dbartsSpec()` and `dbarts(o)` write the model record `dbarts(f, d)`
  writes.

### 2.6 Data objects that carry blocks, and setData

Once `bases =` goes, a data object carries blocks only where a sampler of several forests wrote them: `s$data`
and what is made from it. A sampler of one forest writes none (run: its `s$data@bases` is `NULL`, and so is
that of `dbarts(forests = list(forest()))`). So:

- **`bases = list(NULL)`**, one forest with no basis, cannot be written. Today it is refused at `$setData` with
  a text saying the data "carry forest bases", though they state none, and at `dbarts()` with a text that
  advises `bases = list(NULL, z1, ...)` (run). Today `dbartsData()` makes such an object and the refusal comes
  at the sampler's creation; once `bases =` goes only a sampler writes a record, and the record's check at
  creation keeps refusing one of fewer than two forests, so no path makes one.
- **`$setData` on a sampler of one forest**, given another sampler's data that carry blocks: refused, as today,
  before anything is stored. The text says the data are those of a sampler of several forests, and points to
  the data object they were made from or to `dbarts()`; it no longer says "hand over data without 'bases'".
  Its own data handed back is taken (run).
- **`$setData` on a sampler of several forests**: refused whatever it is given, as today.
- **`dbarts(s$data)`, `dbartsSpec(s$data)` and a reload**: the record re-creates the forests (2.4).

Three texts advise the data door today (the `$setData` refusal, `dbarts()`'s refusal of a declaring list, and
the refusal of a model of one forest with a basis); each is reworded when this lands.

### 2.7 What else (A) needs

- **Where the drop sits.** In the resolution, not only in the sampler's creation, so that `dbartsSpec()`'s
  data, which consumers keep, carries no frame.
- **A data object saved by 0.9-34** already fails on main (run: "no slot of name "bases""). This design changes
  nothing there; reading one would be a call of its own, and no accessor is planned.
- **The slots**: proposed `frame` (the data as given, `NULL` for the matrix interface, after the drop and under
  the opt-out), `rows` (the positions, `NULL` when every row is kept in order), `frame.columns` (the frame's
  column names, kept after the drop and under the opt-out) and the coding record (call 9). `x` is empty where
  `frame` is set and the block where it is not. `dbartsData()` writes them; the resolution makes `x` and
  clears the first two; the help describes them so, read-only (dec-B554). The default show method prints the
  frame, as it prints `x` today.
- **The block's maker**: predict's path from new rows to the block, with its refusal of a missing value in a
  column that had none switched off for the training rows. One helper, called at each door of 2.3; every read
  of `x` before the resolution moves after it.
- **A comment to fix in the build**: the one above the row numbers `dbarts()` hands to the data object says a
  code basis is built only on the rows kept; it is read over every row and then cut (run, 2.1).
- **The record's check**: the positions lie within the rows handed over and are as many as the response; the
  frame has that many rows. The sampler's creation and the class's validity call it, as they call the check of
  the levels.
- **`dbartsData(o)`** handed a data object returns it as it is, frame and record included.
- **The test set** keeps no frame: a model of several forests refuses a test set.
- **Help**: the slots are described as the others a sampler writes are, read-only copies (dec-B554). The help
  shows no slot write as a door.

## 3. The landscape

Run where installed (mgcv 1.9.4, nlme 3.1.171, lme4 2.1.0, brms 2.23.0); stochtree and bcf read from source.

- `lm`, `glm`: a varying slope is a term; names found in `data`, then the formula's environment; the model frame
  keeps every variable named; `predict` fails without one.
- mgcv `s(x2, by = z)`, lme4 `(1 + z | g)`: the pre-built object (`gam(fit = FALSE)`, `lFormula()`) holds the
  model's columns; `predict` needs them.
- nlme: the fit keeps the whole data frame by default (`keep.data = TRUE`), columns the model never names
  included. This is the precedent for K1.
- brms: by-variables must be in `data` or `data2`.
- stochtree, bcf (Hahn, Murray and Carvalho): values only. bartCause's `bcf()`: a column of `data` or a
  vector. BART and bartMachine: no basis.

## 4. What a user meets

- A formula user, a forests-list user beside a frame, a matrix-interface user of `dbarts()`: nothing moves.
- treatSens: nothing moves.
- stan4bart: its data object reads the block's values in about 6 places, so it passes `keep.data = FALSE`
  (one line) and nothing else moves; its fit and bundle sizes are as today.
- bartCause: `bases =` goes, so the second forest takes `basis =` with the full-length factor it already builds,
  and its data objects pass `keep.data = FALSE`, so that its saved fits (`data.rsp`, `data`) carry no frame and
  its about 12 reads of the block keep working (about 10 lines). Draws identical (run: the data door, a list
  on `dbarts()` and on `dbartsSpec()` agree, with and without `subset`); its forest's label becomes the
  basis's text, and it reads no label.
- A `dbartsSpec()` caller with a basis: a name that is a column now finds it. Where it is also a caller's
  variable, the column wins, as at `dbarts()`: the one call whose draws move. A kept-length value on an object
  that dropped rows is refused. Beside a sampler's own data a name of a frame column is refused (call 6).
- A reader of `o@x` on an object built from a frame finds no rows there and calls `extract(o, "predictors")`
  (call 10). In 0.9-34 the slot held the model matrix.
- Anyone who saves a data object built from a frame: it now holds the frame in place of the block, about the
  same size where the frame is mostly the model's columns, larger where it holds others.
- A `keepFits = FALSE` fit of several forests prints its forest block, which it does not today (call 11).

## 5. The cost

- R, about 250 lines, no `src/` or header change:
  - the slots and their check, about 15
  - `dbartsData()` keeping the frame and recording the rows, about 30
  - `bases =` leaving the exported function, about 20
  - `dbartsSpec()` resolving against the frame, with the refusal, about 30
  - `dbarts(o)` taking a declaring list, about 25
  - the drop at the resolution and at the four internal row-cutting paths, about 10
  - the opt-out and the refusal of a frame column's name beside a frameless object, about 15
  - the block made from the frame: the maker, the coding record, the calls at each door, reads of `x` moved
    after the resolution, the validity of an object with no block, about 75
  - `extract(o, "predictors")`, about 10
  - a fit's printed numbers stored per forest, the labels and tree counts on every fit, about 20
- Tests: 83 calls of `dbartsData(bases = )` in 24 files: 45 fixtures become forests lists drawing what they
  drew, 36 pins of the data door go where the forest door pins the same or move there, 2 oracles are rewritten
  on the list's value. About 100 lines in 15 files read `x` off an object built from a formula; they move to
  `extract()` or `keep.data = FALSE`. New: the sibling identity, the row positions under each `subset` and
  na.action above, the frame gone from `s$data` and from `dbartsSpec()`'s data, the block made from the frame
  equal to today's in each case of 2.3, a column changed in `d` after the object was made, full-length values,
  the refusal of a sampler's data beside a declaring list, an opt-out object and a one-forest sampler's data
  refusing a frame column's name, an environment as `data`, and a `keepFits = FALSE` fit's block.
- Help: `dbartsData` (the slots, `keep.data`), `dbartsSpec`, `dbarts`, `dbartsSampler-class` and `extract`.
  Four design notes that name the data door (forest-kind-by-class, public-surface, written-surface,
  multiplier-combiner) change at landing.
- NEWS: `bases =` never shipped and needs nothing. Against 0.9-34, one entry:
  "dbartsData() given a formula and a data frame keeps the frame and the positions of the rows it keeps in
  place of the predictor matrix, which extract(object, "predictors") makes; a forest given to dbarts() or
  dbartsSpec() beside it finds its basis there. keep.data = FALSE keeps the matrix instead, as before."
- bartCause: about 10 lines and its suite. stan4bart: one line and its suite.

## Calls made

For the maintainer:

1. Full-length values at every creation door, a cut one refused (2.5). Recommended. The alternative takes
   either length at `dbartsSpec()`.
2. `dbarts(o)` takes a declaring list (2.5). Recommended. The alternative keeps it refused.
3. A sampler's own data carrying blocks, beside a declaring list, refused at both doors (2.5). Recommended. The
   alternatives are the record replaced whole, or any list refused.
4. The frame is ruled (dec-B670); open are the opt-out and its default. Recommended: `keep.data = FALSE`
   kept, the frame kept by default; the opt-out object keeps the block in place of the frame, as today, and
   stan4bart and bartCause pass it (2.3). The alternatives are no switch, under which bartCause's saved fits
   carry the caller's frame and both consumers move about 18 reads to `extract()`, and the flip, the frame
   left out by default, under which `dbartsSpec(dbartsData(f, d), L)` finds no column `dbarts(f, d, L)` finds
   unless asked.
5. The frame dropped at the resolution, which makes the block (D1, 2.4). Recommended. The alternative is D2.
6. A data object with no frame to read names in (a sampler's own data, an opt-out object) keeps the frame's
   column names, and a basis naming one is refused, pointing to the data object built with its frame; other
   names are the caller's values (2.4). Recommended. The alternatives are widening call 3 to any declaring
   list beside a sampler's data, which leaves an opt-out object reading from the caller, or reading every
   name from the caller, a column silently replaced by a caller's variable of its name.
7. An environment given as `data` is kept until the resolution, so that `dbartsSpec(dbartsData(f, e), L)`
   reads a basis's name where `dbarts(f, e, L)` reads it (run: `dbarts()` takes `e$z` over a local `z`), and
   the block is made at the build and kept beside it: an environment is no frame, and a variable changed in
   it later would otherwise change the predictors. Recommended. Cost: a basis's name is read from the
   environment at the resolution, so a change made after the build is read there; a saved data object writes
   the environment. The alternatives are a copy of its variables taken at the build (read once, as a frame
   is; saved, every variable it holds), its variables' names kept and a basis naming one refused, as in call
   6, or the block too made at the resolution from it.
8. `y`, `weights` and `offset` kept as read at the build (2.3). Recommended, though `y` and a column-named
   weight could be made from the frame: each is one column (0.15 MB of a 2.47 MB object), every door reads the
   row count and the response off them (about 250 reads inside dbarts, about 10 in bartCause), the response's
   refusals stay at the build, and a weight or offset given as a caller's value must be kept anyway. The
   alternatives are `y` made at the resolution too (about 25 lines: the response's expression frozen, coded
   again), and also `weights` and `offset` where they name columns.
9. On an object built from a frame, `x` is empty and the coding record has its own slot (2.3). Recommended:
   a reader of the block's rows fails at once. The alternative keeps the block with no rows in `x`, so
   `colnames(o@x)` and its levels keep working while `nrow(o@x)` gives 0 in silence (bartCause writes
   `n.obs <- nrow(responseData@x)`).
10. `extract(o, "predictors")` makes the block of a data object, as `extract(sampler, "predictors")` gives a
   sampler's (about 10 lines). Recommended. The alternatives are a `model.matrix()` method, base R's name for
   a model's design matrix, or no accessor.
11. Every fit of several forests stores, per forest, the numbers its printed block reads (a numeric basis's
   mean, sd, range and constancy over the rows in, a factor's levels), with the labels and tree counts; the
   basis itself stays only where the fit keeps per-forest fits, for `predict(bases = )` at the training rows
   and `extract()`'s contribution (2.3). Recommended: `print()` and `summary()` then work on every fit,
   `keepFits = FALSE` included, for a few numbers a forest. The alternatives are today's, a `keepFits = FALSE`
   fit printing no block, or the basis kept on every fit, n values a forest, to print four numbers.

Made here:

12. Rows stored as positions, with nothing stored when every row is kept in order. The alternative stores the
   positions always, 4 bytes a row.
13. The four internal row-cutting paths make the block and drop the frame and the record. The alternative
   cuts the record in step.
14. A forest's record names columns, never numbers. The alternative numbers them, which reordered newdata
   breaks.
15. A forest is still chosen by label or position, and the matrix interface takes values only, as in the first
   note.
16. `bases =` gets R's own "unused argument" error, the argument never having shipped. The alternative is a
   message of our own.
17. A record of fewer than two forests stays refused at creation, so `bases = list(NULL)` has no path, and the
   `$setData` refusal is reworded for another sampler's data (2.6). The alternative takes a record of one
   plain forest as no record, which nothing can write once `bases =` goes.
18. `x.test` kept as built: it comes from `test`, not `data`, and a model of several forests refuses a test
   set, so no basis reads it. The alternative keeps the test frame and makes it at the resolution too.
19. `na.action` and `rowNames` kept: the first is not given by the positions alone, the second is saved
   compactly. The alternative makes the training row names from the frame through their one accessor.
20. The block is made by predict's path, not by a second run of the build. The alternative re-runs
   `model.frame()` and the coding at each door, two paths that must agree.
21. A fit keeps its basis vector where it reports per forest; it cannot be made again, having no frame
   (dec-B185), and a kept sampler's block is the sampler's current basis (2.3). The alternative reads it from
   the kept sampler under `keepTrees = TRUE`, wrong after `fit$fit$setForestBasis()`.
22. No accessor for a data object saved by 0.9-34, which already fails on main. The alternative reads one,
   about 5 lines, a call of its own.
