# A forest's basis names a column of the data

Status: PROPOSED, 2026-11-07, revised twice 2026-11-07 for the maintainer's ruling. Follows dec-B663 (a forest's
basis is a reference to a column of the data, and the per-forest bases a data object holds are the package's
record of that resolution) and dec-B666 (the frame kept on the data object is the direction; the numbering and
the drop at the sampler to be examined).

## In short

A data object keeps the data it was built from, whole and as given, and the positions of the rows it kept;
`keep.data = FALSE` leaves the frame out.
`dbartsSpec(o, forests = )` and `dbarts(o, forests = )` read a forest's basis against them exactly as
`dbarts(formula, data, forests = )` reads it. The resolution that makes a sampler's specification drops both,
so a sampler, its data, a copy, a saved sampler and a fit carry neither. `dbartsData(bases = )` goes; the
slots `bases`, `basis.levels` and `basis.scale` stay, written by a sampler.

Examined, run on main:

- Rows are stored as their positions among the rows handed over, in the fit's order. That is the index the
  package already computes to cut a basis. The object records nothing when every row is kept in order.
- Columns: keeping the whole frame costs nothing in a session, since R shares its columns until one side
  changes one, and then copies only that column. Saved, it costs the frame: a saved data object roughly
  doubles where the frame is mostly predictors. Keeping fewer columns, or numbers pointing into the caller's
  frame, either breaks a basis that names a predictor or reads the wrong column after the caller edits the
  frame, and saves nothing in a session. An opt-out, as nlme's `lme(keep.data = )`, serves a consumer that
  keeps its data object: bartCause's saved fits would otherwise carry the caller's whole frame.
- Dropping at the sampler keeps a saved fit free of the frame, as dec-B185 has fits. It also keeps stan4bart's
  saved sampler bundle from growing by the frame.
- A data object with no frame to read names in (a sampler's own data, an opt-out object) keeps the frame's
  column names, and a basis naming one of them is refused, so a column is never read from the caller instead.

About 145 lines of R, nothing under `src/`, 83 test calls in 24 files, four help pages, a NEWS sentence, about
8 lines in bartCause. Seven calls for the maintainer, listed at the end.

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
at the kept length on an object that dropped rows is refused, as `dbarts()` refuses one (2.4).

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
  two together 1.64 MB; with 13 predictors of 14 columns, 1.47, 1.47 and 2.94 MB. So a saved data object
  grows by the frame, about double where the frame is mostly predictors.
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
frame's column names alone, and a basis naming one of them beside it is refused (call 6). bartCause passes it:
its `bartc()` returns the response fit's data object as `data.rsp`, and its `bcf()` returns its own as `data`,
so under K1 alone every saved bartCause fit carries the caller's whole frame, columns the model never names
included, which dec-B185 keeps out of fits; slots being read-only (dec-B554), bartCause cannot strip it. Cost:
an argument and about 10 lines.

Recommended: K1 with the opt-out, kept by default (call 4). In R, keeping the frame already is keeping a
reference: nothing is copied until one side changes it. Numbers into the caller's object would save nothing in
a session, and saved they have nothing to point to. The default keeps `dbartsSpec(dbartsData(f, d), L)` the
same as `dbarts(f, d, L)` for a user who never heard of the argument. The alternative is K1 without a switch,
bartCause's fits carrying the frame.

The forest's record of its basis names columns, not numbers. `predict(newdata)` with the columns in another
order gives the same prediction (run), so a name is what survives new data.

### 2.3 Dropping at the sampler

**D1. Dropped where the specification is resolved**, in the step `dbarts()`, `dbartsSpec()` and `bart()`
share, and again by `new("dbartsSampler", ...)` for an object handed straight to it:

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

### 2.4 The doors, under K1 and D1

- **`dbartsData()`**: keeps the frame and the row record; `bases =` gone; its formula still refuses a
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

### 2.5 Data objects that carry blocks, and setData

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
- **`dbarts(s$data)`, `dbartsSpec(s$data)` and a reload**: the record re-creates the forests (2.3).

Three texts advise the data door today (the `$setData` refusal, `dbarts()`'s refusal of a declaring list, and
the refusal of a model of one forest with a basis); each is reworded when this lands.

### 2.6 What else (A) needs

- **Where the drop sits.** In the resolution, not only in the sampler's creation, so that `dbartsSpec()`'s
  data, which consumers keep, carries no frame.
- **A data object saved by 0.9-34** already fails on main (run: "no slot of name "bases""). This design changes
  nothing there; reading one would be a call of its own, and no accessor is planned.
- **The slots**: proposed `frame` (the data as given, `NULL` for the matrix interface, after the drop and under
  the opt-out), `rows` (the positions, `NULL` when every row is kept in order) and `frame.columns` (the frame's
  column names, kept after the drop and under the opt-out). `dbartsData()` writes them and the resolution
  clears the first two; the help describes them so, read-only (dec-B554). The default show method prints the
  frame, as it prints `x` today.
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
- stan4bart, treatSens: nothing moves. stan4bart's saved fit carries its BART frame once more unless it
  passes `keep.data = FALSE` (2.3).
- bartCause: `bases =` goes, so the second forest takes `basis =` with the full-length factor it already builds,
  and its data objects pass `keep.data = FALSE` so that its saved fits (`data.rsp`, `data`) carry no frame
  (about 8 lines). Without the opt-out, every saved `bartc()` and `bcf()` fit grows by the caller's frame.
  Draws identical (run: the data door, a list on `dbarts()` and on `dbartsSpec()` agree, with and without
  `subset`); its forest's label becomes the basis's text, and it reads no label.
- A `dbartsSpec()` caller with a basis: a name that is a column now finds it. Where it is also a caller's
  variable, the column wins, as at `dbarts()`: the one call whose draws move. A kept-length value on an object
  that dropped rows is refused. Beside a sampler's own data a name of a frame column is refused (call 6).
- Anyone who saves a data object: it now includes the frame, unless built with `keep.data = FALSE`.

## 5. The cost

- R, about 145 lines, no `src/` or header change:
  - the slots and their check, about 15
  - `dbartsData()` keeping the frame and recording the rows, about 30
  - `bases =` leaving the exported function, about 20
  - `dbartsSpec()` resolving against the frame, with the refusal, about 30
  - `dbarts(o)` taking a declaring list, about 25
  - the drop at the resolution and at the four internal row-cutting paths, about 10
  - the opt-out and the refusal of a frame column's name beside a frameless object, about 15
- Tests: 83 calls of `dbartsData(bases = )` in 24 files: 45 fixtures become forests lists drawing what they
  drew, 36 pins of the data door go where the forest door pins the same or move there, 2 oracles are rewritten
  on the list's value. New: the sibling identity, the row positions under each `subset` and na.action above,
  the frame gone from `s$data` and from `dbartsSpec()`'s data, a column changed in `d` after the object was
  made, full-length values, the refusal of a sampler's data beside a declaring list, an opt-out object and a
  one-forest sampler's data refusing a frame column's name, and an environment as `data`.
- Help: `dbartsData`, `dbartsSpec`, `dbarts` and `dbartsSampler-class`. Four design notes that name the data
  door (forest-kind-by-class, public-surface, written-surface, multiplier-combiner) change at landing.
- NEWS: `bases =` never shipped and needs nothing. Against 0.9-34 a user meets one new fact, one sentence:
  "dbartsData() keeps the data frame it is given and the positions of the rows it keeps, so that a forest
  given to dbarts() or dbartsSpec() beside it finds its basis there; a saved data object therefore includes
  the frame, which keep.data = FALSE leaves out."
- bartCause: about 8 lines and its suite.

## Calls made

For the maintainer:

1. Full-length values at every creation door, a cut one refused (2.4). Recommended. The alternative takes
   either length at `dbartsSpec()`.
2. `dbarts(o)` takes a declaring list (2.4). Recommended. The alternative keeps it refused.
3. A sampler's own data carrying blocks, beside a declaring list, refused at both doors (2.4). Recommended. The
   alternatives are the record replaced whole, or any list refused.
4. The whole frame kept, with an opt-out, `keep.data = FALSE`, kept by default (K1, 2.2). Recommended. The
   alternatives are K1 without a switch (bartCause's saved fits carrying the caller's frame), K2, K3 and K4.
5. The frame dropped at the resolution (D1, 2.3). Recommended. The alternative is D2.
6. A data object with no frame to read names in (a sampler's own data, an opt-out object) keeps the frame's
   column names, and a basis naming one is refused, pointing to the data object built with its frame; other
   names are the caller's values (2.3). Recommended. The alternatives are widening call 3 to any declaring
   list beside a sampler's data, which leaves an opt-out object reading from the caller, or reading every
   name from the caller, a column silently replaced by a caller's variable of its name.
7. An environment given as `data` is kept, the environment itself, until the resolution drops it, so that
   `dbartsSpec(dbartsData(f, e), L)` reads what `dbarts(f, e, L)` reads (run: `dbarts()` takes `e$z` over a
   local `z`, which the first note's choice, keeping nothing, would have taken instead). Recommended. Cost: a
   name is read from the environment at the resolution, so a change made to it after the data object was built
   is read; a saved data object writes the environment. The alternatives are a copy of its variables taken at
   the build (read once, as a frame is; saved, every variable it holds) or its variables' names kept and a
   basis naming one refused, as in call 6.

Made here:

8. Rows stored as positions, with nothing stored when every row is kept in order. The alternative stores the
   positions always, 4 bytes a row.
9. The four internal row-cutting paths drop the frame and the record. The alternative cuts the record in step.
10. A forest's record names columns, never numbers. The alternative numbers them, which reordered newdata
   breaks.
11. A forest is still chosen by label or position, and the matrix interface takes values only, as in the first
   note.
12. `bases =` gets R's own "unused argument" error, the argument never having shipped. The alternative is a
   message of our own.
13. A record of fewer than two forests stays refused at creation, so `bases = list(NULL)` has no path, and the
   `$setData` refusal is reworded for another sampler's data (2.5). The alternative takes a record of one
   plain forest as no record, which nothing can write once `bases =` goes.
14. No accessor for a data object saved by 0.9-34, which already fails on main. The alternative reads one,
   about 5 lines, a call of its own.
