# A forest's basis names a column of the data

Status: PROPOSED, 2026-11-07, revised 2026-11-07 for the maintainer's ruling. Follows dec-B663 (a forest's
basis is a reference to a column of the data, and the per-forest bases a data object holds are the package's
record of that resolution) and dec-B666 (the frame kept on the data object is the direction; the numbering and
the drop at the sampler to be examined).

## In short

A data object keeps the data it was built from, whole and as given, and the numbers of the rows it kept.
`dbartsSpec(o, forests = )` and `dbarts(o, forests = )` read a forest's basis against them exactly as
`dbarts(formula, data, forests = )` reads it. The resolution that makes a sampler's specification drops both,
so a sampler, its data, a copy, a saved sampler and a fit carry neither. `dbartsData(bases = )` goes; the
slots `bases`, `basis.levels` and `basis.scale` stay, written by a sampler.

Examined, run on main:

- Rows are stored as their positions among the rows handed over, in the fit's order. That is the index the
  package already computes to cut a basis. The object records nothing when every row is kept in order.
- Columns: keeping the whole frame costs nothing in a session, since R shares its columns until one side
  changes one, and then copies only that column. Saved, it costs the frame. Keeping fewer columns, or numbers
  pointing into the caller's frame, either breaks a basis that names a predictor or reads the wrong column
  after the caller edits the frame, and saves nothing in a session.
- Dropping at the sampler keeps a saved fit at 1.49 MB where keeping would make it 2.86 MB (a 20,000-row
  frame). It also keeps stan4bart's saved sampler bundle independent of the row count, as it is built to be.

About 135 lines of R, nothing under `src/`, 83 test calls in 24 files, four help pages, about 6 lines in
bartCause. Five calls for the maintainer, listed at the end.

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

A later drop of rows. No user-facing door cuts a data object's rows: the sampler's setters keep the row count,
and `setActiveRows` masks rows without dropping them. Four internal paths do cut or expand them: the hurdle's
two parts, its zero part's merge, the hazard's person-period expansion and `xbart()`'s folds. None of them
resolves a basis afterwards, so each drops the frame and the record with the rows. The alternative cuts the
record in step at each, which nothing would read.

How a value is cut: a value has the length of the rows handed over and is cut by the positions. A value already
at the kept length on an object that dropped rows is refused, as `dbarts()` refuses one (2.4).

### 2.2 Columns: what is kept

**K1. The whole frame, as given.** `data` itself, every row and column, when it is a data frame or a list. An
environment, or no `data`, keeps nothing, and names are found where `forest()` was called, as today.

    o <- dbartsData(y ~ x1 + x2, data = d)
    s <- dbarts(o, forests = list(forest(), forest(x2, basis = factor(z))))

- In a session: nothing (run: the slot holds the same vectors as `d`). A later `d$z[1] <- 9` copies `z` alone;
  the other columns stay shared, and the object keeps the old `z`.
- Saved (20,000 rows, 14 columns): the frame 1.19 MB; the data object today 0.45 MB; the two together 1.64 MB.
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

Recommended: K1. In R, keeping the frame already is keeping a reference: nothing is copied until one side
changes it. Numbers into the caller's object would save nothing in a session, and saved they have nothing to
point to.

The forest's record of its basis names columns, not numbers. `predict(newdata)` with the columns in another
order gives the same prediction (run), so a name is what survives new data.

### 2.3 Dropping at the sampler

**D1. Dropped where the specification is resolved**, in the step `dbarts()`, `dbartsSpec()` and `bart()`
share, and again by `new("dbartsSampler", ...)` for an object handed straight to it:

- `s$data` holds the predictors, the response and the record (blocks, levels, scale), with no frame and no row
  record. Its rows are its own.
- `dbarts(s$data)` and `dbartsSpec(s$data)` with no list, or a list of knobs only, re-create the forests from
  the record, as today.
- `dbartsSpec(s$data, forests = )` and `dbarts(s$data, forests = )` with a list declaring a basis: refused at
  both, naming both places (2.4). The original `o` still holds its frame, and a new model starts from it.
- `$setData(o)` with an object that keeps its frame, on a sampler of one forest: taken, the frame dropped.
  Several forests: refused, as today.
- `$copy()`, `saveRDS` and `readRDS` of a sampler, and a `bart` fit with `keepSampler = TRUE`: no frame. Run:
  the fit saves at 1.49 MB.
- stan4bart: `dbartsSpec()`'s data, which its saved sampler bundle keeps with the design matrix stripped so
  that the bundle does not grow with the rows, has no frame. Its fit keeps the data object made before the
  specification (`bartData`), which carries the BART frame: one more copy of columns its fit already saves in
  `frame`, a fraction of a fit its comments put at tens of MB of trees.
- Memory: nothing on any sampler.

**D2. Kept on the sampler:**

- `s$data` carries the frame and the rows. A declaring list beside it can be resolved against the frame, the
  record replaced whole.
- `$setData(o)` brings `o`'s frame.
- A copy shares the frame. A saved sampler and a fit with `keepSampler` write it: 2.86 MB against 1.49 in the
  run, a fit carrying the caller's whole frame, which dec-B185 removed from fits for this reason.
- stan4bart's bundle gains the frame and grows with the rows, against its own comment.
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
- **A sampler's own data beside a declaring list** (call 3). Recommended: refused at both doors, naming both
  places; knobs only keep the record. This is the interim refusal of dec-B663 made permanent and extended to
  `dbarts()`, whose broader refusal narrows to it. The alternatives are the list replacing the record whole,
  where a plain forest in it loses the basis and nothing is left to resolve a name against, or any forests
  list refused beside a record, knobs too.
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
  advises `bases = list(NULL, z1, ...)` (run). The record's check also refuses a record of fewer than two
  forests, so no path makes one.
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
- **A data object saved before these slots existed** is read through accessors that answer "none" for a missing
  slot, as the levels are read today; 0.9-34's objects then load and fit as before.
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
- stan4bart, treatSens: nothing moves. stan4bart's saved fit carries its BART frame once more (2.3).
- bartCause: `bases =` goes, so the second forest takes `basis =` with the full-length factor it already builds
  (about 6 lines). Draws identical (run: the data door, a list on `dbarts()` and on `dbartsSpec()` agree, with
  and without `subset`); its forest's label becomes the basis's text, and it reads no label.
- A `dbartsSpec()` caller with a basis: a name that is a column now finds it. Where it is also a caller's
  variable, the column wins, as at `dbarts()`: the one call whose draws move. A kept-length value on an object
  that dropped rows is refused.

## 5. The cost

- R, about 135 lines, no `src/` or header change:
  - the slots and their check, about 15
  - `dbartsData()` keeping the frame and recording the rows, about 30
  - `bases =` leaving the exported function, about 20
  - `dbartsSpec()` resolving against the frame, with the refusal, about 30
  - `dbarts(o)` taking a declaring list, about 25
  - the drop at the resolution and at the four internal row-cutting paths, about 10
  - the accessors for old objects, about 5
- Tests: 83 calls of `dbartsData(bases = )` in 24 files: 45 fixtures become forests lists drawing what they
  drew, 36 pins of the data door go where the forest door pins the same or move there, 2 oracles are rewritten
  on the list's value. New: the sibling identity, the row positions under each `subset` and na.action above,
  the frame gone from `s$data` and from `dbartsSpec()`'s data, a column changed in `d` after the object was
  made, full-length values, and the refusal of a sampler's data beside a declaring list.
- Help: `dbartsData`, `dbartsSpec`, `dbarts` and `dbartsSampler-class`. Four design notes that name the data
  door (forest-kind-by-class, public-surface, written-surface, multiplier-combiner) change at landing.
- NEWS: 0.9-34 had no bases; at most a clause in the multi-forest item.
- bartCause: about 6 lines and its suite.

## Calls made

For the maintainer:

1. Full-length values at every creation door, a cut one refused (2.4). Recommended. The alternative takes
   either length at `dbartsSpec()`.
2. `dbarts(o)` takes a declaring list (2.4). Recommended. The alternative keeps it refused.
3. A sampler's own data beside a declaring list refused at both doors (2.4). Recommended. The alternatives
   are the record replaced whole, or any list refused.
4. The whole frame kept (K1, 2.2). Recommended. The alternatives are K2, K3 and K4.
5. The frame dropped at the resolution (D1, 2.3). Recommended. The alternative is D2.

Made here:

6. Rows stored as positions, with nothing stored when every row is kept in order. The alternative stores the
   positions always, 4 bytes a row.
7. The four internal row-cutting paths drop the frame and the record. The alternative cuts the record in step.
8. An environment given as `data` is not kept. The alternative keeps it, a reference that sees later changes.
9. A forest's record names columns, never numbers. The alternative numbers them, which reordered newdata
   breaks.
10. A forest is still chosen by label or position, and the matrix interface takes values only, as in the first
   note.
11. `bases =` gets R's own "unused argument" error, the argument never having shipped. The alternative is a
   message of our own.
12. A record of fewer than two forests is refused by the record's check, so `bases = list(NULL)` has no
   path, and the `$setData` refusal is reworded for another sampler's data (2.5). The alternative takes a
   record of one plain forest as no record, which nothing can write once `bases =` goes.
