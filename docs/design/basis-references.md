# A forest's basis names a column of the data

Status: PROPOSED, 2026-11-07, for the maintainer's ruling. Follows dec-B663: a forest's basis is a reference
to a column of the data, resolved by name where a data frame is given, and the per-forest bases a data object
holds are the package's record of that resolution.

## In short

Recommended: a data object keeps the data frame it was built from and the numbers of the rows it kept, so that
a forest handed to `dbartsSpec()` or `dbarts()` beside it finds its basis there, as `dbarts()` finds it when
given the formula and the frame. `dbartsData(bases = )` goes; the slots `bases`, `basis.levels` and
`basis.scale` stay, written by a sampler. Everything else stays as built: the formula's `forest()` term, a
forests list, `$setForestBasis`, `predict`, the matrix interface, copy and reload. No draw moves for a call
taken both before and after, but one: a basis whose name is both a caller's variable and a column of the frame
(section 4). The refusal dec-B663 orders for the meantime is still needed afterwards, narrowed to a sampler's
own data handed back (section 6). About 110 lines of R, 83 test calls in 24 files, four help pages and about 6
lines in bartCause.

One premise of dec-B663 is wrong: a consumer does call `dbartsData(bases = )`. bartCause's `bcf()` builds its
data object with `bases = list(NULL, factor(z))` and hands it to `dbarts()` with a forests list that states
only knobs.

## 1. What exists

Run on main on 2026-11-07. The R code is unchanged between the library run and the branch point.

**The formula's `forest()` term** (`bart()`, `dbarts()`).

    bart(y ~ x1 + x2 + forest(x2, basis = factor(z)), data = d)

The basis is code. Its names are found in `data` and then where the formula was written, when the model is
built. `predict(fit, newdata)` reads `z` from `newdata` by name, and refuses newdata without it, naming the
forest. A caller's vector in a term (`basis = zv`) is a name too: `predict` then asks `newdata` for a column
`zv`.

**A forests list with a data frame** (`dbarts()`).

    dbarts(y ~ x1 + x2, data = d, forests = list(forest(), forest(x2, basis = factor(z))))

Read when `forest()` is called: a name of a column of `data` is that column, and hides a caller's variable of
the same name (run: a caller's `z` beside `d$z`, the forest took `d$z`); any other name is its value where
`forest()` was called. A value covers every row of the data and is cut by `subset` and the na.action, as
weights are.

**The matrix interface** (`dbarts(x, y, forests = )`). A basis is a value, full length, cut by `subset`. A name
that is a column of `x` is not found: `forest("x2", basis = factor(z))` with `z` a column of `x` fails with
"object 'z' not found where forest() was called". `bart()` has no forests argument; its forests come only as
formula terms.

**The data door** (`dbartsData(bases = )`).

    o <- dbartsData(y ~ x1 + x2, data = d, bases = list(NULL, factor(d$z)))

A list by forest position, of values only: a formula, a name and a string are refused. The list fixes the
number of forests, so `bases = list(NULL, NULL)` is a model of two forests with no basis. `dbartsData()`'s
formula refuses a `forest()` term. The object then records each forest's block (a factor's two indicator
columns) in `bases` and its levels in `basis.levels`; `basis.scale` is empty until a sampler is made.

**A data object handed to `dbarts()`.** With bases on it and a forests list that states only knobs, taken. A
list that declares a basis is refused: "a basis declaration cannot reach a pre-built data object through
dbarts()".

**`dbartsSpec()`.** There is no frame, so a basis is the caller's value at the call: with `z` only in `d`,
`forest(basis = factor(z))` fails. A value must have the object's kept-row count: on an object built with
`subset = 1:40` a full-length value is refused, where `dbarts()` refuses the cut one and takes the full one.
Beside an object carrying bases (none, z, w), a list of two declaring a basis gives two forests and a list of
three with a plain third drops w, with no message; a list stating only knobs keeps the object's bases.

**`$setForestBasis(forest, basis)`.** The forest by position or label. The basis a value, or a one-sided
formula read where it was written: `~ factor(z)` with `z` only in the data fails, as does `~ x2` with `x2` a
predictor. A string is refused ("a basis is a column, not its name").

**`predict(fit, newdata, bases = )`.** Reads code bases from `newdata` by name; `bases =` overrides by label or
position, and with no newdata replaces the fit's own at the training rows. A sampler's `$predict` refuses a
sampler of several forests.

**`$setData`** is refused on a sampler of several forests, so it is no door.

**Copy and reload.** `$copy()`, `saveRDS` and `readRDS`, `new("dbartsSampler", control, model, data)` and
`dbarts(s$data)` all carry the bases, the levels and the model's record of each code basis, and run.

**What is recorded.** The data object: `bases`, `basis.levels`, `basis.scale` (the column's sd over the rows
in, for a numeric forest stating no sd). The model: for each code basis its terms, its factor levels, its label
and the data columns it names. A `bart` fit: `bases` (the factor, for a forest of levels) and `basis.terms`.

**0.9-34** had none of this: no forests, no bases, no `dbartsSpec()`. Its `dbartsData()` took formula, data,
test, subset, weights, offset and offset.test.

## 2. The shape, door by door

The rule makes values data and references model: a forest (model) names columns; the data holds them; the
per-forest blocks on a data object are what a sampler wrote when it resolved the names.

### 2.1 What a pre-built data object keeps

The one real choice. A package that builds its data once and then calls `dbartsSpec()` needs the columns its
forests will name.

**A. The frame and the kept rows.**

    o <- dbartsData(y ~ x1 + x2, data = d)
    s <- dbarts(o, forests = list(forest(), forest(x2, basis = factor(z))))

The object keeps `data` as given, every row, and the numbers of the rows it kept. A forest beside it is read
exactly as `dbarts(y ~ x1 + x2, data = d, forests = )` reads it, `scale()` and `poly()` over every row
included, and would draw the same, as a data door and a list on the same rows draw the same today (run). A sampler drops the frame when it is made, so a sampler, a fit and a saved
sampler carry none. Cost: about 15 lines and two slots. In a session the frame costs nothing (run: a frame
stored in a slot is the same memory until either is changed); a saved data object writes it again (10.2 MB for
a frame, 20.3 MB for the frame and an object holding it). stan4bart saves its data object in its fit, and would
carry its BART columns once more. nlme keeps the whole frame in its fit the same way (`keep.data = TRUE`, its default).

**B. The columns named on the object.**

    o <- dbartsData(y ~ x1 + x2, data = d, columns = ~ z)

The object keeps only what it is told; the name `columns` is a placeholder. Cost: a new argument, the column
named twice, a refusal when a forest names one not kept, about 25 lines. The object again states something of
the model: which columns are bases.

**C. Forest terms in `dbartsData()`'s formula.**

    o <- dbartsData(y ~ x1 + x2 + forest(x2, basis = factor(z)), data = d)
    s <- dbartsSpec(o, forests = list(forest(), forest(sd = 0.5)))

The object resolves the bases and records the forests; a list beside it states knobs by label or position,
and a list declaring a basis is refused for good. Reads as `bart()`'s formula, and is how mgcv's `gam(fit =
FALSE)` and lme4's `lFormula()` build their objects, from the whole model formula. Cost: about 30 lines; the
data object states how many forests there are and which column each multiplies, the model dec-B663 sets
aside.

**D. The object's predictors only.**

    o <- dbartsData(y ~ x1 + x2 + z, data = d)
    s <- dbarts(o, forests = list(forest(c("x1", "x2")), forest("x2", basis = factor(z))))

A basis names a predictor column. It must then be a predictor, kept out of each forest's `vars` where it is not
to be split on; it exists on the kept rows only, so `scale()` reads other rows than at `dbarts()`; a factor
stored as categorical codes is rebuilt. About 20 lines. bartCause works as it stands, its treatment being a
predictor column already.

Recommended: A.

### 2.2 What `dbartsData()`'s formula accepts

Under A, B and D, as built: a `forest()` term is refused, its message pointing to a forests list. Under C, the
terms `bart()` takes.

### 2.3 How `dbartsSpec()` and `dbarts(o)` resolve

Under A: the names in a forest's basis are looked up in the kept frame, then where `forest()` was called, and
built on the kept rows, by the function `dbarts()` uses. `dbarts(o, forests = )` with a declaring list is
taken, not refused, so the two doors read one data object alike.

Rows: today `dbartsSpec()` wants a value cut to the object's rows and `dbarts()` refuses one. With the kept
rows recorded, a value is full length at every door and cut as weights are; a cut value is refused as
`dbarts()` refuses one. Cost: a call taken today, a kept-length value on an object that dropped rows, is
refused; where no row was dropped nothing changes. The alternative takes both lengths, two rules for one
argument.

### 2.4 `$setForestBasis` and `predict`

Kept as built: a forest is selected by its label or position, as every per-forest door selects it, and the new
basis is a value or a formula read where it was written. `predict(fit, newdata)` already resolves each code
basis by name in `newdata`, the reference itself; `bases =` stays for a basis not among the columns of
`newdata` and for the counterfactual at the training rows.

The alternative also selects a forest by the column its basis reads: `s$setForestBasis("z", v)` for the forest
of `factor(z)`, whose label is "factor(z)". Cost: a third key, ambiguous when two forests read one column
(`basis = dose` and `basis = I(dose^2)`), about 15 lines. Not recommended. Where the basis is a bare column the
label is already its name.

### 2.5 The matrix interface

No frame: a basis is a value, full length, cut by `subset` and the na.action, as weights are, at
`dbarts(x, y, forests = )` as built and at `dbartsSpec(dbartsData(x, y), forests = )` under 2.3.

    dbartsSpec(dbartsData(x, y, subset = s), forests = list(forest(), forest(basis = factor(z))))

The alternative lets a basis name a column of `x`, which then is a predictor of every forest that does not
leave it out by `vars`. Not recommended.

### 2.6 What a fit and a reload record

Nothing to add. The model's record of a code basis names the data columns it reads; a fit's `basis.terms`
rebuilds the basis at `newdata`; a copy or a reload is made from the model and the data's blocks, bit for bit.
Under A, `dbartsSpec()` and `dbarts(o)` write the same record `dbarts(y ~ ., data = d)` writes.

### 2.7 `dbartsData(bases = )` leaves

The argument goes from the exported function. `bases` and `basis.levels` stay as slots, `NULL` as created and
written by a sampler, as `basis.scale` and `cutPoints` are, so `dbarts(s$data)` keeps working, the record
reused whole (dec-B554). A call naming `bases =` gets R's "unused argument" error: it never shipped, so it gets
no message of its own. Inside, the fitting functions either call an internal builder that still takes the
blocks, or build the data object plain and resolve the forests afterwards through `dbartsSpec()`'s function,
one resolution for every door at about 40 lines more; that is the plan's to choose.

### 2.8 A record beside a declaration

After 2.7 only a sampler's own data carries blocks. Handed back with a list that declares a basis:

- (i) refused at both `dbarts()` and `dbartsSpec()`, naming both places; a list of knobs only keeps the
  record. About 8 lines. Recommended: re-creating from a sampler's data means its forests.
- (ii) the list replaces the record whole. A plain forest in it then has no basis, and the record has no frame
  to resolve a name against.
- (iii) any forests list beside a record refused, knobs too; `$setForest` changes knobs on a live sampler.

## 3. The landscape

Run where installed (R 4.6; mgcv 1.9.4, nlme 3.1.171, lme4 2.1.0, brms 2.23.0); stochtree and bcf are not
installed and were read from their source.

- `lm`, `glm`: a varying slope is a term, `y ~ x1 + x2:z`. Names are found in `data`, then the formula's
  environment; the model frame keeps every variable the formula names; `predict` reads them from `newdata` and
  fails without one ("object 'z' not found"). Weights are a name in `data` too, kept as `(weights)`.
  `lm.fit` and `glm.fit` take a matrix and vectors.
- mgcv: `s(x2, by = z)`, `z` a variable of `data`, kept in the model frame; `gam(fit = FALSE)` returns an
  object whose frame holds it, and `gam(G = )` fits from it. `by = "z"` is an error. `predict` without `z`
  warns that not all required variables were supplied.
- lme4: `(1 + z | g)`; `lFormula()`'s frame holds `z` and `g`; `predict` fails without `z`.
- nlme: `random = ~ z | g`; the fit keeps the whole data frame by default, a column the model never names
  included; `predict` fails without `z`.
- brms: `s(x2, by = z)`, `gp(x2, by = g)` and non-linear formulas (`b * z`, `b ~ x2`) need their variables in
  `data` or `data2`, and a caller's variable is refused.
- stochtree: `bart(X_train, y_train, leaf_basis_train = )`, `bcf(X_train, Z_train, y_train, ...)`;
  `predict(object, X, leaf_basis = )`, `predict(object, X, Z)`. Values only, no formula, no names.
- bcf (Hahn, Murray and Carvalho): `bcf(y, z, x_control, x_moderate, pihat, ...)`,
  `predict(object, ..., z_pred = )`. Values only.
- bartCause's `bcf()`: a formula and `data`, the treatment a column of `data` or a vector.
- BART and bartMachine: no forest multiplied by a basis.

The modeling packages with a formula name the column, and their pre-built objects hold the columns the model
names. The BART packages take values.

## 4. What a user meets

- **A formula user**, `bart(y ~ x1 + x2 + forest(x2, basis = factor(z)), data = d)`: nothing moves.
- **A forests-list user** with a frame, `dbarts(y ~ ., data = d, forests = )`: nothing moves.
- **A package that builds data once.**
  - stan4bart: `dbartsData(bartform, bartfr, factors = "categorical")`, then `dbartsSpec()` with one forest and
    no basis. Nothing moves; under A its saved fit carries `bartfr` once more.
  - treatSens: `dbartsData(x, y)`, the matrix interface, one forest. Nothing moves.
  - bartCause: `dbartsData(..., bases = list(NULL, factor(z)))`, then `dbarts(o, forests = )` stating knobs.
    It breaks when `bases =` goes. The edit gives the second forest `basis = ` the full-length factor it
    already builds and drops the bases line, about 6 lines. Its draws do not move (run: the data door, a list
    on `dbarts()` and a list on `dbartsSpec()` give identical draws, with and without `subset`). Its second
    forest's label becomes the basis's text instead of `forest2`; bartCause reads no label.
  - A caller of `dbartsSpec()` with a basis: a name that is a column of the frame now finds it. Where it is
    also a caller's variable, the column wins, as in `dbarts()`: the one call whose draws move. A value cut to
    the kept rows on an object that dropped rows is refused (2.3).
- **A matrix-interface user**, `dbarts(x, y, forests = )`: nothing moves.

## 5. The cost

- R, about 110 lines with no change under `src/` or to the shipped header: the two slots and their check (about
  10); `dbartsData()`, keeping the frame, recording the kept rows in each branch and dropping `bases` from the
  exported function (about 40); `dbartsSpec()`, resolving against the frame on the kept rows, full-length
  values and the refusal of 2.8 (about 30); `dbarts()`, taking a declaring list beside a data object and
  dropping the frame at creation (about 30).
- Tests: 83 calls of `dbartsData(bases = )` in 24 files. Read one by one: 45 are fixtures that build a sampler
  and become a forests list, drawing what they drew; 36 pin the data door itself (its refusals of a missing or
  infinite value, a short value, a single entry or all entries empty, its row cut, its ignored-argument
  warning, its forest count and labels, its place among the doors of the kind and defaults tests) and go where
  the forest door pins the same, or move there where it does not; 2 are oracles in the basis-terms tests (code
  against the same basis as a value) and are rewritten on the list's value. New: about 6 tests, among them
  `dbartsSpec(dbartsData(f, d), forests = L)` identical to `dbarts(f, d, forests = L)`, the frame gone from a
  sampler's data, full-length values, and the refusal of 2.8.
- Help: `dbartsData` (the `bases` item goes; the slots paragraph says the bases are a sampler's record and
  the frame is kept), `dbartsSpec` (where a basis is found), `dbarts` (the sentence on a declaration with
  nowhere to ride) and `dbartsSampler-class` (one mention). Four design notes (forest-kind-by-class,
  public-surface, written-surface, multiplier-combiner) name the data door; they change when this lands.
- NEWS: 0.9-34 had no bases, so nothing is retracted. The multi-forest item may gain a clause: a data object
  keeps its data, so a forest given to `dbartsSpec()` names its basis as in `dbarts()`.
- bartCause: about 6 lines and a check of its suite.

## 6. Recommendation

A, with 2.2 to 2.7 as written and (i) of 2.8. It is the rule read literally: the basis is part of the data
frame, the forest names it, and `dbartsSpec()` beside a data object reads a forest as `dbarts()` reads it beside
the frame, the sibling calls agreeing to the bit.

Flip to B if a data object should not carry its frame, saved objects and stan4bart's fit size outweighing the
second naming. Flip to C if a pre-built object should be the whole model, as mgcv's and lme4's are, and
`dbartsData()`'s formula should read as `bart()`'s. D only if the cost must be least.

The interim refusal (a data object carrying bases beside a list declaring one, refused at `dbartsSpec()`) is
still needed. Under A it is not subsumed: it becomes (i) of 2.8, the rule for a sampler's own data handed back,
and spreads to `dbarts()`, whose broader refusal narrows to the same case. Built now, it is not wasted; its
message loses the advice to use `dbartsData(bases = )` when this lands. Under C it is the permanent rule as
ordered.

## Calls made in this note

Each user-visible call is marked for the maintainer.

1. A over B, C and D (for the maintainer). Alternatives and costs in 2.1.
2. A value is full length at every door, a cut one refused at `dbartsSpec()` (for the maintainer). The
   alternative takes both lengths at `dbartsSpec()`.
3. `dbarts(o, forests = )` declaring a basis is taken under A (for the maintainer). The alternative leaves
   `dbarts()` refusing and only `dbartsSpec()` resolving, two readings of one data object.
4. A record beside a declaring list refused, (i) of 2.8 (for the maintainer). Alternatives (ii) and (iii).
5. A forest is selected by label or position only (for the maintainer). The alternative adds the column name
   as a key.
6. The matrix interface takes a value only (for the maintainer). The alternative resolves a name against the
   columns of `x`.
7. `bases =` gets R's own "unused argument" error (for the maintainer). The alternative is a message naming
   the forests list, as shipped arguments that retire get.
8. The frame is dropped when a sampler is made. The alternative keeps it on `s$data`, where nothing reads it.
9. The slots keep their names, `bases` and `basis.levels`. The alternative renames them as records, which
   moves every reader of them for no change of meaning.
10. How the fitting functions reach the blocks once the argument goes is left to the plan (2.7).
