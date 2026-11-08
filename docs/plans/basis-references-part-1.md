# basis-references, part 1: the data object and the forest doors

Status: PLANNED 2026-11-08, not built; awaiting its blind critique (slice 2 changes draws) and the maintainer's
answers to the open calls below. Design: [basis-references.md](../design/basis-references.md). Rests on
dec-B663, dec-B666, dec-B670, dec-B674 to dec-B677, dec-B679 to dec-B684 (dec-B682 revises dec-B674's words);
dec-B678 is part 2's (dec-B679).

agent: slices 1 and 2 serial, each one sonnet implementer; slice 1 a sonnet review (R surface), slice 2 an
opus review that runs its mutants (draws move for named calls); slice 3 one sonnet implementer per consumer
repository and one sonnet review.
rng: slice 1 NEUTRAL for every call. Slice 2 NEUTRAL for every call that runs today and keeps its meaning;
POSTERIOR-CHANGING for the calls listed under [Draws that move](#draws-that-move), whose new values are fixed
by the rulings and checked against the sibling door's draws. Expected bitwise against equivalence-38da480c
(55), bcf-equivalence-38da480c (19) and multinomial-equivalence-80b1c8d4 (11) after each slice; a moved
stream is a stop.
budget: about 365 lines of R (slice 1 about 265, slice 2 about 100 added and 60 removed), no `src/` or header
change; tests about 90 to 360 lines moved and 170 new in slice 1, about 250 rewritten and 80 new in slice 2;
six help pages, one NEWS entry; about 8 lines in bartCause and 1 in stan4bart.

## Goal

A data object built from a frame keeps the frame, the positions of the rows it kept and the coding record,
and reads its data through `model.frame(o)` and `nobs(o)`; `o@x`, `o@y`, `o@weights` and `o@offset` there
are stand-ins naming `model.frame(o)`. Every door that builds a model reads a forest's basis against the data
by one function: `dbartsSpec(dbartsData(f, d), forests = L)`, `dbarts(dbartsData(f, d), forests = L)` and
`dbarts(f, d, forests = L)` draw the same. `dbartsData(bases = )` is gone. A sampler's data is today's
block-holding object, with the caller's frame dropped and its column names kept; part 2 makes it the model
frame.

## Context

- The design, sections 2.1 to 2.3 and 2.7 to 2.9 of [basis-references.md](../design/basis-references.md),
  and its [Calls made](../design/basis-references.md#calls-made). Section 2.4's model frame, 2.5's block-holding
  sampler frame and 2.6's sizes are part 2's.
- dec-B679 splits the design; the 305-line sizing it names is the note's fourth revision (D1: the frame
  dropped at the resolution, today's object handed on). This plan is that design with dec-B677's accessors
  in place of `extract()` and dec-B680 to dec-B684 in place of calls 1, 2, 3, 6 and 7.
- Today a basis beside a frame rides [`dbartsData`](../../R/data.R)'s own `bases` argument: a value as itself,
  code as the row numbers of [`basisRowNumbers`](../../R/forestBasis.R), cut there by
  [`validateForestBases`](../../R/data.R) with [`alignSubsetRowsToFrame`](../../R/data.R), and built on the
  kept rows by [`buildFitBases`](../../R/forestBasis.R). [`dbartsSpec`](../../R/spec.R) reads a basis with no
  data ([`readDeclaredBases`](../../R/forestBasis.R) given `NULL`) at the object's kept length.
- Premises, run 2026-11-08 on main at 7a54b2b1d (shipped build) and on 0.9-34:
  - 0.9-34: `dbartsData(y ~ x1 + f, d, subset = 2:10, weights = w)@x` is the 9 by 4 model matrix (`x1 f.a f.b
    f.c`), `o@y` and `o@weights` numeric of 9; `model.frame(o)` and `nobs(o)` fail ("$ operator not defined
    for this S4 class", "no 'nobs' method is available"). Main: `o@x` a `dbartsMixedMatrix` of 9 by 2, the rest
    as 0.9-34.
  - `dbartsSpec(dbartsData(y ~ x1 + x2, d), forests = list(forest(), forest(basis = factor(z))))` reads the
    workspace `z` where `dbarts(y ~ x1 + x2, d, forests = <the same>)` reads `d$z`.
  - Beside a two-forest sampler's data, `dbartsSpec(s2$data)`, a list of tree counts alone and
    `list(forest())` each give two forests keeping the held basis.
  - `dbarts(o, forests = <a declaring list>)` stops ("'forests' declares a 'basis' but 'formula' is already a
    dbartsData"); a 60-row vector beside `dbartsData(..., subset = 1:40)` is refused, a 40-row one taken.
  - `dbartsSpec(dbartsData(y ~ x1 + x2, data = e), forests = <basis w>)` reads the workspace `w`, `e` holding
    another `w`.
  - A probit fit's 0/1 weights ride `s$data@weights` and survive `s$setData(s$data)`, so the help sentence
    "Under such a family the sampler's own data object carries no weights, so setData(sampler$data) drops the
    mask in force" in [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd) is false on main (root TODO,
    sampler-help-mask-sentence).

## Constraints

- No `src/` change, no change to `inst/include/dbarts/dbarts.h`, no engine change.
- Every formula fit, `dbarts(f, d)`, `bart(f, d)` and the matrix interface keeps its draws bitwise: they build
  the block as today and never go through the block maker.
- A sampler's data, `dbartsSpec()`'s returned `data` and every saved sampler hold no caller frame.
- Out of scope, part 2: the sampler's model frame, its setters, the stand-ins on `s$data` (call 36), what
  `model.frame(s$data)` decodes (call 39, but see open call 2), the mask in the frame (call 41), stan4bart's
  bundle and its three reads of `sampler$data`, the lint rule (call 37), the bridge's whole-matrix split
  (call 31). Also out: call 21 (a fit's own basis vector), which is keepfits-forest-records' (dec-B671); the
  examination's other recommendations for bartCause (root TODO, bartcause-data-route) and stan4bart's
  hand-written slots.
- A data object saved by a development build before this lands is owed nothing (call 22); one saved by
  0.9-34 already fails on main.

## The doors

| door | part 1 |
|---|---|
| `dbartsData(f, d)`, `d` a data frame or list | keeps `d` whole, `rows`, the coding record, the frozen response, weights and offset expressions, `frame.columns`; `x`, `y`, `weights`, `offset` stand-ins (slice 1). `bases =` gone (slice 2) |
| `dbartsData(f, d, keep.data = FALSE)` | the block, `y`, `weights`, `offset` as today, plus `rows` and `frame.columns` (slice 1) |
| `dbartsData(f, data = e)`, an environment | the block as today, plus `rows` and the environment's names in `frame.columns` (dec-B683; the part-1 shape of "the shape a sampler's data holds") |
| `dbartsData(f)` with no data, `dbartsData(x, y)` | the block as today, plus `rows`; no `frame.columns` |
| the resolution, at `dbarts(o)`, `dbartsSpec(o)`, `bart(o)`, `rbart_vi(o)`, `xbart(o)`, [`ingestHurdleData`](../../R/bart.R), `new("dbartsSampler", data = o)`, [`dbartsSampler$setData`](../../R/dbarts.R) | a frame-keeping object becomes the block-holding object, made from the frame in the build's order, `frame` and `rows` dropped, `frame.columns` kept; any other object passes as it is (slice 1) |
| `dbarts(f, d, forests = )`, `bart(f, d)` | builds with `keep.data = FALSE` internally; bases read by the shared reader against `d`, cut by `rows` (slice 2); draws bitwise |
| `dbarts(x, y, forests = )` | the shared reader, a value full length cut by `rows` (slice 2); draws bitwise |
| `dbarts(o, forests = )`, `dbartsSpec(o, forests = )` | a declared basis read against `o` by the shared reader (dec-B681), full length (dec-B680), names per dec-B682 and dec-B683, the list stating the forests in full (dec-B684) (slice 2) |
| `model.frame(o)`, `nobs(o)` | new methods on every data object (dec-B677): from the frame where one is kept, from the block otherwise (slice 1; open call 2) |
| `o@x`, `o@y`, `o@weights`, `o@offset` | stand-ins on a frame-keeping object only (dec-B675, dec-B676); real slots everywhere else, `s$data` included (slice 1) |
| `makeTestModelMatrix(o, newdata)` | reads the coding record on a frame-keeping object (slice 1) |
| `$setForestBasis`, `predict(newdata)`, `predict(bases = )` | unchanged |

## Slice 1: the data object keeps its frame

Tier: R surface, no draw change ([Process by risk](README.md#process-by-risk)). Files: `R/A_class.R` (slots,
validity), `R/data.R` ([`dbartsData`](../../R/data.R), the block maker, the stand-ins,
[`makeTestModelMatrix`](../../R/data.R), `model.frame` and `nobs` methods), the door entries in `R/dbarts.R`,
`R/spec.R`, `R/bart.R`, `R/rbart.R`, `R/xbart.R`, `R/bartcore.R`, `NAMESPACE`, `man/dbartsData.Rd`, a new page
for the two methods, `man/makeind.Rd`, `man/dbartsSampler-class.Rd` (`setData` takes a frame-keeping object),
`_pkgdown.yml`, `inst/NEWS.Rd`.

### Steps

1. Slots: `frame`, `rows` (integer positions, `NULL` when every row is kept in order), `frame.columns`, the
   coding record (the block with no rows, with a missing-value flag per column) and the frozen response,
   weights and offset expressions; `x`, `y`, `weights`, `offset` widen to hold a stand-in. Validity per kind:
   a frame-keeping object's positions lie within the frame's rows; a block-holding one is checked as today.
   About 50 lines with the stand-ins.
2. The stand-in: one small classed object per slot; `dim`, `dimnames`, `length`, both subset operators, `as.matrix`,
   `as.vector`, `mtfrm`, `Ops`, `Math` and `Summary` stop with "the data are read with model.frame(o)";
   `print` and `format` show one line, so `show(o)` and `str(o)` work; `attr()` gives `NULL` (dec-B676).
3. [`dbartsData`](../../R/data.R) takes `keep.data = TRUE`. It builds and checks everything as today, every
   refusal staying at the build (dec-B675), then on a data frame or list stores the frame, the positions
   (the subset value it already evaluates once, less the na.action's rows), the coding record and the frozen
   expressions (a column stays a name, a caller's value its value, dec-B185), and drops the block and
   values. Row names only where the frame cannot give them (call 19). An environment keeps the block and
   `ls(e)`. About 70 lines.
4. The block maker and the resolution: from a frame-keeping object, evaluate the frozen terms on the whole
   frame, cut to the positions, code against the record's levels, lift sparse columns as the build does, and
   remake `y`, `weights` and `offset` the same way (dec-B675); the result is today's object, `frame` and
   `rows` dropped. One function, called first at each door in the table; reads of a data object's slots
   before that point (in `bart()`, `dbarts()`, the hurdle's entry) move after it. Internal builders
   (`dbarts()`, `xbart()`, `rbart_vi()`, the hurdle) pass `keep.data = FALSE`. About 100 lines.
5. `model.frame()` and `nobs()` methods, registered in `NAMESPACE` (dec-B677). On a frame-keeping object,
   the maker's variables over the rows in, factors as factors, `(weights)` and the `offset()` and `(offset)`
   columns, as `model.frame()` of an `lm` fit; on a block-holding object, per open call 2. `nobs()` counts
   every row in, weight-0 rows included; the help says how this differs from `lm`. About 45 lines.
6. [`makeTestModelMatrix`](../../R/data.R) reads the coding record and its flags. About 5 lines.
7. Tests and help; the reads of the old slots in tests move to `model.frame()`, to the made object or to
   `keep.data = FALSE`. A static scan finds 83 lines in 17 files reading `@x`, `@y`, `@weights` or `@offset` of a
   named formula-built object; the note's bound by name is about 360.

### Pins (fail before, pass after)

- The sibling identity, for each case of the note's section 2.3 list: `dbarts(dbartsData(f, d, ...))` and
  `dbarts(f, d, ...)` draw identically, and the resolution's block, `y`, `weights` and `offset` are
  `identical()` to those of `dbartsData(f, d, ..., keep.data = FALSE)`.
- The positions under each `subset` and na.action of the note's section 2.1, a repeated row and a character
  `subset` among them.
- The stand-ins: each method above stops naming `model.frame(o)`; `attr(o@x, "factor.levels")` is `NULL`;
  `print(o)` and `str(o)` run.
- `model.frame(o)`: `model.response()`, `model.weights()` and `model.offset()` equal those of the `lm` fit on
  the same call; `nobs(o)` keeps rows of weight 0.
- `d$z[1] <- 9` after the build leaves `o@frame$z` and the fit unchanged; a later change of `e$x1` leaves an
  environment's object unchanged.
- `keep.data = FALSE` and the matrix interface: every slot reads as today; `frame.columns` holds `names(d)`.
- Every door in the table takes a frame-keeping object, and `s$data` and `dbartsSpec(o)$data` hold no frame.
- `makeTestModelMatrix(o, newdata)` on a default object equals it on `keep.data = FALSE`.
- Saved size: a default object against `keep.data = FALSE`, and a sampler within 0.001 MB of main's.

Mutants for the review: the maker evaluating terms on the kept rows (predict's path; the caller's vector
under `subset` fails); the na.action cut skipped (`na.omit` fails); one door left without the resolution
(`$setData(o)` fails).

### Gates

Lint (`lintr::lint_package()` against the slice's library), air, codoc, win-drift, doc-freshness; the touched
test files; the four snapshot files on the reference build (none builds a data object with a forests list, so
none should move); the equivalence trio bitwise on the reference build; the exact gates with `quick` (a fit
now carries new data slots); `pkgdown::check_pkgdown(".")`; NEWS parse; `R CMD check --as-cran` from a clean
tarball. CI runs the full suite, check and sanitizers.

## Slice 2: the forest doors

Tier: changes draws, for the named calls only; no engine change. Files: `R/data.R` (`bases` out of
[`dbartsData`](../../R/data.R), [`validateForestBases`](../../R/data.R)'s cut moved),
`R/forestBasis.R` (the shared reader beside [`readForestBasis`](../../R/forestBasis.R);
[`basisRowNumbers`](../../R/forestBasis.R) and [`fitBasisRows`](../../R/forestBasis.R) retired), `R/model.R`
([`resolveFormulaBasisSubset`](../../R/model.R)), `R/dbarts.R`, `R/spec.R`, `R/bart.R`, `man/dbartsData.Rd`,
`man/dbartsSpec.Rd`, `man/dbarts.Rd`, `man/forest.Rd` where it names the data door,
`man/dbartsSampler-class.Rd`, `inst/NEWS.Rd`.

### Steps

1. One reader, `resolveForestBases(forests, data, rows, held)`: reads each declared basis against the data's
   frame (a frame-keeping object's, or the `d` `dbarts(f, d)` holds), then the held columns (open call 1),
   then where `forest()` was called; a value has one entry per row handed over and is cut by `rows`, code is
   built on the kept rows (dec-B680). `dbarts(f, d)`, `dbarts(x, y)`, `dbarts(o)` and `dbartsSpec(o)` call
   it; the row-number plumbing through `bases` and the restating of "'bases'" as "'basis'" go. About 30 lines
   added, 60 removed.
2. `dbarts(o, forests = )` takes a declaring list, read exactly as `dbartsSpec(o, forests = )` reads it
   (dec-B681). About 25 lines.
3. Beside an object with no caller frame (`s$data`, `keep.data = FALSE`), a basis naming a column the object
   holds (a block column, or a held basis, open call 1) is read from it; a name in `frame.columns` that is no
   such column is refused, naming the data object built with its frame; any other name is the caller's value
   (dec-B682). Beside an environment's object, the same, the refusal naming a data frame or list (dec-B683).
   About 15 lines.
4. Beside any data, `s$data` included, a forests list states the forests in full: a forest declaring no
   basis has none, a list of tree counts alone gives plain forests, no list gives one forest (dec-B684). A
   declared basis takes its scale and levels from the column it reads, as a declaration does today. dec-B663's
   interim refusal is not built. About 10 lines.
5. `bases =` leaves [`dbartsData`](../../R/data.R): R's own "unused argument" (call 16). The three texts that
   advise the data door are reworded (the `$setData` refusal, `dbarts()`'s refusal, which goes, and the
   refusal of a one-forest model with a basis). About 20 lines.
6. Tests: the 83 calls of `dbartsData(bases = )` in 24 files become forests lists (45 fixtures that draw what
   they drew, 36 pins moved to the forest door, 2 oracles on the list's value). A forest built from a declared
   basis is labeled by the basis's text where `bases =` gave `forest2`; tests reading labels change with it.

### Draws that move

- D1. Beside a frame-keeping object, a basis name that is both a column of the frame and a caller's
  variable: the column, as at `dbarts(f, d)`; today the caller's.
- D2. Beside `s$data` or a `keep.data = FALSE` object, a name that is a held column and a caller's variable:
  the held column (dec-B682); today the caller's.
- D3. Beside a sampler's data of several forests: no list gives one forest and a list of tree counts plain
  forests (dec-B684); today the held bases are kept.

No test file and no snapshot moves. A parse of every test file finds 62 door calls with a forests list beside
a data object: the frame-keeping ones name no column of their frame, the one sampler's data of several
forests (`original$data` in test-forest-scale-record.R) takes a list stating both forests, and the one name
collision, ["plainData <- dbartsData(y ~ x1 + x2 + x3, frame)"](../../inst/tinytest/test-control-doors.R)
beside `factor(z)`, has the workspace `z` equal to `frame$z` and is refused for its tree count. The doors
handed a sampler's data with no list (16, in test-control-doors.R, test-standardization-record.R and
test-leaf-standardization-rows.R) are all of one forest. The four snapshot files build no data object with a
forests list. The implementer confirms each by running the touched files against main's library and the
slice's.

### Pins (fail before, pass after)

- D1, D2 and D3 each, with the caller's variable unlike the column, the new draws `identical()` to the
  sibling door's: `dbarts(f, d, forests = L)` for D1, a fit naming the held values for D2, and
  `dbarts(f, d, forests = <the list as stated>)` for D3.
- dec-B680 at `dbarts(o)`, `dbartsSpec(o)` and a `keep.data = FALSE` object under `subset = 1:40`: a 60-row
  vector taken, equal to `dbarts(f, d, subset = 1:40, forests = )`; a 40-row one refused with the length
  expected.
- dec-B681: `dbarts(o, forests = <declaring>)` taken, identical to the `dbartsSpec()` route and to
  `dbarts(f, d, forests = )`.
- dec-B682 and dec-B683: each refusal's text, and the caller's value for a name of neither kind.
- `dbartsData(bases = )` gives "unused argument".

Mutants for the review: the reader skipping the frame beside a frame-keeping object (D1 fails); a held basis
kept beside a list of tree counts (D3 fails); a value cut by the kept length (dec-B680 fails).

### Gates

Slice 1's, with the equivalence trio bitwise the gate for D's absence elsewhere (bart2twoforest and the 19
bcf scenarios read their bases through the new reader), and the full suite once, locally, in one process.

## Slice 3: the consumers

Tier: R surface, in each consumer's repository, prepared against slice 2's library and pushed to the
consumer's main right after dbarts' merge (the patches name an argument main does not have before it).
Under dec-B674's examination clause, from the examination recorded in the note's section 2.5:

- bartCause passes `keep.data = FALSE` at its two data call sites (three objects): it reads the block
  (`nrow()`, `colnames()`, `x[, z]`, `as.matrix(x)`), `@y` and `@weights`, writes `@x.test`, and its saved
  fits keep `data.rsp`, which would otherwise carry the user's whole frame. `bcf()` drops `bcf.bases` and
  gives the treatment forest `basis =` the full-length `factor(z)` it builds, cut by the object's `rows` under
  `subset`. Its reads of `attr(fit$fit$data@x, "term.labels")` off kept samplers stay (part 1's `s$data` holds
  the block). About 8 lines. Pins: `bartc()` and `bcf()` fits identical before and after, with and without
  `subset`; its suite.
- stan4bart passes `keep.data = FALSE` in `glFormula`: its fit's `bartData` reads (`@x`, `@varTypes`,
  `factor.levels`, `makeTestModelMatrix()`), its slot writes and its bundle stay as today. Its reads of
  `sampler$data@y`, `@weights` and `data@offset` stay in part 1. 1 line. Pin: its suite, a fit identical.
- treatSens uses the matrix interface and reads no slot: nothing.

## Records

At the landing: ledger entries for the calls made here, this plan's Landing note and Status, the design's
Status line, TODO, the plans and design indexes; the four design notes that name the data door
(forest-kind-by-class, public-surface, written-surface, multiplier-combiner). NEWS, one entry against
0.9-34: "A data object made by dbartsData() from a formula and a data frame keeps the frame and the
positions of the rows it keeps; its x, y, weights and offset slots, which held the model matrix and the
response in 0.9-34, no longer hold them: model.frame(object) returns them and nobs(object) counts the rows.
A forest given to dbarts() or dbartsSpec() beside a data object finds its basis in the frame.
keep.data = FALSE keeps the matrix and the response in the object instead, as before." Part 2 extends it.

## Verification

- `R_LIBS=<lib> Rscript -e 'tinytest::run_test_file(<each touched file>)'` with `library(dbarts)` attached:
  0 failed; the new pins fail against main's library.
- `Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-38da480c.rds --bitwise
  --strict-coverage`, and the bcf and multinomial harnesses against their current files, on a `--preclean`
  reference install: 55, 19 and 11 "identical draws" lines, no "max |z|" line.
- The four snapshot files on the reference build: 27 results, 0 failed.
- The exact gates of `.github/workflows/exact-gates.yaml` with `quick`: all pass.
- Each gate of [Process by risk](README.md#process-by-risk) on its own exit status.

## What part 2 reshapes

The resolution's product (today's block-holding object becomes the model frame, about 30 lines of slice 1's
100); `model.frame()` on a sampler's data (call 39) and on the other block-holding objects (open call 2,
about 20 lines); the held basis's name (open call 1, beside call 28's `(basis.<label>)` columns, about 10);
an environment's object (dec-B683's model frame, about 5); the stand-ins and their attributes on `s$data`
(call 36); stan4bart's three reads and its bundle; bartCause's `term.labels` reads, under call 36. Between
the two landings main holds both shapes (dec-B679).

## Open calls

Flagged for the maintainer; none is decided here.

1. **How a held basis is named beside a sampler's data.** dec-B682 reads "a basis column" the object holds
   and dec-B684 keeps a held basis by naming it, `forest(basis = z)`. Today `s2$data@bases` is by forest
   position with the coded columns (`factor(z)0`, `factor(z)1`); the name `z` is in the model's forest
   record, not the data. Shown: (a) by the basis's written text, recorded with the held basis
   (`forest(basis = z)` takes a basis written `z`, `forest(basis = factor(z))` one written `factor(z)`; a
   source column named otherwise is refused under dec-B682 since it was a column of the caller's frame),
   about 10 lines, recommended, as it is the spelling dec-B684 names; (b) by its source columns, `s$data`
   keeping each column a code basis read (`z`), so any expression of them is read again, one column of n
   doubles each; (c) by the forest's label, as part 2's call 28 names its columns, so `forest(basis = z)` on a
   forest labeled `treat` finds nothing.
2. **What `model.frame()` of a block-holding object shows in part 1** (`s$data`, `keep.data = FALSE`, the
   matrix interface, an environment's object). dec-B677 has it read "from the stored table otherwise"; this is
   call 39's question, reached a landing early. Shown: (a) decoded where exact (the response to its recorded
   type and levels, a factor column as a factor, a matrix term's columns joined back), the rest as held, call
   39's recommendation, about 20 lines; (b) as held. Recommended: rule call 39 before slice 1 is built and
   build its answer for both kinds; failing that, (a), which part 2 keeps.
3. **Objects with no frame that dec-B682 does not name** (the matrix interface's, and a formula's with no
   `data`). Recommended: a basis name read as their building door reads it, where `forest()` was called,
   so `dbartsSpec(dbartsData(x, y), forests = L)` draws as `dbarts(x, y, forests = L)`; their samplers' data
   follow dec-B682. Alternative: dec-B682's held-column rule for every frameless object, under which the two
   differ where a basis names a column of `x` that is also a workspace variable.
4. **A formula term that draws** (`I(x + rnorm(n))`) is evaluated again at the resolution, so
   `dbarts(dbartsData(f, d))` fits other values than `dbarts(f, d)` and consumes the generator. Recommended:
   a help sentence, as `predict.lm` re-evaluates such a term on newdata. Alternative: keep each term's
   evaluated columns on the object, the block's size again, against dec-B670.
5. **The note's agent calls this part builds on**, unmarked until the session's close-out (dec-B629): 12
   (positions, nothing stored when every row is kept in order), 13, 14 (a forest's record names columns), 15,
   16 (`bases =` gets "unused argument"), 17, 18 (`x.test` kept as built), 19 (row names only where the frame
   cannot give them), 20 (`model.frame(o)` in the build's order) and 22.

## Calls made in planning

Agent calls, each with its alternative, for the close-out:

- `keep.data = FALSE`, an environment's and the matrix interface's objects keep `rows`, since dec-B680's
  length rule at every building door needs it (bartCause's `bcf()` under `subset`). Alternative: none that
  keeps dec-B680.
- `dbarts(f, d)` builds with `keep.data = FALSE` internally and never runs the block maker, so no formula
  fit depends on it. Alternative: `dbarts(f, d)` as `dbarts(dbartsData(f, d))`, one path, every formula fit
  then remaking its block (0.004 s at 20,000 rows) and depending on the maker.
- A declared basis beside a sampler's data takes its scale from the column it reads, as a declaration does
  today, not the held `basis.scale` by position, since a list stating the forests in full may reorder them.
  Alternative: the held scale carried to the forest naming it.
- Slice 2 takes the changes-draws tier for its three named moves (a slice spanning two tiers takes the
  higher). Alternative: the R-surface tier, the sibling-door oracle standing in for the opus review.
- The consumer patches land after dbarts' merge, pushed at once. Alternative: patches that test
  `"keep.data" %in% names(formals(dbarts::dbartsData))`, mergeable first, a branch on a development build.
