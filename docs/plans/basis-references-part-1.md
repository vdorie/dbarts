# basis-references, part 1: the data object and the forest doors

Status: PLANNED 2026-11-08, revised after its blind critique; not built. Awaiting the maintainer's answers to
the four open calls below and a word on the size (about 450 lines of R against dec-B679's 305). Design:
[basis-references.md](../design/basis-references.md). Rests on dec-B663, dec-B666, dec-B670, dec-B674 to
dec-B677, dec-B679 to dec-B686 (dec-B682 revises dec-B674's words; dec-B685 and dec-B686 also shape this
part); dec-B678, dec-B687, dec-B688 and dec-B692 are part 2's (dec-B679).

agent: slices 1 and 2 serial, each one sonnet implementer; slice 1 a sonnet review (R surface) that runs its
mutants, slice 2 an opus review that runs its mutants (draws move for named calls); slice 3 one sonnet
implementer per consumer repository and one sonnet review.
rng: slice 1 NEUTRAL for every call. Slice 2 NEUTRAL for every call that runs today and keeps its meaning;
POSTERIOR-CHANGING for the calls listed under [Draws that move](#draws-that-move), whose new values are fixed
by the rulings and checked against the sibling door's draws. Expected bitwise against equivalence-38da480c
(55), bcf-equivalence-38da480c (19) and multinomial-equivalence-80b1c8d4 (11) after each slice; a moved
stream is a stop.
budget: about 450 lines of R added (slice 1 about 315, slice 2 about 135) and about 270 removed (slice 2), no
`src/` or header change; tests about 90 to 360 lines moved and 250 new in slice 1, about 250 rewritten and
120 new in slice 2; seven help pages, one NEWS entry; about 8 lines in bartCause and 1 in stan4bart.

## Goal

A data object built from a frame keeps the frame, the positions of the rows it kept and the coding record,
and reads its data through `model.frame(o)` and `nobs(o)`; `o@x`, `o@y` and, where given, `o@weights` and
`o@offset` there are stand-ins naming `model.frame(o)`. Every door that builds a model reads a forest's basis
against the data by one function: `dbartsSpec(dbartsData(f, d), forests = L)`, `dbarts(dbartsData(f, d),
forests = L)`, `dbarts(f, d, forests = L)` and a formula's `forest()` term draw the same. `dbartsData(bases = )`
is gone. A sampler's data is today's block-holding object, with the caller's frame dropped and its column
names kept; part 2 makes it the model frame.

## Context

- The design, sections 2.1 to 2.3 and 2.7 to 2.9 of [basis-references.md](../design/basis-references.md),
  and its [Calls made](../design/basis-references.md#calls-made). Section 2.4's model frame, 2.5's block-holding
  sampler frame and 2.6's sizes are part 2's.
- dec-B679 splits the design; the 305-line sizing it names is the note's fourth revision. This plan is that
  design with dec-B677's accessors, dec-B680 to dec-B684 in place of calls 1, 2, 3, 6 and 7, and the doors the
  fourth revision did not count (the formula term door, the response kinds, hazard, `model.frame()` of every
  kind); hence about 450 lines.
- Rulings that reach this part beyond its own: dec-B685 and dec-B677 settle what `model.frame()` of a
  block-holding object shows (decoded where exact), built here for every kind; dec-B686 keeps the coding
  record's attributes on `s$data@x`, so the reads of `attr(s$data@x, "term.labels")` stay in both parts;
  dec-B694 (a hazard sampler's `$setData` takes survival data, root TODO hazard-setdata-surv) needs the
  resolution to carry the Surv record, which it does; dec-B700 (a single string names its column, root TODO
  basis-string-name) resolves as dec-B682 does, so whichever lands second puts the string read in the shared
  reader.
- Today a basis beside a frame rides [`dbartsData`](../../R/data.R)'s own `bases` argument: a value as itself,
  code as the row numbers of [`basisRowNumbers`](../../R/forestBasis.R), cut there by
  [`validateForestBases`](../../R/data.R) with [`alignSubsetRowsToFrame`](../../R/data.R), and built on the
  kept rows by [`buildFitBases`](../../R/forestBasis.R). A formula's `forest()` term takes the same route:
  [`ingestFormulaTerms`](../../R/formulaTerms.R) reads its rows with [`fitBasisRows`](../../R/forestBasis.R)
  and hands them on as `basisRows`. [`dbartsSpec`](../../R/spec.R) reads a basis with no data
  ([`readDeclaredBases`](../../R/forestBasis.R) given `NULL`) at the object's kept length.
- Premises, run 2026-11-08 on main at 7a54b2b1d (shipped build) and on 0.9-34:
  - 0.9-34: `dbartsData(y ~ x1 + f, d, subset = 2:10, weights = w)@x` is the 9 by 4 model matrix (`x1 f.a f.b
    f.c`), `o@y` and `o@weights` numeric of 9; `model.frame(o)` and `nobs(o)` fail ("$ operator not defined
    for this S4 class", "no 'nobs' method is available"). Main: `o@x` a `dbartsMixedMatrix` of 9 by 2, the rest
    as 0.9-34. 0.9-34 had no families, no forests and no bases.
  - `dbartsSpec(dbartsData(y ~ x1 + x2, d), forests = list(forest(), forest(basis = factor(z))))` reads the
    workspace `z` where `dbarts(y ~ x1 + x2, d, forests = <the same>)` reads `d$z`. Beside a matrix object,
    a formula object with no `data` and the matrix sampler's data, a basis name reads the workspace.
  - Beside a two-forest sampler's data, `dbartsSpec(s2$data)`, a list of tree counts alone and
    `list(forest())` each give two forests keeping the held basis. A one-forest creation on data that carry
    bases is refused ("the data carry forest bases but no basis forest was configured").
  - After `s2$setForestBasis(dose / 10, forest = 2)` the held scale record stays sd(dose), 0.605;
    `dbartsSpec(s2$data)` carries it (the new forest's coefficient sd 2.52), while a list declaring the
    column in force takes 25.2.
  - `dbarts(o, forests = <a declaring list>)` stops ("'forests' declares a 'basis' but 'formula' is already a
    dbartsData"); a 60-row vector beside `dbartsData(..., subset = 1:40)` is refused, a 40-row one taken; the
    same beside `dbartsData(X, y, subset = 1:40)`, while `dbarts(X, y, subset = 1:40, forests = )` takes 60.
  - Hazard, 60 subjects and 1830 period rows: the formula route refuses `forest(basis = dose)` ("60 rows
    where the fit has 1830 observations"); the matrix route refuses 60 values and takes 1830. `aft` takes no
    second forest.
  - `dbarts(dbartsData(Surv(t, s) ~ x1 + x2, d), family = "aft")` and `family = "hazard"` build a sampler,
    reading the object's `survivalStatus`, `survivalTime` and `survivalTimeOmitted`.
  - A zero-weight build warns once at `dbartsData()`; `dbarts(o)` and `dbartsSpec(o)` are silent. An object
    with no weights or offset has `NULL` in those slots. `dbartsData(o, keep.data = FALSE)` is "unused
    argument".
  - `dbartsSpec(dbartsData(y ~ x1 + x2, data = e), forests = <basis w>)` reads the workspace `w`, `e` holding
    another `w`.
  - A probit fit's 0/1 weights ride `s$data@weights` and survive `s$setData(s$data)`, so the help sentence
    "Under such a family the sampler's own data object carries no weights, so setData(sampler$data) drops the
    mask in force" in [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd) is false on main (root TODO,
    sampler-help-mask-sentence).

## Constraints

- No `src/` change, no change to `inst/include/dbarts/dbarts.h`, no engine change.
- Every formula fit, `dbarts(f, d)`, `bart(f, d)` and the matrix interface keeps its draws bitwise: they build
  the block as today and never go through the resolution.
- A sampler's data, `dbartsSpec()`'s returned `data` and every saved sampler hold no caller frame.
- Out of scope, part 2: the sampler's model frame, its setters, the stand-ins on `s$data` (dec-B686), the
  mask in the frame (dec-B688), stan4bart's bundle (dec-B687) and its three reads of `sampler$data`,
  `extract(s, "predictors")` (dec-B692), the lint rule (call 37), the bridge's whole-matrix split (call 31).
  Also out: call 21 (a fit's own basis vector), which is keepfits-forest-records' (dec-B671); the
  examination's other recommendations for bartCause (root TODO, bartcause-data-route); stan4bart's
  hand-written slots (root TODO, consumer-slot-writes); `rbart_vi(o, group.by = )` reading the caller's
  variable (root TODO, rbart-group-by-frame); dec-B694's survival `$setData` and dec-B700's string names.
- A data object saved by a development build before this lands is owed nothing (call 22); one saved by
  0.9-34 already fails on main.

## The doors

| door | part 1 |
|---|---|
| `dbartsData(f, d)`, `d` a data frame or list | keeps `d` whole, `rows`, the coding record, the response's record (its coding, and a Surv response's status, time and omitted rows), the frozen response, weights and offset expressions, `frame.columns`; `x`, `y` stand-ins, `weights` and `offset` stand-ins where given and `NULL` otherwise (slice 1). `bases =` gone (slice 2) |
| `dbartsData(f, d, keep.data = FALSE)` | the block, `y`, `weights`, `offset` as today, plus `rows` and `frame.columns` (slice 1) |
| `dbartsData(f, data = e)`, an environment | the block as today, plus `rows` and `ls(e, all.names = TRUE)` in `frame.columns` (dec-B683; the part-1 shape of "the shape a sampler's data holds") |
| `dbartsData(f)` with no data, `dbartsData(x, y)` | the block as today, plus `rows`; no `frame.columns` |
| `dbartsData(o, ...)`, a data object | returned as it is; `keep.data` joins the arguments it takes quietly |
| the resolution, at `dbarts(o)`, `dbartsSpec(o)`, `bart(o)`, `rbart_vi(o)`, `xbart(o)`, [`ingestHurdleData`](../../R/bart.R), `new("dbartsSampler", data = o)`, [`dbartsSampler$setData`](../../R/dbarts.R) | a frame-keeping object becomes the block-holding object `dbartsData(f, d, keep.data = FALSE)` would have made, `frame.columns` kept; any other object passes as it is; no warning the build gave is given again (slice 1) |
| `dbarts(f, d, forests = )`, `bart(f, d)` | builds with `keep.data = FALSE` internally; bases read by the shared reader against `d`, cut by `rows` (slice 2); draws bitwise |
| a formula's `forest()` term, `bart(y ~ forest(x2, basis = z), d)` and `dbarts()` | [`ingestFormulaTerms`](../../R/formulaTerms.R) reads its bases by the shared reader; `fitBasisRows` and the `basisRows` hand-on go (slice 2); draws bitwise |
| `dbarts(x, y, forests = )` | the shared reader, a value full length cut by `rows` (slice 2); draws bitwise |
| `dbarts(o, forests = )`, `dbartsSpec(o, forests = )` | a declared basis read against `o` by the shared reader (dec-B681), full length (dec-B680), names per dec-B682, dec-B683 and open call 2, the list stating the forests in full (dec-B684), a held basis's scale per open call 1 (slice 2) |
| a hazard fit's basis, at `dbarts(f, d)`, `dbarts(x, y)`, `dbarts(o)` | its length per open call 4 (slice 2); `dbartsSpec()` refuses every hazard family, as today |
| `model.frame(o)`, `nobs(o)` | new methods on every data object (dec-B677): from the frame where one is kept, decoded from the block where exact otherwise (dec-B685's rule, for every kind) (slice 1) |
| `o@x`, `o@y`, `o@weights`, `o@offset` | stand-ins on a frame-keeping object only (dec-B675, dec-B676); real slots everywhere else, `s$data` included (slice 1) |
| `makeTestModelMatrix(o, newdata)` | reads the coding record on a frame-keeping object (slice 1) |
| `$setForestBasis`, `predict(newdata)`, `predict(bases = )` | unchanged |

## Slice 1: the data object keeps its frame

Tier: R surface, no draw change ([Process by risk](README.md#process-by-risk)); the review runs the mutants,
since no bitwise gate reaches the resolution (see Gates). Files: `R/A_class.R` (slots, validity), `R/data.R`
([`dbartsData`](../../R/data.R), the resolution, the stand-ins, [`makeTestModelMatrix`](../../R/data.R),
`model.frame` and `nobs` methods), the door entries in `R/dbarts.R`, `R/spec.R`, `R/bart.R`, `R/rbart.R`,
`R/xbart.R`, `R/bartcore.R`, `NAMESPACE`, `man/dbartsData.Rd`, a new page for the two methods,
`man/makeind.Rd`, `man/dbartsSampler-class.Rd` (`setData` takes a frame-keeping object), `_pkgdown.yml`,
`inst/NEWS.Rd`.

### Steps

1. Slots: `frame`, `rows` (integer positions, `NULL` when every row is kept in order), `frame.columns`, the
   coding record (the block with no rows, with a missing-value flag per column), the response's record and
   the frozen response, weights and offset expressions; `x`, `y`, `weights`, `offset` widen to hold a
   stand-in. Validity per kind: a frame-keeping object's positions lie within the frame's rows; a
   block-holding one is checked as today. About 75 lines with the stand-ins.
2. The stand-in: one small classed object per slot, set only where a value exists (`weights` and `offset`
   stay `NULL` where none was given, since `is.null()` cannot be intercepted); `dim`, `dimnames`, `length`,
   both subset operators, `as.matrix`, `as.vector`, `mtfrm`, `Ops`, `Math` and `Summary` stop with "the data
   are read with model.frame(o)"; `print` and `format` show one line, so `show(o)` and `str(o)` work;
   `attr()` gives `NULL` (dec-B676).
3. [`dbartsData`](../../R/data.R) takes `keep.data = TRUE`. It builds and checks everything as today, every
   refusal and warning staying at the build (dec-B675), then on a data frame or list stores the frame, the
   positions (the subset value it already evaluates once, less the na.action's rows), the coding record,
   the response's record and the frozen expressions (a column stays a name, a caller's value its value,
   dec-B185), and drops the block and values. Row names only where the frame cannot give them (call 19). An
   environment keeps the block and `ls(e, all.names = TRUE)`. `dbartsData(o, keep.data = )` on a data object
   takes the argument quietly. About 70 lines.
4. The resolution: [`dbartsData`](../../R/data.R)'s own formula branch run again on the stored frame, with
   the frozen expressions, `subset` set to the positions and `keep.data = FALSE`, then the stored na.action
   record put back; so the block, the sparse lifting, `y`, `weights`, `offset` and the Surv attributes
   (`survivalStatus`, `survivalTime`, `survivalTimeOmitted`) come from the build's own code, in the build's
   order (dec-B675, call 20). A flag the branch reads skips the warnings and checks the build already gave.
   One function, called first at each door in the table; reads of a data object's slots before that point
   (in `bart()`, `dbarts()`, the hurdle's entry) move after it. Internal builders (`dbarts()`, `xbart()`,
   `rbart_vi()`, the hurdle) pass `keep.data = FALSE`. About 100 lines; remaking the block by hand instead
   would be 200 or more.
5. `model.frame()` and `nobs()` methods, registered in `NAMESPACE` (dec-B677). On a frame-keeping object,
   the resolution's variables over the rows in, factors as factors, `(weights)` and the `offset()` and
   `(offset)` columns, as `model.frame()` of an `lm` fit. On a block-holding object (`s$data`,
   `keep.data = FALSE`, the matrix interface, an environment's), decoded where exact, under dec-B685's rule
   and dec-B677's "factors kept as factors": the response to its recorded type and levels, a Surv response
   rebuilt, a factor column as a factor, a matrix term's columns joined back; the rest as held, the help
   listing each case. `nobs()` counts every row in, weight-0 rows included; the help says how this differs
   from `lm`. About 65 lines.
6. [`makeTestModelMatrix`](../../R/data.R) reads the coding record and its flags. About 5 lines.
7. Help: the data object's slots, `keep.data`, the two methods; a frame-keeping object reloaded in another
   session evaluates its formula's terms again at its first fit, so the packages its formula calls must be
   installed; a term that draws (`I(x + rnorm(n))`) is evaluated again there, as `predict.lm` evaluates such a
   term on newdata.
8. Tests: the reads of the old slots move to `model.frame()`, to the made object or to `keep.data = FALSE`. A
   static scan finds 83 to 87 lines in 14 to 17 files reading `@x`, `@y`, `@weights` or `@offset` of a named
   formula-built object; the note's bound by name is about 360.

### Pins (fail before, pass after)

- The new slots, `keep.data = FALSE` and `dbartsData(o, keep.data = FALSE)` taken.
- The positions under each `subset` and na.action of the note's section 2.1, a repeated row, a character
  `subset` and a `subset` holding `NA` (lm keeps the padded row for the na.action) among them.
- The stand-ins: each method above stops naming `model.frame(o)`; `attr(o@x, "factor.levels")` is `NULL`;
  `is.null(o@weights)` where no weights were given; `print(o)` and `str(o)` run.
- `model.frame(o)`: `model.response()`, `model.weights()` and `model.offset()` equal those of the `lm` fit on
  the same call; `nobs(o)` keeps rows of weight 0. On each block-holding kind, a factor response and a
  factor predictor come back as factors.
- `d$z[1] <- 9` after the build leaves `o@frame$z` and the fit unchanged; a later change of `e$x1` leaves an
  environment's object unchanged; a dotted name of `e` is in `frame.columns`.
- `makeTestModelMatrix(o, newdata)` on a default object equals it on `keep.data = FALSE`.
- Saved size: a default object against `keep.data = FALSE`, and a sampler within 0.001 MB of main's.

### Guards (pass before and after)

- The sibling identity, for each case of the note's section 2.3 list and each response kind (gaussian; a
  factor and a logical response under probit; an ordinal factor; a multinomial factor through `dbarts(o)`;
  `aft` and `hazard` from a Surv response): `dbarts(dbartsData(f, d, ...))` and `dbarts(f, d, ...)` draw
  identically, and the resolution's block, `y`, `weights`, `offset` and Surv attributes are `identical()` to
  those of `dbartsData(f, d, ..., keep.data = FALSE)`. The same for `bart(o)`, `xbart(o)`, `rbart_vi(o)` and
  the hurdle beside their formula siblings, and for `new()` and `$setData(o)`.
- No repeated warning: under `withCallingHandlers`, a zero-weight build warns once at `dbartsData()` and
  each door in the table then gives none.
- Every door in the table takes a frame-keeping object; `s$data` and `dbartsSpec(o)$data` hold no frame.

Mutants, which the reviewer runs: the resolution evaluating terms on the kept rows (predict's path; the
caller's vector under `subset` fails); the na.action record not put back (`na.omit` fails); one door left
without the resolution (`$setData(o)` fails); the Surv attributes not remade (`aft` fails); the build's
warnings given again (the warning guard fails).

### Gates

Lint (`lintr::lint_package()` against the slice's library), air, codoc, win-drift, doc-freshness; the full
suite once, locally, in one process (138 test files build a data object, and no bitwise gate reaches the
resolution: none of the four snapshot files calls `dbartsData(`, the equivalence trio's only data object is
a matrix-interface one at `setData`, and no exact gate builds `dbartsData(f, d)`); the four snapshot files on
the reference build; the equivalence trio bitwise on the reference build; the exact gates with `quick` (a
fit now carries new data slots); `pkgdown::check_pkgdown(".")`; NEWS parse; `R CMD check --as-cran` from a
clean tarball. CI runs the full suite, check and sanitizers.

## Slice 2: the forest doors

Tier: changes draws, for the named calls only; no engine change. Files: `R/data.R` (`bases` out of
[`dbartsData`](../../R/data.R), [`validateForestBases`](../../R/data.R)'s cut moved,
[`alignSubsetRowsToFrame`](../../R/data.R) and [`restrictBasesToRows`](../../R/data.R) retired),
`R/forestBasis.R` (the shared reader beside [`readForestBasis`](../../R/forestBasis.R);
[`basisRowNumbers`](../../R/forestBasis.R) and [`fitBasisRows`](../../R/forestBasis.R) retired),
`R/formulaTerms.R` ([`ingestFormulaTerms`](../../R/formulaTerms.R)), `R/model.R`
([`resolveFormulaBasisSubset`](../../R/model.R) retired), `R/dbarts.R`, `R/spec.R`, `R/bart.R`,
`man/dbartsData.Rd`, `man/dbartsSpec.Rd`, `man/dbarts.Rd`, `man/forest.Rd` where it names the data door,
`man/dbartsSampler-class.Rd`, `inst/NEWS.Rd`.

### Steps

1. One reader, `resolveForestBases(forests, data, rows, held)`: reads each declared basis against the data's
   frame (a frame-keeping object's, or the `d` `dbarts(f, d)` and a formula term hold), then the columns the
   object holds (its block's, and its held bases under open call 3), then where `forest()` was called (open
   call 2 for objects with no frame that dec-B682 does not name); a value has one entry per row handed over
   and is cut by `rows`, code is built on the kept rows (dec-B680); a hazard fit's basis per open call 4.
   `dbarts(f, d)`, a formula's `forest()` term, `dbarts(x, y)`, `dbarts(o)` and `dbartsSpec(o)` call it; the
   row-number plumbing through `bases` and `basisRows` and the restating of "'bases'" as "'basis'" go. About
   60 lines added; about 270 removed (the `dbarts()` plumbing about 95, the `dbartsSpec()` block about 55,
   the five helpers about 90, the `bases` sites in `dbartsData` about 30).
2. `dbarts(o, forests = )` takes a declaring list, read exactly as `dbartsSpec(o, forests = )` reads it
   (dec-B681). About 25 lines.
3. Beside an object with no caller frame (`s$data`, `keep.data = FALSE`), a basis naming a column the object
   holds (a block column, or a held basis, open call 3) is read from it; a name in `frame.columns` that is no
   such column is refused, naming the data object built with its frame; any other name is the caller's value
   (dec-B682). Beside an environment's object, the same, the refusal naming a data frame or list (dec-B683).
   About 15 lines.
4. Beside any data, `s$data` included, a forests list states the forests in full: a forest declaring no
   basis has none, a list of tree counts alone gives plain forests, no list gives one forest (dec-B684). Held
   bases no forest names are dropped from the data handed on, with their `basis.levels` and `basis.scale`
   entries, since a creation of fewer basis forests than the data carry is refused today. A forest naming a
   held basis takes its scale per open call 1; a basis given as a new value takes its scale and levels from
   that value, as a declaration does today. dec-B663's interim refusal is not built. About 15 lines.
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
- D4. Beside a sampler's data whose basis `$setForestBasis` swapped, a forest naming the held basis: under
  open call 1 (a) the scale recorded at the first creation, as today's no-list route gives, no prior moving;
  under (b) the column's spread now, the coefficient's prior moving by the ratio of the two spreads.

Calls that run today and stop: beside a frame-keeping, `keep.data = FALSE` or matrix-interface object under
`subset`, a basis of the kept length (dec-B680); under open call 4 (a), a hazard basis of one value per period
row on the matrix route. A formula hazard fit with a basis, refused today, fits under (a) or (b).

No test file and no snapshot moves. A parse of every test file finds 62 door calls with a forests list beside
a data object: the frame-keeping ones name no column of their frame, the one sampler's data of several
forests (`original$data` in test-forest-scale-record.R, swapped) takes a list stating both forests, its
second a new value (`bigDose`), which takes its scale from that value under either answer to open call 1,
and the one name collision, ["plainData <- dbartsData(y ~ x1 + x2 + x3,
frame)"](../../inst/tinytest/test-control-doors.R) beside `factor(z)`, has the workspace `z` equal to
`frame$z` and is refused for its tree count. The doors handed a sampler's data with no list (16, in
test-control-doors.R, test-standardization-record.R and test-leaf-standardization-rows.R) are all of one
forest. The four snapshot files build no data object with a forests list. The implementer confirms each by
running the touched files against main's library and the slice's.

### Pins (fail before, pass after)

- D1, D2, D3 and D4 each, with the caller's variable unlike the column, the new draws `identical()` to the
  sibling door's: `dbarts(f, d, forests = L)` for D1, a fit naming the held values for D2,
  `dbarts(f, d, forests = <the list as stated>)` for D3 (each building a sampler, so the dropped held bases
  are exercised), and for D4 the forest's recorded scale after a swap as open call 1 rules.
- The held basis kept by name: `dbartsSpec(s2$data, forests = <a list naming the held basis>)` on an
  unswapped sampler identical to today's `dbartsSpec(s2$data)`, the name as open call 3 rules.
- dec-B680 at `dbarts(o)`, `dbartsSpec(o)`, a `keep.data = FALSE` object and a matrix-interface object
  under `subset = 1:40`: a 60-row vector taken, equal to `dbarts(f, d, subset = 1:40, forests = )` and
  `dbarts(X, y, subset = 1:40, forests = )`; a 40-row one refused with the length expected.
- A hazard basis, formula and matrix routes, as open call 4 rules, each refusal naming the length.
- dec-B681: `dbarts(o, forests = <declaring>)` taken, identical to the `dbartsSpec()` route and to
  `dbarts(f, d, forests = )`.
- dec-B682 and dec-B683: each refusal's text, and the caller's value for a name of neither kind.
- `dbartsData(bases = )` gives "unused argument".

Mutants for the review: the reader skipping the frame beside a frame-keeping object (D1 fails); a held basis
kept beside a list of tree counts (D3 fails); a value cut by the kept length (dec-B680 fails); the term door
left on `fitBasisRows` (the equivalence trio's two-forest formula scenarios fail).

### Gates

Slice 1's, with the equivalence trio bitwise the gate for D's absence elsewhere (bart2twoforest, the
two-forest formula scenarios and the 19 bcf scenarios read their bases through the new reader), and the full
suite once, locally, in one process.

## Slice 3: the consumers

Tier: R surface, in each consumer's repository, prepared against slice 2's library and pushed to the
consumer's main right after dbarts' merge (the patches name an argument main does not have before it).
Under dec-B674's examination clause, from the examination recorded in the note's section 2.5:

- bartCause passes `keep.data = FALSE` at its two data call sites (three objects), on its merits: its saved
  fits keep `data.rsp`, which would otherwise carry the user's whole frame, and it reads the block
  (`nrow()`, `colnames()`, `x[, z]`, `as.matrix(x)`), `@y` and `@weights`, and writes `@x.test`. `bcf()`
  drops `bcf.bases` and gives the treatment forest `basis =` the full-length `factor(z)` it builds, cut by
  the object's `rows` under `subset` (dec-B680). The other route, naming the treatment column
  (`forest(basis = factor(<treatment>, levels = c(0, 1)))`, read from the block under dec-B682), would drop
  that vector, but on bcf's matrix route it works only if open call 2 takes (b); the value route works on
  both routes under either answer, so it is taken, and is examined again if (b) is ruled. Its reads of
  `attr(fit$fit$data@x, "term.labels")` stay (the block here, dec-B686 in part 2). About 8 lines. Pins:
  `bartc()` and `bcf()` fits identical before and after, with and without `subset`; its suite.
- stan4bart passes `keep.data = FALSE` in `glFormula`, because it reads the coded block (`as.matrix(@x)`,
  `@varTypes`, `factor.levels`, `makeTestModelMatrix()`), which no accessor gives from a frame-keeping object
  (keeping lme4's frame would cost it almost no size), and because its hand-written slots (`@weights`,
  `@n.cuts`, `@useQuantiles`, `@x.test`, `@testUsesRegularOffset`) would meet stand-ins; the second couples
  this line to root TODO consumer-slot-writes. Its bundle is dec-B687's, in part 2. 1 line. Pin: its suite, a
  fit identical.
- treatSens uses the matrix interface and `dbartsSpec(o, forests = list(forest(n.trees = )))`, reading no
  slot; under dec-B684 one forest is still one forest: nothing.

## Records

At the landing: ledger entries for the calls made here, this plan's Landing note and Status, the design's
Status line, TODO, the plans and design indexes; the four design notes that name the data door
(forest-kind-by-class, public-surface, written-surface, multiplier-combiner). NEWS, one entry against
0.9-34: "A data object made by dbartsData() from a formula and a data frame keeps the frame and the
positions of the rows it keeps; its x, y, weights and offset slots, which held the model matrix and the
response in 0.9-34, no longer hold them: model.frame(object) returns the model's variables, factors as
factors, and nobs(object) counts the rows. A forest given to dbarts() or dbartsSpec() beside a data object
finds its basis in the frame. keep.data = FALSE keeps the matrix and the response in the object instead, as
before." Part 2 extends it.

## Verification

- `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts", lib.loc = <lib>)'`, once, in one process: 0
  failed; each touched file by `tinytest::run_test_file()` with `library(dbarts)` attached; the new pins fail
  against main's library.
- `Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-38da480c.rds --bitwise
  --strict-coverage`, and the bcf and multinomial harnesses against their current files, on a `--preclean`
  reference install: 55, 19 and 11 "identical draws" lines, no "max |z|" line.
- The four snapshot files on the reference build: 27 results, 0 failed.
- The exact gates of `.github/workflows/exact-gates.yaml` with `quick`: all pass.
- Each gate of [Process by risk](README.md#process-by-risk) on its own exit status.

## What part 2 reshapes

The resolution's product (today's block-holding object becomes the model frame, about 30 lines of slice 1's
100); `model.frame()` on a sampler's data, kept to dec-B685's cases; the held basis's name (open call 3,
beside call 28's `(basis.<label>)` columns, about 10); an environment's object (dec-B683's model frame,
about 5); the stand-ins on `s$data` (dec-B686); the mask (dec-B688); stan4bart's three reads and its bundle
(dec-B687). Between the two landings main holds both shapes (dec-B679).

## Open calls

Flagged for the maintainer; none is decided here. Ordered by which bites first.

1. **A held basis's scale after `$setForestBasis`** (D4). Under dec-B684, naming the held basis is the only
   way to keep it. Shown: (a) the forest naming it keeps the scale recorded at the first creation, matched by
   name, as today's no-list route gives, about 5 lines, recommended: the prior rule is "taken once, moved by
   no later call"; (b) every declared basis takes its scale from the column it reads, the prior moving by the
   ratio of the spreads (2.52 to 25.2 in the premise).
2. **Objects with no frame that dec-B682 does not name** (the matrix interface's, and a formula's with no
   `data`). Shown: (a) a basis name read where `forest()` was called, as their building door reads it, so
   `dbartsSpec(dbartsData(x, y), forests = L)` draws as `dbarts(x, y, forests = L)`, recommended; such an
   object and its own sampler's data then read a column name differently (the workspace, the sampler's
   column); (b) dec-B682's rule for every frameless object, a column of `x` read from it, `dbarts(x, y)`
   differing unless its door changes too; bcf's matrix route could then name its treatment column. `lm`
   refuses a matrix as `data`.
3. **How a held basis is named beside a sampler's data.** dec-B684's example `forest(basis = z)` rules out
   naming by the forest's label. Shown: (a) by the expression it was written as, recorded with the held basis
   (`forest(basis = factor(z))` keeps one written so; a source column of the caller's frame is refused under
   dec-B682, the message naming the written form), about 10 lines, recommended; (b) by its source columns,
   `s$data` keeping each (`z`), so any expression of them is read again, about 15 lines and one column a
   source. Goes with part 2's call 28, `(basis.<label>)`, so both parts answer one way.
4. **A hazard basis's length.** dec-B680 says one entry per row handed over; a hazard caller hands over
   subjects, which dbarts expands to period rows. Shown: (a) one value per subject, repeated over its periods
   as weights and offsets are, the formula route starting to work and the matrix route's period-length
   vector refused, about 10 lines, recommended; (b) both lengths, a period-length vector a basis that varies
   over time; (c) as built, one per period row, the formula route unable to take a basis.

Brought to the maintainer, not open calls: the size, about 450 lines of R against dec-B679's 305 (the
formula term door, the response kinds, hazard and `model.frame()` of every kind were not in that count);
and the cost of a formula term that draws, evaluated again at the resolution, so `dbarts(dbartsData(f, d))`
fits other values than `dbarts(f, d)` and consumes the generator, a cost the rulings did not weigh (the help
sentence is an agent call; keeping each term's columns instead would revise dec-B670 and dec-B675).

## Calls made in planning

Agent calls, each with its alternative, for the session's close-out (dec-B629); the landing does not wait
for them:

- The note's agent calls this part builds on: 12 (positions, nothing stored when every row is kept in
  order), 13, 14 (a forest's record names columns), 15, 16 (`bases =` gets "unused argument"), 17, 18
  (`x.test` kept as built), 19 (row names only where the frame cannot give them), 20 (the resolution in the
  build's order) and 22, each with the alternative the note gives.
- `keep.data = FALSE`, an environment's and the matrix interface's objects keep `rows`, since dec-B680's
  length rule at every building door needs it (bartCause's `bcf()` under `subset`). Alternative: none that
  keeps dec-B680.
- `dbarts(f, d)` builds with `keep.data = FALSE` internally and never runs the resolution, so no formula
  fit depends on it. Alternative: `dbarts(f, d)` as `dbarts(dbartsData(f, d))`, one path, every formula fit
  then remaking its block (0.004 s at 20,000 rows) and depending on the resolution.
- The resolution re-runs the build's formula branch on the stored frame rather than remaking the block by
  hand. Alternative: a separate maker, 200 lines or more, which can drift from the build.
- An absent weights or offset slot holds `NULL`, the stand-in only where a value was given, since NULL means
  absent and `is.null()` cannot be intercepted. Alternative: a stand-in in every slot, `is.null(o@weights)`
  then FALSE with no word.
- The build's warnings are not given again at a door; the resolution skips them. Alternative: muffling every
  warning at the resolution, which would hide a new one.
- An environment's names are `ls(e, all.names = TRUE)`, so a dotted name is refused under dec-B683 rather
  than read as the caller's. Alternative: `ls(e)`.
- Slice 1 keeps the R-surface tier with a sonnet review that runs its mutants, since no draw is meant to
  move and its guards are sibling-door identities. Alternative: the opus review, as the critique offered,
  since no bitwise gate reaches the resolution.
- Slice 2 takes the changes-draws tier for its named moves (a slice spanning two tiers takes the higher).
  Alternative: the R-surface tier, the sibling-door oracle standing in for the opus review.
- bartCause's `bcf()` takes the value route for its treatment basis (Slice 3). Alternative: naming the
  treatment column, open to its matrix route only under open call 2 (b).
- The consumer patches land after dbarts' merge, pushed at once. Alternative: patches that test
  `"keep.data" %in% names(formals(dbarts::dbartsData))`, mergeable first, a branch on a development build.
