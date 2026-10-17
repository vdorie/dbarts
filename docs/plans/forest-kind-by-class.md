# forest-kind-by-class: what a forest multiplies is decided by the class of its basis, at every door, fixed at creation, and is one numeric column or two levels

Status: PLANNED, revised from its blind critique on 2026-10-11 with no call open
([Revised from the blind critique](#revised-from-the-blind-critique)); written to the tree of 37b59ea0,
the final tip of the sd unit slice. Rulings: dec-B260; dec-B254 for the kind as structure; dec-B289 and
dec-B308 for a factor's levels; dec-B346 for a column of zeros; dec-B281 and dec-B282 for the one size;
dec-B452 for a swap that never changes the kind; dec-B453 for numbers refused for a factor forest at
predict; dec-B231 and dec-B418 for a state that carries no data. No call is open
([Open calls](#open-calls)). Follows [forest-sd-unit.md](forest-sd-unit.md), whose tip that is and
which lands first, and [forest-defaults-by-kind.md](forest-defaults-by-kind.md) and
[written-surface.md](written-surface.md), landed.

agent: two pushes, both of the tier that changes draws ([Process by risk](README.md#process-by-risk)):
push 1 because four call sequences fit another model after it, push 2 because it rewrites and re-records
a baseline scenario. The plan has had its blind critique; each push takes an implementer and one opus
review that runs the mutants and checks the fix round. Push 1 (the rule): opus implementer for the R
code, sonnet for the help once the code is fixed. Push 2 (one size a forest), in two commits and a
baseline commit: sonnet implementer for the first, which respells; opus for the R code of the second,
which refuses, and sonnet for its tests and help. The reason for opus on push 1: it moves no arithmetic,
and every slip in it is a forest read as one kind at one door and another at the next, or a coefficient
multiplying another level's rows after a swap, with no message. The tip has both faults today (Context),
and the three reviews of written-surface's second push each found a fault of that shape one step to the
side of the last. On the refusals: a door left out is a model 1.0-0 was ruled not to fit, accepted in
silence. Nothing under src/ changes, so no new numerics are reachable and no local sanitizer run is owed.
rng: stated per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- Push 1. NEUTRAL for every fit that both builds accept, the four sequences below apart. Read off the
  tip and the rule, no build of the rule existing: of the 259 cells of Context's first table (37 kinds of
  object at 7 creation doors) 156 are created today and keep their kind, their block and so their draws,
  78 are refused and stay refused, 10 are created where they are refused, 13 are refused where they are
  created (6 of them the character matrix of one column, 3 a factor with NA as a level handed over as a
  value), and 2 are the first sequence below. Of the 28 cells of a number and a factor with a missing
  value on a row the fit drops (the 7 doors such a row reaches, under `subset` and under a missing
  response) 15 are created today and keep their draws, and 13 are created where they are refused: 9 for
  the rows a missing value is judged on, 4 the factor at the data door. Three refusals have no cell,
  their objects being none of the 37: `$setData` on a sampler of one forest handed a data object that
  carries bases; at the data door a basis with a row whose norm is not representable and a matrix of no
  columns, which the forest doors refuse already; and a factor matrix, at every door.
  POSTERIOR-CHANGING for four sequences, none in a test, a benchmark, a help example or a consumer
  (traced over the test files that name a basis, and searched elsewhere):
  1. a logical vector handed to `dbartsData(bases = )`: one 0/1 numeric column before, two levels after,
     which is what the same vector is at every other door. The oracle is the base build's fit of it
     through a `forests` list: identical draws.
  2. `$setForestBasis` on a forest created with a factor or a character vector, given one whose levels
     stand in another order or lack one of them: matched by position before, by name after (dec-B289).
     The oracle is the base build given the same vector as `factor(v, levels = <the levels at creation>)`.
  3. `predict(bases = )` with such a vector: the same, with the same oracle.
  4. a factor or a character vector of three or more levels handed over as a value, where `subset` or the
     na.action leaves one of the levels with no row and two or more with one: a column of zeros with a
     coefficient of its own before, the levels that have a row after (dec-B308). The oracle is the base
     build's fit of the same factor written as code, which drops the level already: identical draws.
- Push 2, first commit. No draw moves: each respelled call is the same block of 0s and 1s (at the tip
  `cbind(1 - z, z)` at the seven doors and `factor(z)` at the six that take it are one sum of draws).
  One recorded scenario is rewritten, because the second commit refuses its model:
  ["bart2twoforest"](../../benchmarks/R/equivalence.R) of the gaussian equivalence baseline multiplies
  a forest by a factor of three levels. It is given two, recorded afresh and merged into a copy of the
  baseline; the other 54 scenarios are bitwise.
- Push 2, second commit. No draw of a model both builds accept moves. Refused where accepted: a factor,
  character or logical basis of three or more levels; a basis of several numeric columns, however
  written; a swap or a prediction that gives a forest either. At the traced commit the tests create 222
  such forests in 11 files (Context); each is rewritten onto a shape that stays or goes with the test of
  a shape that does not.
window: pre-release, after the sd unit and before the multiplier slice, which gives a factor and a number
different priors: until it lands a two-level factor and its indicator columns handed over as numbers are
one fit, so this is the only tip on which every such block in the tests, the benchmarks and bartCause
can be respelled and the respelling proved bit for bit. Serial with any other work in
[`expandForestBasis`](../../R/model.R), [`basisColumns`](../../R/forestBasis.R),
[`buildFitBases`](../../R/forestBasis.R), [`replacementForestBasis`](../../R/forestBasis.R),
[`validateForestBases`](../../R/data.R), [`setForestBasis`](../../R/dbarts.R),
[`resolveForestBases`](../../R/generics.R), [`extractForest`](../../R/generics.R),
[`packageBartResults`](../../R/bart.R), [`bartcoreSamplerSetData`](../../R/bartcore.R) or the
multi-forest block of [`resolveSamplerSpec`](../../R/spec.R). No engine, bridge or header file
changes. bartCause's edit lands between the two pushes, after the sd unit's own edit to the same files
there.
budget: ~4350 lines planned, forecast 6600 to 8700; stop at ~8700, or at ~885 of R code
([Budget and stops](#budget-and-stops)). Push 1 ~1965: R ~410, tinytest ~1135, of which ~840 in one new
file, benchmarks ~5, help ~185, records ~230. Push 2 ~2385: its first commit ~445 (tinytest ~380,
benchmarks ~40, help ~25); its second ~1925 (R ~180; tinytest ~1375, of which ~425 in one new file and
~950 rewritten or removed in 11, test-forest-basis-terms.R most of it; help ~210; design notes, TODO and
the two indexes ~160); its baseline commit ~15 beside the recorded file.

## Goal

A forest multiplies nothing, the levels of a factor, or a numeric column, and which of the three follows
from the class of what the caller handed over and from nothing else: a factor, a character vector or a
logical vector is levels; a number is numbers, whatever its values. One object is one kind at every door,
`dbartsData(bases = )` included. The kind and a factor's levels are recorded when the sampler is created
and hold for its life: the levels are those with a row among the rows the fit keeps, in every spelling;
`$setForestBasis` and `predict(bases = )` take the class the forest was created with, a factor is
matched to the recorded levels by name, and a forest created without a basis takes none. What a row the
fit drops holds, a missing value included, decides nothing, at any door. A fit stores for a forest of
levels the factor itself, so that what it stores is what `predict(bases = )` takes. In 1.0-0 what a
forest multiplies has one size: a basis is one numeric column, or a factor, a character vector or a
logical vector of two levels. Three or more levels, and several numeric columns however they are
written, are refused by name at every door, with the form that works: a forest for each level but the
first, or for each column.

## Context

Measured on the shipped build of 37b59ea0, R 4.6.1, the trace of the tests apart, which ran on the same
slice before its rebase and is rerun as step 1.0: 150 rows; x1, x2; z a 0/1 double with 66 zeros and
84 ones, and the same column as an integer, a logical, a factor, an ordered factor and a character
vector; g a character vector of three levels; dose, age. "Draws" is the sum of 10 seeded train draws of
a two-forest sampler, 15 trees; the same number is the same fit. A basis is written as code where a call
of `forest()` names it, `forest(basis = B)` with B a variable or a column, and is a value where it is
handed over as an object, `do.call(forest, list(basis = B))`.

- The kind at creation, today. Seven doors: a `forest()` term of a formula at `dbarts()` and at
  `bart()`; a `forests` list on a formula with its basis written as code, and with it handed over as a
  value; the list on the matrix interface; the list in `dbartsSpec()`; and `dbartsData(bases = )`, which
  `dbarts()`, `bart()` and `dbartsSpec()` then take, and `bartBT()` where a response is handed over
  beside it. The first three were run with the basis as code, the next three with it as a value. 37
  kinds of object, 259 cells, 88 of them refused.

  | handed over | as code | as a value | the data door |
  |---|---|---|---|
  | a factor; an ordered factor; `I()` of one; a character vector | levels: 2 columns, draws 1844.053557 | the same fit | refused: "'bases' must be numeric or logical" |
  | a logical vector; `I()` of one | levels, the same fit | the same fit | one numeric 0/1 column: draws 1929.997248, another model |
  | a 0/1 integer; a 0/1 double; `I(z)`; `matrix(z)`; `as.numeric()` of a logical; a time series | numbers, one column, draws 1929.997248 | the same | the same |
  | two values that are not 0 and 1 (-0.5, 0.5; 1, 2); three values; dose | numbers, one column | the same | the same |
  | `cbind(1 - z, z)`, held in a variable | numbers, 2 columns; draws 1844.053557, the factor's own | the same | the same |
  | `model.matrix(~ g - 1)` | numbers, 3 columns; draws 1823.835081, those of g as a factor | the same | the same |
  | `model.matrix(~ g)[, -1]`, a level omitted; `cbind(dose, age)` | numbers, 2 columns | the same | the same |
  | g as a factor; as a character vector | levels, 3 columns, draws 1823.835081 | the same | refused, as not numeric |
  | a logical matrix, of one column or of two | refused: "'basis' (B) has a term that is a logical matrix: a basis is a factor, a character or logical vector, or numeric" | refused: "a 'basis' must be numeric, a factor, or a character vector" | numbers, 1 and 2 columns |
  | a character matrix of one column | levels, the factor's fit | the same | refused, as not numeric |
  | a data frame, of one column and of two; a list | refused: "invalid type (list) for variable 'B'" | refused, as the logical matrix is | refused, as not numeric |
  | a Date; a complex | refused: "has a term that is a Date" | refused, as the logical matrix is | refused, as not numeric |
  | a factor with a level no row has | created: 2 columns, the level dropped | refused: "a 'basis' factor level with no observations contributes nothing to a forest: '2'; drop it with droplevels()" | refused, as not numeric |
  | a factor of one level | refused: "a 'basis' factor must have at least two levels" | the same text | refused, as not numeric |
  | a logical that is `TRUE` on every row | refused: "'basis' (B) is TRUE on every row the fit keeps, so its other level has no observations and the forest would be multiplied by a constant" | refused with the droplevels text, for 'FALSE' | one constant numeric column |
  | a number with a missing value | refused: "a 'basis' cannot be NA" | the same text | refused: "'bases' values must all be finite" |
  | a factor with NA as a level and a row at it, `factor(x, exclude = NULL)` | refused: "'basis' (B) has the level 'NA' at a new row, which no row of the fit had" | created: 3 columns, one of them NA's | refused, as not numeric |
  | `rep(1, n)` | numbers, one constant column | the same | the same |
  | `rep(0, n)` | refused: "a 'basis' column of all zeros contributes nothing to a forest" | the same text | accepted: a forest multiplied by zero |

  So a 0/1 number is one numeric column at every door today, and is not the factor's fit; dec-B260
  changes nothing for it. What the tip decides by something other than class is the data door, and
  between code and a value what an unused level does. Code names its columns as `model.matrix()` names
  them (`B0`, `B1`; `Bdose`, `Bage`), a value keeps the names every column has or none, and the data
  door keeps what it is given, `cbind(1 - z, z)` arriving with one column named.
- What code gives, by shape. `gf`, a factor of three levels, 3 columns (`gfa`, `gfb`, `gfc`);
  `dose + age` 2; `splines::ns(dose, 3)` 3; `poly(dose, 2)` 2; `1 + dose` 2 (`(Intercept)`, `dose`); and
  one column each for `dose`, `0 + dose`, `dose - 1`, `scale(age)`, `poly(dose, 1)`, `dose:age` and
  `I(dose * age)`. A matrix written in place, `forest(basis = cbind(dose, age))`, is refused by the
  grammar ("the columns of a basis are separated by '+'; write dose + age"); a matrix held in a variable
  and named in code is one term of several columns.
  `forest(x1, basis = gf == "b") + forest(x1, basis = gf == "c")` beside a forest with no basis is
  created: three forests, each of the two a logical's two levels. So is
  `forest(x1, basis = dose) + forest(x1, basis = x2)`. A code basis travels through the data object as
  one column of row numbers ([`basisRowNumbers`](../../R/forestBasis.R)) and is built on the rows kept
  afterwards ([`buildCodeBasis`](../../R/forestBasis.R)); its record for new rows holds the levels kept.
- Rows. Under a `subset` that leaves one level of g with no row, and under a missing response on those
  rows:

  | basis | door | result |
  |---|---|---|
  | the factor column, `factor(g)` or the character column, written as code | a formula term; a list on a formula; a list on the matrix interface | 2 columns, the level dropped; draws 1300.647320 |
  | the factor, or the character vector, handed over as a value | a list on a formula; the matrix interface | 3 columns, one of them all zeros, accepted; draws 1308.622837 |
  | the factor already cut to the kept rows, handed over as a value | `dbartsSpec()` over a data object made under the `subset` | refused, the droplevels text |
  | its indicator columns | the data door | 3 columns, one all zeros, accepted; draws 1308.622837 |
  | a level emptied by a missing response | code; a value; the data door | as under `subset`: 2 columns; 3 with one all zeros; the same |
  | a factor of two levels that the rows leave with one | code; a value | refused: "a 'basis' factor must have at least two levels" |
  | a number left all zero by `subset` | code; a value in a list or on the matrix interface | refused, the zero-column text |
  | the same | the data door, under `subset` or a missing response | accepted: one column of zeros |
  | a block of several numeric columns of which `subset` empties one | written as code | refused, the zero-column text |
  | the same | a value in a list | accepted: 3 columns, one all zeros |
  | a logical left all `TRUE` | code; a value | refused, the text for a logical, with and without the basis's label |
  | the same | the data door | one constant numeric column |

  A factor of two levels whose rows at one level all have weight zero is created with both levels, as
  code and as its indicator columns through the data door: a row of weight zero is a row the fit keeps.
- A missing value on a row the fit drops, today. A number and a factor, each with one missing value on a
  row that `subset` drops, and on a row a missing response drops, at the seven doors such a row
  reaches. `dbartsSpec()` has no such row: a value handed to it has the kept rows' length, and one of
  the data's full length is refused for its length. 28 cells, 13 of them refused.

  | door | the number, either way | the factor, either way |
  |---|---|---|
  | written as code: a term at `dbarts()` and at `bart()`, a list on a formula | created, the fit of the complete column | created |
  | a value in a list on a formula; on the matrix interface | refused: "a 'basis' cannot be NA" | refused, the same text |
  | the data door on a formula | created | refused, as not numeric |
  | the data door on the matrix interface | created under `subset`; refused under a missing response: "'bases' values must all be finite" | refused, as not numeric |

  Code reads its model frame on the rows kept; a value at a forest door is checked whole before any row
  is cut; the data door checks after `subset` and the na.action on a formula, and on the matrix
  interface after `subset` and before the rows of a missing response go. With the missing value on a
  row the fit keeps every door refuses.
  `lm` judges a variable's missing values after `subset`. At `predict` a missing value in `bases` on a
  row of `newdata` that the na.action drops is passed by and one on a row it keeps is refused, as code
  is; a swap has no dropped row.
- NA as a level. Besides the row of the first table: `addNA()` of a complete factor, the level NA with
  no row on it, is created as code with the two levels that have rows, and as a value is refused with
  the droplevels text, as any level no row has. At a swap a factor with a row at its NA level is taken
  and the forest has a third column. `lm` fits a coefficient for the level and predicts at it.
- After creation. `$setForestBasis` takes any class on any forest: 8 kinds of forest by 17 kinds of value
  were tried, and every pair is accepted except a factor of one level, a logical matrix, and a
  replacement of another width on a held forest
  ([`refuseHeldWidthChange`](../../R/model.R)). A forest created with no basis takes one, of any kind,
  wherever it stands; a factor forest takes dose, or `cbind(dose, age)`, or a column of zeros
  (dec-B346, built); a numeric forest takes a factor and has two columns afterwards; the width follows
  the value. A one-sided formula is read as code, and `~ dose` on a factor forest is taken. A factor
  forest takes its own stored block, `data@bases[[f]]`, handed back. The sd unit slice pins that a
  forest created with no basis and given a column keeps the prior it was created with, and records
  `hasBasis` for it.
- A factor's levels are matched by position. A held character forest (control, treated), its
  coefficients 0 and 1: after a swap to the same vector as a factor with levels (treated, control) the
  block's second column, the one the coefficient 1 multiplies, is 1 on the control rows.
- `predict(bases = )` expands what it is given by its own levels and checks the width ("'bases' gives
  forest 2 1 column; its amplitudes take 2"). On a fit created with `factor(z)`: the same factor with
  its levels reversed is accepted and predicts 1.0739 where the aligned one predicts 1.2976; a two-level
  character vector unrelated to z is accepted; `cbind(1 - z, z)` predicts what the factor does;
  `cbind(dose, age)` is accepted (mean prediction 23.1). On a fit created with a 0/1 number a factor of
  one level is accepted as a column of ones, and a column of zeros is accepted, the prediction at
  z = 0. With no new rows `bases` is refused by name. `bases` is a list, or the bare value where one
  forest has a basis. `predict(fit, x, bases = fit$bases)` at the training rows is taken and is the
  stored training fit to 3e-15, the two summing in another order.
- What is recorded. A data object has `bases`, bare numeric matrices, and no record of levels. The
  forests' record on the control has `hasBasis`, `labels` and, for a basis written as code,
  `basisTerms`, whose entry holds the levels kept. A fit carries `bases` and `basis.terms`
  ([`packageBartResults`](../../R/bart.R), [`replayForestBasis`](../../R/model.R)); its `bases` is the
  data object's list of blocks, read in R by [`predict.bart`](../../R/generics.R),
  [`extractForest`](../../R/generics.R) and [`resolveForestBases`](../../R/generics.R). A sampler has no
  method that returns a forest's basis or its levels. A sampler made again from its own control, model
  and data, a copy and a reload carry the same blocks, also after a swap that left a level or a numeric
  column with no row.
- What is checked when a sampler is made. The class's validity runs where
  [`dbartsData`](../../R/data.R) makes the object and nowhere else: the sampler's `initialize` stores the
  data object it is handed and creates, and `copy()` goes through `initialize`; a reload's re-creation
  ([`recreatePointer`](../../R/dbarts.R)) does not. So a data object edited by hand to three numeric
  columns is created by `new("dbartsSampler", control, model, data)` and by `dbarts(data)`, and one with
  a slot removed, which the validity refuses, is created by both and by `dbartsSpec(data)`. `$setData`
  and `$setModel` refuse a sampler of several forests whole. On a sampler of one forest `$setData`
  takes a data object that carries bases with no message: the run is the plain model's, and `copy()`
  then fails with "the data carry forest bases but no basis forest was configured".
- What the tests create, traced on the sd unit slice as it stood at e88b78eb, before response-scale-rows'
  third slice and the sd unit's review round were under it; step 1.0 runs the trace again, and every count
  from here to "Elsewhere" is of that earlier tree. On 37b59ea0 one more file names a basis, 73 with
  test-update-scale-rows.R, whose forests multiply one numeric column, as code and through the data
  door; `$setForestBasis` is named on 86 lines of the same 15 files, the two new ones in
  test-forest-sd-unit.R, one of them a call, a number given to a numeric forest; test-bcf-creation.R
  has 15 more lines and test-forest-sd-unit.R 120. The trace: the 72 test files that name a basis, each
  run in its own process with the functions that read a basis wrapped; 10788 results. 1412 models of
  several forests in 67 files (1177 of two forests, 221 of three, 14 of four). Per forest with a basis,
  1747 in all: two levels 447,
  of which 5 handed over as a value (3 factors, 2 logical vectors); one numeric column 987 (a 0/1
  column 360, a constant column 31, another 596); two indicator columns handed over as numbers 91, of
  which 78 through the data door; three levels 89; three indicator columns as numbers 27; two other
  numeric columns 73; two columns with one nonzero entry a row that are not 0s and 1s 12; three numeric
  columns 19; four 2. The last six, 222 forests, are what dec-B281 and dec-B282 refuse whatever class
  they came in as: 160 in test-forest-basis-terms.R (78 of three levels, 8 of three indicator columns,
  74 of two, three or four numeric columns), 30 in test-bcf-family.R (19 three-column blocks and 11
  scaled two-column blocks through the data door), 12 in test-forest-capture.R, 5 in
  test-forest-labels.R, 4 each in test-forest-basis-r5.R and test-multiforest-leaf-prior-writer.R, 3 in
  test-formula-terms.R, and one each in test-argument-surface.R, test-bcf-creation.R,
  test-predict-blend.R and test-sampler-residuals.R. The 91 blocks of two indicator columns are in 14
  files: test-bcf-family.R 39, test-bcf-creation.R 16, test-bcf-loglik.R 7, test-forest-basis-subset.R
  7, test-forest-defaults.R 6, test-bcf-reporting.R 4, and one or two each in eight more. 70 forests are
  held: two levels as the second forest 38, a forest with no basis 27, two numeric columns as the
  second forest 5. No logical vector reaches the data door, no block is created with a column no row
  has, and no character matrix is handed over.
- What the tests do after creation, from the same trace. `$setForestBasis` is named on 84 lines of 15
  test files and reaches the function that reads a replacement 57 times in 12. Under push 1's rule 27
  of those meet a refusal: 17 put numbers on a forest of levels, 2 of them by a formula
  (test-bcf.R 5, test-forest-arguments.R 4, test-forest-basis-r5.R 4, and one each in
  test-bcf-mutation-pins.R, test-forest-labels.R, test-multi-forest-seam.R and
  test-multiforest-leaf-prior-writer.R), and 10 give a basis to the forest created with none
  (test-forest-basis-r5.R 5, and one each in test-bcf-r5-surface.R, test-forest-arguments.R,
  test-forest-defaults.R, test-forest-sd-unit.R and test-multiforest-leaf-prior-writer.R). The 9
  factors swapped onto a forest of levels all carry the forest's levels in its order. `predict` is given
  `bases` 85 times in 6 files, always numbers: 29 of them hand indicator columns or other numbers to a
  forest of levels (test-predict-blend.R 21, test-forest-basis-terms.R 6,
  test-generics-multithreaded.R 2), and 12, all in test-forest-basis-terms.R, hand two to four columns
  to a numeric forest. Swaps of two numeric columns onto a numeric forest: 7, in test-forest-labels.R
  and test-forest-basis-r5.R.
- The trace again, step 1.0, on a15b05ca, the tree push 1 is built on, before any edit: the 73 test
  files that name a basis, 11494 results. 1444 models of several forests in 68 files (1194 of two
  forests, 236 of three, 14 of four). Per forest with a basis, 1794 in all: two levels 463, of which 5
  handed over as a value; one numeric column 1018 (a 0/1 column 374, a constant column 31, another
  613); the six shapes of a refused size 222 in 11 files, the blocks of two indicator columns 91 in 14
  files, 78 of them through the data door, and the 70 held forests, each file's count as above.
  `$setForestBasis` reaches the function that reads a replacement 58 times in 12 files, 27 of them
  against push 1's rule, file for file as above; the one more is a number given to a numeric forest in
  test-forest-sd-unit.R. `predict` is given `bases` 85 times in 6 files, 29 of them numbers for a
  forest of levels, of which 10 hand the fit's own `bases` back (8 in test-predict-blend.R, 2 in
  test-generics-multithreaded.R); one more hands back a numeric forest's. A fit's `bases` is read
  outside the package 115 times in 8 files, counted by execution and not by line
  (test-forest-basis-terms.R 48, test-argument-surface.R 32, test-bcf-loglik.R 14,
  test-predict-blend.R 12, test-bcf-forest-channel.R 3, test-forest-capture.R 3,
  test-generics-multithreaded.R 2, test-predict-forest.R 1). No count that sizes a step is over the
  earlier tree's by a fifth.
- Who reads a fit's `bases`, by search of 37b59ea0. In R/ the three functions named above. In the tests
  21 lines of 8 files: handed back to `predict` in test-predict-blend.R
  (["blendAtTraining"](../../inst/tinytest/test-predict-blend.R), whose first use holds the result to
  the stored training fit within 1e-12 on a forest created with `~ factor(z)`) and three times in
  test-generics-multithreaded.R; multiplied into the coefficients in test-bcf-forest-channel.R,
  test-bcf-loglik.R, test-argument-surface.R and, through
  [`recombine`](../../inst/common/recombine.R), test-predict-forest.R; and pinned for its dimensions,
  its column names or its values in test-bcf-forest-channel.R, test-predict-blend.R,
  test-forest-basis-terms.R and test-forest-capture.R. No benchmark script, help example, vignette or
  consumer reads one.
- Elsewhere. benchmarks/R: every basis at creation is a factor of two levels, written `~ factor(z)`
  (the four BCF exact gates, `bcf-equivalence.R`, `sbc.R`, `move-census.R`, `rows-pair.R`, and
  [`bcfForests`](../../benchmarks/R/composition-matrix.R) of the composition matrix's runner), but for
  four. Three are in the runner and multiply by one 0/1 number: its `bart()` door for its `bcf` row
  ([`doBart`](../../benchmarks/R/composition-matrix.R)), its `xbart()` recipe for that row
  (`basis = ~yBin`) and its survival check (`forest(x1, basis = ~z)`). The fourth is the scenario
  ["bart2twoforest"](../../benchmarks/R/equivalence.R), which
  multiplies by a factor of three levels, one of the 55 recorded scenarios of the gaussian baseline.
  `bcf-equivalence.R` swaps `cbind(1 - z2, z2)` onto a factor forest once, in the scenario
  ["set_treatment"](../../benchmarks/R/bcf-equivalence.R). No script hands a basis to the data door.
  The runner has no cell and no probe on a basis's kind, a swap of one or a prediction's bases, and
  [feature-matrix.md](../design/feature-matrix.md) states none. Help: the example of man/bart.Rd
  multiplies by a 0/1 integer and calls it the Bayesian causal forest shape; man/forest.Rd and
  man/bart.Rd describe `basis = a + b` as two columns of one forest, with `1 + dose`, `poly()` and
  `ns()` beside it, and man/forest.Rd says a value's unused level is refused and its emptied one kept
  as a column of zeros; the last example of man/forest.Rd, `twoColumns`, fits
  `basis = z + I(w / 50)`; man/bartBT.Rd describes a fit's `bases` as each forest's expanded matrix and
  gives `bases = cbind(1 - zstar, zstar)` as the prediction under an arm; the forests item of
  inst/NEWS.Rd says `forest(x1 + x2, basis = z)` declares the Bayesian causal forest. One vignette
  names `forest(basis = factor(z))` in prose; none fits a model of several forests.
- Consumers, read and not run. bartCause's `bcf()` (dbarts-1.0, 70ad798) hands `cbind(1 - z, z)` to
  `dbartsData(bases = )` under its own `subset`, having refused a treatment with a missing value and
  one that the kept rows leave with one value; its treatment forest takes `amplitude = fixed()` when
  `update.b = FALSE`, which on two columns as the second forest holds (0, 1); it names its coefficients
  itself, calls no swap, passes no `bases` to `predict` and reads no fit's `bases`. Two of its test
  files build the same sampler by hand from the same block (tests/testthat/test-14-bcf.R in three
  places, test-03-responseFit.R in one). stan4bart (bartcore, ac80108) refuses `forests` in its
  `bart_args` by name, makes its data object with no `bases`, and in its multivariate fitter hands
  `dbarts()` a `forests` list of one forest with no basis; treatSens (dbarts-1.0, babfaa6) and bairrtt
  (main, 3f57f61) declare one forest with `forest(n.trees = )` and no basis. The flat C header has no
  entry that takes a basis.

## The rule

At creation, at every door. The bracket is the push that makes the row so.

| handed over | kind | coefficients |
|---|---|---|
| nothing | none | one, on the forest itself |
| a factor, ordered or not; a character vector; a logical vector; `I()` of any of them: of two levels [1] | levels | one per level, in the factor's own order of levels, the sorted values of a character vector, `FALSE` then `TRUE` |
| an integer or double vector, or a matrix of one column, whatever its values; `I()` of one; one numeric term of one column [1] | numbers | one |
| a factor, a character vector or a logical vector of three or more levels | levels until push 2; then refused by name [2] | |
| a numeric matrix of several columns; several numeric terms; one term of several columns, a spline or a polynomial among them; `1 +` beside a column | numbers until push 2; then refused by name [2] | |
| a matrix of logicals or of characters, or a factor matrix, of one column or of several; a data frame; anything else | refused by name [1] | |
| a factor with NA as a level and a kept row at it | refused by name [1] | |

- The class is looked at by one function, which every reader of a basis calls. It answers levels,
  numbers, or the name of the class it does not take, and refuses nothing itself: the reader that
  called it words the refusal, a code door in the words a basis written as code is refused in, the
  others in the texts below. The data door hands it what it is given.
- The content is never looked at to decide a kind: a 0/1 number is a number, and a factor's indicator
  columns handed over as numbers are several numeric columns, refused as such.
- A basis written as code is the class its code gives in the fit's model frame; the same vector handed
  over as a value is that class too. The same column written as code and handed over as a value is one
  fit, under any `subset` and any na.action.
- The rows (dec-B308). A basis of levels is the levels that have a row among the rows the fit keeps,
  after `subset` and the na.action, in every spelling and at every door, with no message. A row of
  weight zero, or one a mask leaves out of the likelihood, is a row the fit keeps. Fewer than two levels
  left is refused, in the tip's words.
- A missing value is judged on the rows the fit keeps, after `subset` and the na.action, at every door,
  as a basis written as code is judged today and as `lm` judges a variable. One on a kept row is
  refused, in the tip's words; one on a row the fit drops is not looked at, with no message. A value
  that is not finite goes the same way. The length alone is checked against the data as handed over:
  a value has a row for each row of the data before `subset`, and loses the rows the fit drops. A
  basis takes no part in deciding which rows are kept. At `predict` the rows kept are those of
  `newdata` that the prediction keeps, as today; a swap has every row.
- NA as a level, `factor(x, exclude = NULL)` or `addNA()`, is judged on the kept rows too. With a kept
  row at it, it is refused at every door, code included, with one text that says to give the level a
  name. With none, it is dropped as any level no row has. A swap and a prediction refuse a row at it
  with the same text and pass by one that is declared and on no row.
- A numeric column that is zero on every kept row is refused at creation, at every door (dec-B346 has
  creation refuse it), the data door included.
- The size is judged after the rows: a factor of three levels that `subset` leaves with two is a basis
  of two levels, as code and as a value.
- A hold changes nothing here: a basis of three levels or of several columns is refused drawn or held,
  in the same words.
- The record. A factor's levels are kept with the data, per forest; a forest's kind is read from the
  data by one function and from nowhere else. The data object's block is the engine's form of a basis,
  a factor's indicator columns among them, and stays so.
- What a fit stores. A fit's `bases` holds, for a forest of levels, the factor: one value for each row
  the fit kept, its levels those recorded at creation. For a numeric forest it holds the block, and
  `NULL` for a forest with none. Each entry is what `predict(bases = )` takes for its forest, so
  `predict(fit, x, bases = fit$bases)` at the training rows is the stored training fit. A fit's kind
  and levels are read from that entry by the function that reads a data object's.

After creation, at `$setForestBasis` and at `predict(bases = )` with new rows, by one function, which
every spelling reaches: a value, a one-sided formula at a swap, a list and a bare value at `predict`.

| forest created as | takes | refuses |
|---|---|---|
| none | nothing | any basis (dec-B452 at a swap; refused at predict today), at both with one text |
| levels | a factor, a character vector or a logical vector, matched to the forest's levels by name (dec-B289); a level of the forest that no row has is a column of zeros | numbers, indicator columns and the sampler's own stored block included (dec-B260 at a swap, dec-B453 at predict); a row at a level the forest does not have (dec-B289, dec-B308); a row at an NA level |
| numbers | a numeric vector or a matrix of one column, a column of zeros among them: at a swap by dec-B346, at predict as today, a prediction at zero being the plain use of a 0/1 number | a factor, a character vector, a logical vector (dec-B452 at a swap); from push 2, several columns |

- Levels are matched by the values the rows carry, as `predict()` matches a factor for `lm`: a level
  the value declares and no row has is passed by, so the factor a forest was created from, an unused
  level and all, is taken again at a swap. A name matches when it is the forest's level character for
  character: "A" is not "a".
- A levels forest has the same columns for life: its width cannot change, coefficient j multiplies
  level j's rows after any swap, and the recorded levels are creation's whatever a swap declares.
- The sampler keeps the block and the levels and not the factor. A caller that will hand a basis back,
  to restore a proposal it rejected, keeps the factor it handed over; the forest's levels are read at
  `data@basis.levels[[f]]`, as its block is at `data@bases[[f]]`.
- A state carries no data (dec-B231, dec-B418). Installed in a sampler whose factor has the same levels
  in another order it goes by position, coefficient j to the recipient's level j, with no message; the
  recipient's levels and block are its own before and after.
- From push 2 no forest changes width, held or drawn: a numeric forest is one column. Between the
  pushes a drawn numeric forest takes another width as the tip does, and a held forest refuses one by
  forest-defaults-by-kind's interim text.
- A basis written as code and rebuilt at new rows from the fit's own record is untouched.
- Nothing of this moves a draw of a sampler that is swapped to a block the tip would have built the same
  way.

When a sampler is made. A data object's record is checked by one plain function: the count of a
forest's recorded levels is its block's width, and from push 2 the size is one this version takes. It
reads the levels through the accessor that answers `NULL` for an object saved without the slot, so such
an object passes. The sampler's `initialize` calls it, which is `new("dbartsSampler", control, model,
data)`, `dbarts()` and every `copy()`; `dbarts()` and `dbartsSpec()` call it on a data object they are
handed; the class's validity calls it too. Both conditions hold through every swap the table takes, so
a copy never meets them. The rows are another matter and are judged at a creation only, where
`dbartsData()` has the rows the fit keeps and where `dbartsSpec()` installs a declaration: a level or a
column that a swap left with no row is in a sampler's data for good (dec-B346), and `copy()`, a reload
and a sampler made again from that data take it as it is. A reload's re-creation is not checked at all.

`$setData` on a sampler of one forest refuses a data object that carries bases, from push 1: the
sampler has no forest to multiply them.

Accepted where refused: a factor, an ordered factor and a character vector at the data door; a factor
with a level no row has, handed over as a value, the level dropped; a number or a factor with a
missing value on a row the fit drops, as a value at a forest door and, the number under a missing
response, at the data door on the matrix interface (13 cells of Context's 28); a factor of one value
at a swap and at predict, where its level is one of the forest's (every new row treated). Changed with
no refusal: a logical vector at the data door, and a value's level that the rows empty, the first and
fourth sequences of the `rng:` line; and what a fit stores for a forest of levels. Refused where
accepted: at the data door a matrix of logicals, of one column or of several, a column of zeros, a
logical vector with one value, a basis with a row whose norm is not representable and a matrix of no
columns, the last two refused as the forest doors refuse them already; a character matrix of one column
at the forest doors, code included (a vector is levels, a matrix is numbers or nothing), and a factor
that has dimensions at the forest doors, at a swap and at predict, its one column read as levels
before; a factor with NA as a level and a row at it, as a
value at a forest door and at a swap; a value of several numeric columns of which the rows empty one,
which is refused as code today; `$setData` on a sampler of one forest handed bases; every cell of the
second table's last column, the sampler's own stored block handed back to a factor forest among them;
and, with push 2, three or more levels and several numeric columns.

Before and after, per door.

| door | before | after push 1 | after push 2 |
|---|---|---|---|
| a `forest()` term at `bart()` and `dbarts()`; a `forests` list on a formula, on a matrix, in `dbartsSpec()` | by class; a value's unused level refused and its emptied level kept as zeros, where code drops both; a value's missing value refused wherever it lies | by class; the level dropped in every spelling; a missing value judged on the rows kept | one numeric column or two levels; anything larger refused |
| `dbartsData(bases = )` | numbers, a logical vector and a matrix of logicals among them; a factor refused; an emptied column kept | by class, and by the rows as every door | the same sizes; two numeric columns, a factor's indicators among them, refused |
| `$setForestBasis` | any class, any forest, levels by position | the table | the table; no forest changes width |
| `predict(bases = )` with new rows | any class at the fit's width, levels by position | the table | the table |
| a fit's `bases` | the data object's blocks | a factor for a forest of levels | the same |
| `predict` from a stored record of a code basis; `fitted`; `extract` | | unchanged | unchanged |
| `$getLeafPrior`, `$setLeafPrior`, `print`, `show` | say nothing of a kind | unchanged | unchanged; the reader's entry and the printed block are the multiplier slice's |
| `new("dbartsSampler", control, model, data)`, `dbarts(data)`, `dbartsSpec(data)`, `copy()` | carry the block, unchecked | carry the block and its levels: the same kind; levels of another count than the block's columns refused | a data object edited to a refused size is refused |
| a reload | carries the block, unchecked | carries the block and its levels, unchecked | the same |
| `$setData` on a sampler of one forest, the data carrying bases | taken; the plain model runs and `copy()` fails | refused | refused |
| a state stored or installed | carries no basis; coefficients by position | unchanged | unchanged |
| the engine | a coefficient for each column of any block | unchanged | unchanged: nothing is narrowed there |

## Refused forms, with their texts

Base R's style, as in written-surface, a value in single quotes as the tip's texts have it. `<arg>` is
`basis` at a forest door and at `$setForestBasis`, `bases` at the data door and at `predict`. It stands
in a text as the quoted word only and never inside code: at the data door and at `predict` it is a
list, and [`dbarts`](../../R/dbarts.R) restates the data object's quoted word and nothing else when it
hands a refusal on. `<f>` is `forest 2`, or `forest 2 ('dose')` where the forest's label
([`forestLabels`](../../R/forestBasis.R)) is not its position's own; the data door has no labels and
says `forest 2`. Every new text names the forest and says what to do at the door that raised it, in
words wherever the door has no code of the caller's to show. The push that adds each is in brackets.

    [1] '<arg>' gives <f> a logical matrix; a basis is a factor, a character vector, a logical vector, or numbers. Give its column as a logical vector for two levels, or as numbers for one coefficient            (a character matrix: "Give its column as a character vector"; a factor that has dimensions, named "a factor matrix": "Give its column as a factor"; a matrix of several columns: "Give one of its columns as ...")
    [1] '<arg>' gives <f> a data frame; a basis is a factor, a character vector, a logical vector, or numbers. Give one of its columns            (a list the same; a Date, a complex: named by its class, the first sentence alone)
    [1] 'basis' (cm) has a term that is a character matrix: a basis is a factor, a character or logical vector, or numeric; write its column as a vector, as cm[, 1]            (a logical matrix and a factor matrix of one column the same; of several columns, the tip's text as it is)
    [1] '<arg>' gives <f> a factor with NA as a level, and a level of a basis has a name: give that level one            (written as code it begins "'basis' (B) is a factor with NA as a level")
    [1] '<arg>' gives <f> a basis, and <f> was created without one: what a forest multiplies is fixed at creation. Create the model again, with the basis, to give it one
    [1] <f> multiplies a factor with the levels '0', '1', and '<arg>' gives it numbers. Give a factor, a character vector or a logical vector with those levels; if the numbers are a factor's indicator columns, give the factor itself. To multiply numbers, create the model again with them
    [1] <f> multiplies numbers, and '<arg>' gives it a logical vector. Give numbers: as.numeric() of a logical vector is its 0/1 column. To multiply two levels, create the model again with the vector
    [1] <f> multiplies numbers, and '<arg>' gives it a factor. Give the column the forest was created with, at these rows; to multiply a factor's levels, create the model again with the factor            (a character vector the same)
    [1] '<arg>' gives <f> the level 's', which it does not have: its levels are 'p', 'q', fixed at creation. Give every row one of them
    [1] 'basis.levels' gives forest 2 3 levels and its block in 'bases' has 2 columns; dbartsData() and dbarts() write the two together            (a count of one is "1 level", "1 column")
    [1] $setData: the data carry forest bases, and this sampler was created with one forest and none; a sampler's forests are fixed at creation. Make a new sampler from the data with dbarts(), or hand over data without 'bases'
    [2] <f> has a basis of 3 levels ('a', 'b', 'c'); in this version a factor, a character or a logical basis has two levels. Give each level but the first a forest of its own, as forest(x1, basis = arm == "b") + forest(x1, basis = arm == "c")            (a value: "Give each level but the first a forest of its own, its basis a logical vector that is TRUE on that level's rows")
    [2] <f> has a basis of 2 numeric columns ('dose', 'age'); in this version a basis is one numeric column, or a factor, a character or a logical vector of two levels. Give each column a forest of its own, as forest(x1, basis = dose) + forest(x1, basis = age)            (a value, or one term of several columns: "If the columns are a factor's indicators, give the factor itself; otherwise give each column a forest of its own")
    [2] <f> has a basis of 2 numeric columns ('(Intercept)', 'dose'): '1 +' in a basis adds a constant column, and in this version a basis is one numeric column. Write the column alone, as basis = dose
    [2] <f> multiplies one numeric column, and '<arg>' gives it 2; in this version a basis is one numeric column. Give the one column

What each advises runs where it is raised. All but the third and the first three [2] texts say it in
words, a value door having no code of the caller's and `<arg>` being a list at two of them. No text
advises `as.numeric()` of a factor, which gives its codes: the seventh names it for a logical vector
alone. The sixth names the forest's levels where it asks for a factor, and it and the second [2] text
tell the caller whose numbers are a factor's indicator columns, `cbind(1 - z, z)`, to hand the factor
over. The third is the code door's: the tip's text for a term that is a logical matrix, the class it
names being the one function's answer, with the column written in the caller's own words, which the
grammar takes (`forest(basis = cm[, 1])`, run). The fourth is one text at every door; a code door
begins it as it begins every text, with the basis as written. The fifth is the table's for a forest
with none, at a swap and at `predict` alike. In the first two [2] texts the example is code only where
the basis is: the forest's predictors as written, the factor's name, each term that is one column; a
quoted level inside it is R's. Where the basis is a value, at a forest door or the data door, or one
term of several columns, the words in brackets stand in its place and no `forest()` call is shown. The
third [2] text is chosen by how the basis was written, `1 +` in its code, and never by a column of
ones: `cbind(1, dose)` handed over takes the second. The last [2] text is the swap's and the
prediction's.

The tenth is the plain function's, raised by the sampler, by `dbarts()` and `dbartsSpec()`, and by the
data object's validity. Replaced by push 1: "'bases' must be numeric or logical" and "a 'basis' must be
numeric, a factor, or a character vector"; the data door's "'bases' values must all be finite", by the
forest doors' "a 'basis' cannot be NA" and "a 'basis' must be finite" with `<arg>`, on the rows kept;
at `predict`, the text for a forest given a basis it has none of, by the fifth; and, for a factor with
NA as a level written as code, "has the level 'NA' at a new row, which no row of the fit had", by the
fourth. Gone with push 1, its case at creation no longer a refusal (dec-B308): "a 'basis' factor level
with no observations contributes nothing to a forest: ...; drop it with droplevels()". Replaced by
push 2: [`refuseHeldShape`](../../R/model.R)'s text for a held basis of three or more columns and
[`refuseHeldWidthChange`](../../R/model.R)'s for a held forest's width at a swap; "'bases' gives forest
2 2 columns; its amplitudes take 1" where the forest is numeric;
[`refuseReorderedBasisNames`](../../R/forestBasis.R)'s refusal of a swap whose columns stand in another
order, which has no case left. Kept in the tip's words and extended to the doors that lack them: "a
'basis' column of all zeros contributes nothing to a forest", "a 'basis' factor must have at least two
levels" and the text for a logical with one value, each now met at the data door too, with `<arg>`;
so are "a 'basis' row's norm is not representable" and "a 'basis' must have at least one column", which
take `<arg>` with them.
Kept as they are: the texts about a basis's length, the at-least-two-forests text,
[`refuseHeldOneColumn`](../../R/model.R)'s refusal of a held single numeric column (dec-B353),
forest-defaults-by-kind's two interim texts by position, at `predict` the refusal of `bases` with no
new rows, and every other text of a basis written as code (a logical matrix of several columns as a
term, a factor beside other terms, a level at a new row that the fit did not have, an operator the
grammar does not take), which comes first where both apply.

## Constraints

- Every fit that both builds accept, the four sequences of the `rng:` line apart, draws what it draws
  on the base build. The gaussian baseline is re-recorded for the one scenario push 2 rewrites and for
  nothing else; no snapshot file is regenerated.
- No change to src/, to the flat C header, to the stored state or to what the bridge reads. The engine
  keeps a coefficient for each column of whatever block it is handed, and is not told what 1.0-0
  refuses. The forests' record on the control keeps `hasBasis`, `labels` and `basisTerms` as the sd unit
  and written-surface leave them.
- The class is read in one function at creation and the table applied by one function afterwards; the
  size is refused by one function; no door has a copy of any of the three. A caller of the first words
  its own refusal and tests no class.
- `data@bases[[f]]` stays a bare numeric matrix, with the names it has today and no other attribute.
- Of a fit, only the `bases` entry of a forest of levels changes; `forestFits`, `glue`, `basis.terms`,
  the labels and every draw are what they are, and `extract` and `predict` return what they return.
- The rows are judged at a creation and nowhere else: no check of a level's or a column's rows stands
  where `copy()`, a reload or a sampler made again from its own data would run it.
- No prior moves: until the multiplier slice a factor of two levels and a 0/1 number keep the fits they
  have.
- What written-surface fixed for later stays fixed (dec-B272): how the columns of a basis of several
  terms are named, and every refusal of a spelling that could later state a size for each column. A
  basis of several terms is still read and named, and then refused.
- The reader, `extract` and `print` are untouched. The unit and the default of an sd are untouched.
- Nothing of `updateBasisScale` rides along.
- Push 1 adds no test on a factor of three or more levels or on several numeric columns, the rows of
  its loop apart that push 2 turns.
- Each push leaves the help saying what the code does. Base R calls stay within DESCRIPTION's R floor.

## Pushes

1. The rule. Changes the four sequences; the design note lands here. After it the data door takes a
   factor, which is what lets bartCause hand one over.
2. One size a forest. Two commits gated together, and the baseline in a commit of its own after review.
   - First: every call in the tests, the benchmarks and the help that hands a two-level factor's
     indicator columns over as numbers is written as the factor it stands for, bit for bit; the one
     baseline scenario on three levels is rewritten.
   - Second: a basis of three or more levels, and of several numeric columns, is refused at every door
     (dec-B281, dec-B282), and the tests of those shapes are rewritten or go.
   - The baseline: the rewritten scenario recorded and merged, and the gates repinned on the new file.
   - The three go up in one push. cpp-tests' bitwise compare fails on the first two by design, the
     scenario having moved and the baseline not yet; a per-push run is of the push's last commit.

bartCause's line lands on its own branch after push 1 and before push 2: its `cbind(1 - z, z)` through
the data door is two numeric columns, which push 2's second commit refuses. No tip breaks it: push 1's
tip takes both spellings, and push 2 is not pushed until bartCause's suite has passed on it with the
factor.

Push 1 alone leaves a coherent tip: its own respelling is the swaps and predictions it refuses, 27 and
29 at the traced commit less the predictions that hand a fit's own `bases` back, which it carries. The
two commits of push 2 are cut where the review changes: the first is read against "the same draws" and
nothing else, the second against "which door did it miss" and "what coverage went with a deleted
test". The pushes are neither split further nor put in another order: the respelling can
be proved only on a tree that takes both spellings, which push 1 makes of the data door, and the size is
refused at every door at once.

## Steps

"Fails today" is what the base build does where the test expects otherwise. New names are for the
reader's sake. The fixture is Context's, with `arm` a factor of three levels beside g. A value is handed
over through `do.call(forest, list(basis = v))`, or directly at the data door, at `$setForestBasis` and
at `predict`: a name written in a call of `forest()` is code.

### Push 1: the rule

1.0 The trace, again. Before any edit, on the tree the push is built on, Context's trace is run again:
    each test file that names a basis in its own process, the functions that read a basis wrapped in
    the loaded namespace, nothing installed. Counted per file: the creations by door, class, kind and
    width; each `$setForestBasis` that reaches [`replacementForestBasis`](../../R/forestBasis.R), with
    its forest's kind and its value's class; each `predict(bases = )` likewise, and whether it hands the
    fit's own `bases` back; each read of a fit's `bases`. The counts go in the landing note beside
    Context's, and every count in the steps below, and every stop that names one, is read as the
    rerun's. A count that sizes a step (the creations of a refused size, the blocks of two indicator
    columns, the swaps and the predictions the rule refuses) more than a fifth over Context's is
    reported before anything is built.
1.1 One reader of a class. The class of a basis is tested in five places today, each written out:
    [`expandForestBasis`](../../R/model.R) (a value, and the factor a code basis hands it),
    [`basisColumns`](../../R/forestBasis.R) (the variables of a code basis's model frame),
    retired: [`refuseEmptiedValueBasis`](../../R/forestBasis.R) (a value against the rows kept; gone
    since step 1.3 was built),
    [`replacementForestBasis`](../../R/forestBasis.R) (a formula given to a swap) and
    [`validateForestBases`](../../R/data.R) (the data door's "numeric or logical"). One function
    (`basisClass`) answers levels, numbers, or the name of a class it does not take ("a logical
    matrix", "a character matrix", "a data frame", "a Date"), and refuses nothing. The five call it and
    test no class themselves; each words its own refusal from the answer.
    [`basisColumns`](../../R/forestBasis.R) uses the tip's text for a term, the third of "Refused
    forms", so a character matrix of one column written as code is refused there; the readers of a
    value use the first two. The data door reads a value as a forest door does, and reads it where it
    has the rows the fit keeps (step 1.3): [`validateForestBases`](../../R/data.R) checks the length
    against the data as handed over and carries the value, cut with the rows
    ([`restrictBasesToRows`](../../R/data.R), [`alignForestBasisToSubset`](../../R/model.R)), to that
    point, so the data door has no rule of its own. `dbarts()` and `dbartsSpec()` hand a value on to
    the data object as it was given and no longer expand it first
    (retired: [`expandValueBasis`](../../R/forestBasis.R), gone since), so that the data object's reading is the one reading
    of a value. A code basis's column of row numbers passes as it does. Tests, a new file
    test-forest-kind.R, built as one loop that asserts against the FIRST door, so that no door can be
    right alone.
    - The objects. Each object of Context's first table is taken through the seven creation doors. Per
      object: the kind, the width, the levels, `unname(data@bases[[2]])` `identical()` across doors, and
      the 10 seeded draws `identical()` across the sampler doors. A refused object is refused at every
      door a value reaches with one text but for `<arg>` and `<f>`, and at a code door with the code
      door's. Fails today: at the data door for a factor, an ordered factor, a character vector and
      a logical vector, which it refuses or reads as a number, and for a matrix of logicals and
      `rep(0, n)`, which it takes; as a value for a factor with an unused level, which is refused, and
      for one with NA as a level, which is created; as code for the character matrix of one column,
      which is created, and for the NA level, which takes another text.
    - The rows a missing value is judged on. A number and a factor with a missing value on a row that
      `subset` drops, and on a row that a missing response drops, at the seven doors such a row
      reaches: the three code doors, a value in a list on a formula and on the matrix interface, the
      data door on a formula and on the matrix interface. Each is created, block and draws
      `identical()` to the first door's, which is the fit of the complete column (fails today in 13 of
      the 28 cells). The same two with the missing value on a row the fit keeps are refused at all
      seven with one text. `dbartsSpec()` over a data object made under the `subset` creates from the
      value at the kept rows and refuses the value at the data's length for its length, as today.
    - `addNA()` of a complete factor is created at every door with the factor's two levels (fails
      today as a value and at the data door).
    - The loop's table has one column for what push 2 expects, and the objects of three levels and of
      several columns sit in it as created here.
    - Pinned beside it, so that a reading by content fails: a 0/1 integer, a 0/1 double, `I(z)` and
      `matrix(z)` are numbers of one column and do NOT draw what `factor(z)` draws; `cbind(1 - z, z)`
      handed over as a value is numbers of two columns. The same column written as code and handed
      over as a value is one kind and one fit.
1.2 The record. A slot `basis.levels` on [`dbartsData`](../../R/A_class.R), beside `response.levels`:
    `NULL`, or a list with one entry per forest, `NULL` for a forest with no basis or with numbers and
    the levels for one with a factor. It is written wherever `bases` is: where
    [`dbartsData`](../../R/data.R) finishes its rows, with the levels following every restriction of
    rows and never left on the stored block; where `dbarts()` puts the blocks of its code bases on the
    data object, from the levels [`buildCodeBasis`](../../R/forestBasis.R) returns in its record; and
    where a declaration's bases are put on a data object ([`resolveSamplerSpec`](../../R/spec.R)),
    replacing the object's blocks and levels whole when any forest of the declaration gives a basis,
    as the base build replaces the blocks: a list that states settings and no basis leaves the
    object's blocks and their levels as they are, a numeric declaration for a forest that had levels
    leaves it `NULL`, and a forest the declaration gives no basis has none afterwards, whatever the
    object carried at its place. The slot is read through
    an accessor (`dataBasisLevels`) that answers `NULL` for an object saved without it, as
    [`dataRowNames`](../../R/data.R) and [`dataMissingSeen`](../../R/data.R) do. One function,
    `forestBasisKind`, answers "none", "levels" or "numeric" for a forest of a data object and for an
    entry of a fit, and is the only reader of a kind after creation. The control's `basisTerms` keeps
    its own copy of a code basis's levels for the stored record, written in the same call.
    - The check. One plain function (`checkBasisRecord`, of a data object) refuses recorded levels
      whose count is not their block's width, with the tenth text, reading the slot through the
      accessor. Its callers: the sampler's `initialize`; `dbarts()` and `dbartsSpec()` where they are
      handed a data object; and the class's validity. No creation door calls `validObject()` on a data
      object: R's own validity refuses an object saved without a slot, which must go on creating. A
      reload's re-creation ([`recreatePointer`](../../R/dbarts.R)) calls nothing and is not checked.
    - The fit. [`packageBartResults`](../../R/bart.R) stores for a forest of levels the factor, made
      from the block and the recorded levels, and the block otherwise. One function
      (`storedBasisBlock`) gives the block of an entry back, a factor expanded on its own levels, and
      the three readers of R/generics.R go through it: [`predict.bart`](../../R/generics.R) for the row
      it hands on when no row of `newdata` is left, [`extractForest`](../../R/generics.R) for a
      forest's contribution, [`resolveForestBases`](../../R/generics.R) for a forest's kind, levels
      and width. A fit saved before holds blocks and reads as numbers.
    - Tests, the slot: at each of the seven doors for a factor, a character and a logical vector, a
      number and no basis; for a code basis `identical()` to the levels of the control's record; the
      same under `subset` at the formula, list, matrix and data doors, and after a dropped missing
      response (a record lost with the rows makes the forest numeric there);
      `attributes(data@bases[[2]])` holds the dimensions and names and nothing else.
    - Tests, `dbartsSpec()` over a data object, both ways: one that carried a numeric block for the
      forest records the declaration's levels; one that carried levels keeps them, and its block,
      under a list that states `n.trees` and no basis, and has `NULL` under a numeric declaration.
    - Tests, by hand, each at a named door. An object with the slot removed
      (`attr(data, "basis.levels") <- NULL`) creates at `new("dbartsSampler", control, model, data)`,
      at `dbarts(data)` and at `dbartsSpec(data)`, runs, and reads as numbers; `validObject()` refuses
      it, in R's words, and none of the three calls it. An object whose levels are edited to another
      count is refused with the tenth text at each of the three and by `validObject(data)`. No test
      is written on a reload of an edited object: it is not checked, and nothing is promised of it.
    - Tests, the fit: `fit$bases[[2]]` of a fit on a factor, on a character vector and on a logical
      vector is a factor over the fit's rows with the recorded levels, from a `forest()` term, from a
      `forests` list and from `bart(data)` over an object made at the data door (fails today: a
      matrix); under a `subset` that empties a level its levels are the two kept; for a number it is
      the block of today; `storedBasisBlock` of each entry is `identical()` to the kept sampler's
      `unname(data@bases[[2]])`, and `extract(fit, type = "forest", contribution = TRUE)` equals the
      product computed in the test from that block.
1.3 The rows at creation. One function (`finishBasisRows`) finishes a basis on the rows the fit keeps,
    in this order. A missing value, or one that is not finite, on a kept row is refused in the forest
    doors' words. A factor with a kept row at an NA level is refused with the fourth text. Of a levels
    block each column no kept row has goes, with its level (dec-B308), an NA level among them, and
    fewer than two left is refused in the tip's words, a factor's and a logical's. Of a numeric block
    a column of all zeros is refused in the tip's words. It is called where the data object has its
    final rows, after `subset` and the na.action, on both interfaces, and where `dbartsSpec()` installs
    a declaration. A code basis is built on the kept rows already: its frame's check of a missing
    value is the same rule and stays, and its factor's NA level takes the fourth text there. The
    function is called nowhere else: not from the sampler's `initialize`, and so not from `copy()`, and
    not from [`recreatePointer`](../../R/dbarts.R). retired: [`refuseEmptiedValueBasis`](../../R/forestBasis.R),
    gone since, folds into it, and at creation [`expandForestBasis`](../../R/model.R) no longer refuses an unused
    level or looks for a missing value: it is handed the kept rows. Tests, each at the doors a value
    reaches (a list on a formula, the matrix interface, `dbartsSpec()`, the data door), with one text
    or one block: a factor value with a level never used is created with two levels (fails today:
    refused); a factor of three levels of which `subset` empties one, and one emptied by a missing
    response, are two columns, block and draws `identical()` to the same factor written as code (fails
    today: three columns, one of zeros, and other draws); a factor of two levels left with one is
    refused (at the data door today as not numeric); a logical left with one value is refused with the
    logical's text (fails today: a constant column at the data door, the droplevels text as a value
    with no `subset`); a number left all zero by `subset`, and `rep(0, n)`, are refused (fail today at
    the data door); a factor whose rows at one level all have weight zero keeps both levels at every
    door; a factor with a kept row at its NA level is refused at every door, code included, with the
    fourth text, and with those rows dropped by `subset` is created with its two named levels. The
    missing value's rows are step 1.1's.
1.4 After creation. One function (`conformBasis`: a value, the forest's kind, its levels, its width,
    its label) returns the block for the engine or refuses by "The rule";
    [`setForestBasis`](../../R/dbarts.R) and [`resolveForestBases`](../../R/generics.R) both call it,
    the sampler reading kind and levels from its data object and a fit from its stored entry. In
    `$setForestBasis` nothing is stored before it returns, and the recorded levels are never written.
    A one-sided formula is read first, as today, and what it gives is conformed as a value is: the
    factor of its model frame, or its numeric columns. At `predict` a list's entry and a bare value
    both reach it. In this push it adds no rule about a numeric forest's width: a drawn one takes what
    the tip takes, and a held forest's width is forest-defaults-by-kind's interim refusal. Tests,
    every factor in them of two levels:
    - The table: each kind of forest (none; a factor, a character and a logical forest; a 0/1 number;
      dose) by each kind of value, at `$setForestBasis` as a value and as a one-sided formula
      (`~ dose`, `~ factor(z)`), and at `predict` in a list and as the bare value. Each cell is
      accepted or refused with its text (fails today: at a swap every kind of value but a factor of
      one level and a logical matrix is taken by every drawn forest, `~ dose` by a factor forest among
      them). A column of zeros is taken by a numeric forest at a swap and at `predict`, where the
      prediction is the fit's with that forest left out. A forest with none refuses with the fifth text
      at both. After every refused swap the data object, `getForestAmplitudes()` and the next 5 sweeps
      are `identical()` to an untouched twin's.
    - Which model, by name. A held two-level character forest swapped to the same vector with its
      levels reversed keeps its term on the treated rows (fails today: it moves to the control rows). A
      drawn two-level forest swapped to a vector left with one value keeps 2 columns, one of them zero,
      and its coefficients; swapped back, the two coefficients are `identical()` to those before (fails
      today: refused). A row at a level the forest lacks is refused; the same value with that level
      declared and on no row is taken, and `data@basis.levels` is `identical()` to creation's
      afterwards. "A", "a " and "tr" for a forest whose levels are "a" and "treated" are each refused
      with the ninth text (fail today: taken by position). A row at an NA level is refused with the
      fourth text (fails today: a third column). The factor a forest was created from, a level no row
      has among its levels, is taken at a swap and gives the block of creation. The sampler's own
      `data@bases[[2]]` handed back to a factor forest is refused with the sixth text (fails today:
      taken).
    - `predict`: a factor with its levels reversed predicts `identical()` to the aligned one (fails
      today: 1.0739 against 1.2976); every row at one level, written `rep("1", n)`, predicts what
      `factor(rep(1, n), levels = 0:1)` does today; a renamed level is refused (fails today: taken);
      numbers on a factor forest and a factor on a numeric forest are refused (fail today: taken where
      the width fits). A fit made by `bart(data)` over an object from the data door predicts by name
      as one from a `forest()` term does.
    - The training rows, the one identity here that a build cannot satisfy by being wrong the same way
      twice. On a fit created with a factor whose levels are not in sorted order,
      `factor(zc, levels = c("treated", "control"))`, `predict(fit, x, bases = fit$bases)` at the
      training rows is the stored training fit within 1e-12, the bar the first oracle of
      test-predict-blend.R holds the same call to, the blend and the engine summing in another order;
      so is `predict` given the creation factor, and given the same vector as a character.
    - With the forest with no basis second in a list, as forest-defaults-by-kind accepts it, that
      forest refuses a basis and the first takes one: the kind is not a position.
1.5 The same sampler again. Tests: after `copy()`, after `storeState()`, `saveRDS`, `readRDS` and a
    first use, and on `new("dbartsSampler", control, model, data)` from the sampler's own three, the
    kind and the levels are those of the original, the by-name swap of step 1.4 gives the same block,
    and a refused swap is refused. After a swap that leaves a level of a factor forest with no row, and
    after a swap of a numeric forest to a column of zeros, each of the three is made, holds the swapped
    block and runs 5 sweeps (dec-B346; passes today, and is what a check of the rows in the wrong place
    breaks). A state: two samplers on one character vector, the first created with
    `factor(v, levels = c("control", "treated"))` and the second with the levels reversed; the first
    runs and stores, the second installs the state. The install returns what it returns today; the
    second's coefficients are then the first's, column for column; and its `data@basis.levels` and
    block are `identical()` to its own before the install. A build that matches the coefficients to
    the level names, hands on the donor's levels, or refuses, fails it.
1.6 Respell what the rule refuses, by step 1.0's trace; run the suite first and repair what it shows.
    - Swaps of numbers onto a forest of levels, 17 at the traced commit. A block of indicator columns
      becomes the factor (`cbind(1 - z2, z2)` becomes `factor(z2)`). A test that needs numbers moves
      to a numeric forest of ONE column wherever one column serves, so that it stands through push 2:
      a scaled block, `4 * cbind(1 - z, z)` and `7 * zBasis` in test-forest-basis-r5.R, becomes a
      forest created on `z` and swapped to `4 * z`; a column of zeros, a forest on `z` swapped to
      `rep(0, n)`. A swap of another width cannot: it moves to a numeric forest of several columns,
      stands on this push's tip, and goes in step 2.7 with the refusal. No creation of two indicator
      columns as numbers is added for any of them.
    - Predictions, 29 at the traced commit. Those that hand the fit's own `bases` back stand as
      written, the stored entry being the factor:
      ["blendAtTraining"](../../inst/tinytest/test-predict-blend.R) and the calls of
      test-generics-multithreaded.R. The rest hand the factor for the indicator columns.
    - Reads of a fit's `bases`, 21 lines in 8 files. Arithmetic on an entry goes through
      `storedBasisBlock`: [`recombine`](../../inst/common/recombine.R) in its one place, which serves
      test-predict-forest.R, and the loops of test-bcf-forest-channel.R, test-bcf-loglik.R and
      test-argument-surface.R. A pin of a factor forest's entry by its dimensions becomes one of the
      factor's length and levels (test-bcf-forest-channel.R, test-predict-blend.R); a pin of its
      values or its column names reads the block through the function, or the kept sampler's
      `data@bases`, which has the names (test-forest-basis-terms.R). Pins on a numeric forest's entry
      stand (test-forest-capture.R, test-predict-blend.R).
    - The 10 swaps that give a basis to the forest created with none become a pin of the fifth text
      where the test is about that, and otherwise the forest is created with the basis. Among them
      the swap of `dose` onto forest 1 in test-forest-sd-unit.R: what it pinned, the prior of a forest
      created with none and given a column, has no case left from R, and the test says so.
    - The three pins of [`refuseHeldWidthChange`](../../R/model.R)'s text in test-forest-arguments.R
      are among the swaps the rule refuses, two of them numbers on a held factor forest and one a
      column on a held forest with none: each meets the rule first and takes its text, and the
      function is left the one case of two numeric columns held as the second forest until push 2.
    - The pin of
      retired: ["factor level with no observations contributes nothing"](../../inst/tinytest/test-forest-basis-r5.R),
      gone since, becomes a creation with two levels. A pin of the data door's "must all be finite" or of "must
      be numeric or logical" takes the text that replaces it. `bcf-equivalence.R`'s swap becomes
      `factor(z2)`; its scenario is unchanged to the bit.
1.7 `$setData` on a sampler of one forest. [`bartcoreSamplerSetData`](../../R/bartcore.R) refuses a
    data object whose `bases` is not `NULL` where the sampler carries no amplitudes, with the eleventh
    text, beside its refusal of a sampler that carries them
    ([`refuseAmplitudeMutation`](../../R/bartcore.R)) and before anything is stored, and runs the
    record's check on the object there, so that one whose `bases` were taken off by hand and whose
    levels stayed is refused as a sampler made from it is. Test: the text;
    afterwards the sampler's data, its next 5 sweeps and a `copy()` of it are `identical()` to an
    untouched twin's (fails today: taken, and the copy fails); the same data object without its
    `bases` is taken, as today.
1.8 Help and records. man/forest.Rd, the `basis` item: which class is which kind; that a 0/1 number is
    one column with one coefficient and not two levels, and `factor(z)` is how levels are written; that
    a level no kept row has gives no column, as code and as a value; that a missing value counts on
    the rows the fit keeps and on no other; that a level is not NA; its `sd` item loses the sentence on
    a forest created with no basis and given one. man/dbartsData.Rd, `bases`: the same classes as
    `forest(basis = )`, and the same rows.
    [`dbartsSampler$setForestBasis`](../../man/dbartsSampler-class.Rd) and the method's docstring: the
    table; that levels are matched by name; that the sampler keeps the block and the levels and not
    the factor, so a caller who will hand a basis back keeps the factor it handed over, and that the
    levels are read at `data@basis.levels[[f]]`; and one sentence, that a state carries no data and
    goes by position into a sampler whose factor has the same levels in another order.
    [`dbartsSampler$setData`](../../man/dbartsSampler-class.Rd): the refusal of step 1.7.
    man/bartBT.Rd: the `bases` item of `predict` takes the table, and its Bayesian causal forest
    example hands the arm as a factor with the fit's levels where it has `cbind(1 - zstar, zstar)`;
    the `bases` entry of the fit's Value says a factor for a forest of levels, that it is what
    `predict(bases = )` takes back, and how its block is had. man/na.keepPredictors.Rd says any
    missing value in `bases` is an error at `predict`; run, one on a row the na.action drops is passed
    by, and the page is made to say so. docs/design/forest-kind-by-class.md with its index row: the
    rule, the record, what a fit stores, the four changed sequences with their oracles, what a sampler
    saved before reads as. docs/design/public-surface.md and docs/architecture.md where they describe
    a basis, the data object's slots and a fit's elements. TODO: `forest-prior-args` names this push
    landed, and `written-surface-leftovers` loses its two sentences on a swap's reordered levels and
    on a zero column at the data door.
1.9 Mutations (Verification): apply each, install, run the named test, record the failing count, revert,
    `touch` the file.

### Between the pushes: bartCause

B.1 Its own commit on dbarts-1.0, after push 1 has landed on bartcore and after the sd unit's edit
    there, which touches the same two test files. R/bcf.R: `basis <- cbind(1 - z, z)` becomes
    `basis <- factor(as.integer(z), levels = 0:1)`; `as.integer` because `bcf()` takes a logical
    treatment, which `factor(z, levels = 0:1)` would turn into missing values. The comment above it
    loses its sentence about the column order of `cbind`, the order now being that of the levels 0, 1.
    Its hand-built samplers hand the factor over in the same commit, so that they stay the comparator of
    `bcf()` under the multiplier slice: tests/testthat/test-14-bcf.R in its three places (the builder's
    `bases = `, the aligned block, and the block already cut to the kept rows, which keeps its refusal
    for its length) and test-03-responseFit.R in its one. Its pin of the second column of the aligned
    block against the kept treatment stands: a factor's block is still stored as its indicator columns.
    No other file there reads a basis.
B.2 Gate: bartCause's suite whole, on push 1's tip, before and after the commit, none failing; seeded
    `bcf()` fits identical before and after (the default, a binary and a logistic response, `subset`, a
    logical treatment, weights with an offset, stated sds, each coefficient held and both). Not on any
    earlier tip: before push 1 the data door refuses a factor.

### Push 2: one size a forest

First commit, the blocks written by hand.

2.1 Tests. Each creation from a two-level factor's indicator columns handed over as numbers (91 blocks
    in 14 files at the traced commit, 78 through the data door and 13 at a forest door) is read at its
    fixture. Where a block of two columns stands for a factor (the treatment blocks of
    test-bcf-family.R, test-bcf-creation.R, test-bcf-loglik.R, test-bcf-reporting.R,
    test-forest-basis-subset.R and the others the trace names) the fixture hands the factor over, and
    a pin of `data@bases` against the block is kept, against the block. Not touched in this commit: a block of three columns, a block
    times a constant and a constant column; the creations on a 0/1 number, which are one numeric column
    before and after; and the blocks that the second commit turns into a refusal's test, those of
    test-forest-defaults.R held second among them. A swap or a prediction that gives numbers to a
    forest this step respells as a factor is respelled with it, the factor at those rows, in the same
    edit; the tests step 1.6 moved to a numeric forest of one column are not touched. Every edited
    file's results are those of the base build, assertion for assertion.
2.2 Benchmarks and help. The composition matrix runner's `bart()` door for its `bcf` row
    ([`doBart`](../../benchmarks/R/composition-matrix.R)) is written `basis = ~ factor(z)`, as its
    sampler doors are; man/bart.Rd's example of the Bayesian causal forest shape is written with
    `factor(z)`, and the line of its section "Formula Terms" that glosses a numeric z says one
    coefficient, not a treatment; the forests item of inst/NEWS.Rd, which gives
    `forest(x1 + x2, basis = z)` as the Bayesian causal forest, is written with `factor(z)`. These
    change the model the text fits, which nothing pins: the runner checks that a cell constructs or
    refuses, and no draw. The runner's other two bases on a 0/1 number, its `xbart()` recipe's and its
    survival check's, are one numeric column before and after and are left.
2.3 The baseline scenario. ["bart2twoforest"](../../benchmarks/R/equivalence.R) draws its factor from
    two values where it drew from three, with nothing else of its text changed, and its comments say
    two levels; its response moves with the factor, the scenario's data being drawn after it from one
    seed of its own, and no other scenario's data does. The file is made in its own commit after
    review, as d5b637ee made its: the scenario is recorded alone on a `--preclean` reference install at
    the push's newest package-code commit, which names the file
    (`EQUIVALENCE_SCENARIOS=bart2twoforest`, `EQUIVALENCE_CORES=1`), and merged into a copy of
    `equivalence-308c45be.rds` in its scenario order; .github/workflows/cpp-tests.yaml,
    .github/workflows/equivalence.yaml, benchmarks/R/mutation-battery.R and
    [feature-matrix.md](../design/feature-matrix.md) are repinned on it. The
    [MANIFEST](../../benchmarks/baselines/MANIFEST) row says what moved and why (the scenario's own data
    and model, not the sampler: the other 54 are bitwise on the same build), and names as oracle what
    gates a two-level factor forest, the BCF exact gates in `quick` on this build.

Second commit, the refusal.

2.4 One function refuses a size (`refuseBasisSize`: the built block, its levels or none, the forest's
    label, the basis's text and its terms where it is code, the forest's own predictors as written).
    More than two levels: the first [2] text. More than one numeric column: the second, or the third
    where the basis's code has `1 +`, read from its terms and never from a column of ones. It is called wherever a creation's block is
    final: after step 1.3's function where the data object has its rows and where `dbartsSpec()`
    installs a declaration, so that a level the rows empty is dropped first; where
    [`buildFitBases`](../../R/forestBasis.R) puts a code basis's block, built on the rows kept, in its
    place; and in step 1.2's plain function beside the count of the levels, on the recorded levels
    and the block's width and never on the rows, so that the sampler's `initialize`, `dbarts()` and
    `dbartsSpec()` handed a data object, and the class's validity, refuse an object edited by hand,
    in the words for a value. Every sampler this version makes passes it for life, no forest changing
    width from this commit on, so a `copy()` passes it; a reload's re-creation is not checked. It is
    called before anything is read of a hold, so [`refuseHeldShape`](../../R/model.R)'s refusal of a
    held basis of three or more columns can no longer be reached and goes.
2.5 After creation. `conformBasis` refuses a numeric value of several columns for a numeric forest, with
    the last [2] text, at `$setForestBasis` and at `predict`; with it no forest changes width, and
    [`refuseHeldWidthChange`](../../R/model.R) goes, as does
    [`refuseReorderedBasisNames`](../../R/forestBasis.R).
    [`validateForestSd`](../../R/model.R)'s text for a longer `sd` speaks of one sd and of rescaling a
    column in the basis, and stands as it is.
2.6 Tests, a new file test-forest-basis-one-size.R.
    - The loop. Each refused object: a factor, an ordered factor and a character vector of three
      levels; `dose + age`; `ns(dose, 3)`; `poly(dose, 2)`; `1 + dose`; a matrix of two numeric columns
      held in a variable and named in code; and, as values, a matrix of two numeric columns,
      `cbind(1 - z, z)`, `cbind(1, dose)`, `model.matrix(~ g - 1)` and `model.matrix(~ g)[, -1]`.
      `cbind(1, dose)` takes the second [2] text, the third being for `1 +` written in code. Each at
      every door it
      can be written at: code as a formula term at `dbarts` and at `bart`, in a list on a formula, on
      a matrix and in `dbartsSpec`; a value in a list on a formula and on a matrix, in `dbartsSpec`
      and at `dbartsData(bases = )`. Each drawn and with `amplitude = fixed()`: one text, the basis's,
      asserted against the first door's but for `<arg>` and `<f>` (fail today: created). The message's
      example, taken from the text and run as written where it is the caller's own words, creates a
      sampler; a value's text shows no `forest()` call, at a forest door or the data door, and tells
      the caller of `cbind(1 - z, z)` to hand the factor over.
    - What stays. A factor, a character and a logical vector of two levels; `arm == "b"`; one numeric
      column as a vector, as a matrix of one column, and as `poly(dose, 1)`, `scale(age)`, `dose:age`,
      `I(dose * age)`, `0 + dose` and `dose - 1`: each created at every door, with the block step 1.1
      pins.
    - After the rows. `arm` under a `subset` that leaves two of its levels is created with those two,
      as code and as a value. Under a `subset` that keeps all three it is refused with the first [2]
      text at both. `cbind(1 - z, z)` through the data door under a `subset` that leaves one arm is
      refused, with the zero-column text, the rows being judged before the size: a build that drops a
      numeric column no row has, as it drops a level, makes a forest on one constant column there.
    - The forms the messages show. `y ~ forest(x1 + x2) + forest(x1, basis = arm == "b") +
      forest(x1, basis = arm == "c")` is three forests with those labels; `forest(x1, basis = dose) +
      forest(x1, basis = age)` beside a forest with no basis is three forests.
    - After creation. A numeric forest swapped to a matrix of two columns is refused, the sampler
      `identical()` to an untouched twin's afterwards and over the next 5 sweeps; a matrix of one column
      is taken. `predict(bases = )` with two columns for a numeric forest is refused in the same words.
      A held forest's width has no text of its own left.
    - Edited by hand, each at a named door. A data object whose `bases` entry is set to two numeric
      columns, and one whose recorded levels are set to three with a block to match, are refused at
      `new("dbartsSampler", control, model, data)`, at `dbarts(data)`, at `dbartsSpec(data)` and by
      `validObject(data)`, with the second and the first [2] text in the words for a value (fails
      today: the first is created by all three). An object with the slot removed creates at the three
      with a block of one column, and with a factor's block of two is two numeric columns and is
      refused. A `copy()` of a sampler whose level a swap emptied, and of one swapped to a column of
      zeros, is made and runs: the size is the record's and never the rows'.
    - The names, for later (dec-B272). The shapes written-surface pins against
      `names(coef(lm(y ~ 0 + <basis>)))` are read from [`basisColumns`](../../R/forestBasis.R), called
      on the model frame with no sampler, so the names of several columns stay pinned though no model
      has them; the refusal's text quotes those names.
2.7 Rewrite and remove, by step 1.0's trace; run the suite first and repair what it shows. The counts
    below are the traced commit's.
    - test-forest-basis-terms.R, 160 creations and its 12 predictions of several columns. Its fixture
      on g, three levels, moves to a column of two. Where a block tests the grammar or the names of
      several terms it reads the builder, as step 2.6's last test does, and fits nothing. Where it
      tests what a basis of several columns does in a fit or at new rows (the two-column rows of its
      blocks on columns, on rows and on new rows; the texts `dose + age` and `1 + dose`) it becomes
      that shape's refusal or goes, and the one-column rows stand. `scale()`, `poly(dose, 1)` and `I()`
      at new rows stand.
    - test-bcf-family.R, 30 creations through the data door: its arms on three groups and on a block
      times a constant go. They pin the tip's row norm on shapes that no longer exist; what they pin of
      a two-level block is held by the two-column arms the first commit made factors. Its 29 constant
      columns are one numeric column and stay until the multiplier slice.
    - test-forest-capture.R 12 and test-forest-labels.R 5 creations, with the 7 swaps of two numeric
      columns: a basis of `dose + age` used only as some basis becomes `I(dose + age)` or one of its
      columns, the label and the capture under test unchanged; the test of a swap's column order goes
      with the refusal.
    - test-forest-basis-r5.R 4, test-multiforest-leaf-prior-writer.R 4, test-formula-terms.R 3, and
      one each in test-argument-surface.R, test-bcf-creation.R, test-predict-blend.R and
      test-sampler-residuals.R: a factor of three levels becomes one of two, and a basis of several
      numeric columns one column, where the test is not about the size, and the test goes where it is.
    - test-forest-kind.R: the loop's column for push 2 is switched on.
    - test-forest-defaults.R: the held table's row for three or more columns takes the basis's text,
      and `cbind(1 - z, z)` through the data door, held second, turns from created to refused.
2.8 Help and records. man/forest.Rd: the `basis` item and its list of forms say what a basis is in
    this version, one numeric column or two levels, give the two refusals with the form each shows, and
    lose the lines that describe several columns, a constant beside a column, splines and polynomials of
    a higher degree as accepted; `sd` and `amplitude` lose what they say of several columns and of a
    width at a swap; its last example, `twoColumns`, which fits `basis = z + I(w / 50)` and would stop
    the page's examples at the refusal, becomes the form the refusal shows, a forest for each column,
    and reads each forest's column names as it reads the two today. man/bart.Rd, "Formula Terms": the
    line for `basis = a + b` becomes the two-forest form. man/dbartsData.Rd, `bases`;
    [`dbartsSampler$setForestBasis`](../../man/dbartsSampler-class.Rd) with its docstring;
    man/bartBT.Rd's `bases`: one column or two levels. docs/design/forest-kind-by-class.md: the two
    refusals, where the check sits, what is read and named and then refused, and what was kept for
    later with the addition each protects. docs/design/written-surface.md and public-surface.md where
    they say a basis of several columns is accepted. TODO: `factor-basis-refusal` and
    `several-column-basis-refusal` close; `factor-basis-per-level` and `several-column-basis-meaning`
    stay; `written-surface-leftovers` loses its sentence on the column names the data door keeps, one
    column or two levels leaving it no case.
2.9 Mutations (Verification).

## Verification

Every push, against the slice's own library (`R CMD INSTALL --preclean -l <lib> .`, `R_LIBS=<lib>` on
every call; check `dbarts:::buildInfo()$mode` and that the install postdates the source), run in series,
at most two cores (`MAKEFLAGS=-j2`, `EQUIVALENCE_CORES=2`):

- `cd tests/cpp && make && ./test_bartcore`: unchanged and passing (nothing under src/ moves).
- The full tinytest suite on the shipped build, in one process, counted per file: no failure, no file
  stopping. Both pushes: the landing note gives, file by file, the assertions removed, with a swap or
  a prediction the rule refuses in push 1 and with a refused shape in push 2, and those written in
  their place; a file whose count falls is named with what it no longer covers.
- On a reference build (`--preclean --configure-args=--enable-reference-build`): the four
  `test-reproducibility-*.R` files pass unchanged, and the three compares are bitwise, every scenario
  reporting identical draws, counted per scenario with no `max |z|` line: 55 against the gaussian
  baseline, 15 against `bcf-equivalence-1b7d730c.rds` (its respelled swap among them), 11 against
  `multinomial-equivalence-80b1c8d4.rds`. In push 1 the gaussian baseline is
  `equivalence-308c45be.rds` and nothing is re-recorded. In push 2 the two commits are compared with
  it, 54 identical and `bart2twoforest` not, and the baseline commit with the merged file, 55.
- Every command `.github/workflows/exact-gates.yaml` runs, 42 at the tip, in `quick`, unchanged: a data
  object and a fit carry another record, and push 1's class is posterior-changing.
- `inst/include/dbarts/dbarts.h` and src/ have no diff.
- The pair script (below), old side on the base build, new side on the slice's.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD build` with every vignette rebuilt and `R CMD check --as-cran` on
  a tarball from a clean copy. `Rscript benchmarks/R/composition-matrix.R` at push 2, no disagreement.
- `Rscript benchmarks/R/mutation-battery.R run m45,m46,m47` at both pushes, each of the three killed:
  the rows name test-forest-basis-r5.R and test-multiforest-leaf-prior-writer.R among their killers,
  and both files are rewritten here, so resolving anchors does not show that they still kill.
- inst/NEWS.Rd parses at push 2, by the reviewer's gate
  ([Reviewer checklist](README.md#reviewer-checklist)).
- The consumers, each suite whole against a private install, none failing; the landing note gives each
  count beside its count on the base build. Push 1: bartCause on dbarts-1.0 as it stands, and again
  with step B.1's edit. Push 2: bartCause with the edit; and once as it stood before it, which must fail
  at `bcf()` with the second [2] text and nowhere else, the proof that the edit is what the push needs.
  stan4bart on bartcore, treatSens on dbarts-1.0 and bairrtt on main at both pushes, unedited.

The pair script. Each row is one model fitted on the base build as written there and on the slice's as
written here, the same seed: sampler fits compare the train draws, sigma and the coefficients, `bart`
fits `yhat.train`, sigma and `predict` at 40 new rows. `identical()` on every row but those marked.
Rows 01 to 32 and 39 to 44 are push 1's against the tip before it; rows 33 to 38 push 2's against
push 1.

| | base build | slice |
|---|---|---|
| 01 to 06 | a formula term on a factor, an ordered factor, a character, a logical, a 0/1 number, dose | the same text |
| 07 to 10 | a list on a formula, as values: a factor, a logical, `cbind(1 - z, z)`, two numeric columns | the same text |
| 11 to 13 | the matrix interface; `dbartsSpec()`; indicator columns through the data door | the same text |
| 14 | `dbartsData(bases = list(NULL, cbind(1 - z, z)))` | `bases = list(NULL, factor(z))` |
| 15 | a logical vector in a `forests` list | the same vector through the data door: the first changed sequence against its oracle |
| 16 | a logical vector through the data door | the same text: NOT identical |
| 17 | a factor of three levels written as code, under a `subset` that leaves two | the same factor as a value: the fourth changed sequence against its oracle |
| 18 | that factor as a value under that `subset` | the same text: NOT identical |
| 19 to 22 | probit; logistic; three forests of the three kinds; `subset`, weights of which ten are zero, and an offset | the same text |
| 23, 24 | a swap to `cbind(1 - z2, z2)` on a factor forest, 25 sweeps either side; a swap to `factor(z2)` | `factor(z2)`; the same text |
| 25 | a swap to `factor(v, levels = <those at creation>)` | the same vector with its levels in another order: the second sequence against its oracle |
| 26 | a swap to that vector in the other order | the same text: NOT identical |
| 27, 28 | `predict(bases = )` with indicator columns; with the aligned factor | with the factor; with its levels reversed |
| 29, 30 | a held factor through a list and its block through the data door; `bart` with a term and predict from the stored record | the factor at both; the same text |
| 31 | a copy after a swap, 10 sweeps | the same text |
| 32 | bartCause's `bcf()` as it stands, the settings of step B.2 | with step B.1's line |
| 33 to 35 | the three kinds that stay, at a formula term, a list and the data door | the same text |
| 36 | a model of three forests, `arm == "b"` and `arm == "c"` | the same text |
| 37 | `cbind(1 - z, z)` through the data door, drawn and held second | `factor(z)` there |
| 38 | a factor of three levels; `dose + age` | the same text: refused |
| 39, 40 | a number, and a factor, with a missing value on a row `subset` drops, written as code | the same handed over as a value in a list |
| 41 | a number with a missing value on a row a missing response drops, written as code | the same through the data door on the matrix interface |
| 42 | `predict(fit, x, bases = fit$bases)` at the training rows of a factor fit, the block handed back | the same text, the factor handed back |
| 43 | `extract(fit, type = "forest", contribution = TRUE)` and `fitted` on a factor fit | the same text |
| 44 | a sampler of one forest, 10 sweeps | the same after a refused `$setData` of a data object with bases |

Mutations, each expected to fail the named test and no gate before it:

- push 1, a block of 0s and 1s with one 1 a row read as levels at the data door: step 1.1's pins beside
  the loop;
- push 1, a logical vector coerced to a number at the data door; at a forest door: step 1.1's loop;
- push 1, a matrix of logicals of one column read as its vector at the data door; a character matrix of
  one column read as levels where it is written as code: step 1.1's loop at that door;
- push 1, the data door keeping a rule of its own for a factor (another order of levels): step 1.1's
  `identical()` blocks;
- push 1, a missing value judged on every row handed over, two mutations (a value in a list; the data
  door on the matrix interface): step 1.1's rows of a missing value at that door;
- push 1, a missing value on a kept row let through where the rows are final: step 1.1's refusals on
  a kept row;
- push 1, NA taken as a level at a value door; left with the tip's text at a code door; refused where
  no kept row is at it: step 1.3's two tests of the NA level, and step 1.1's `addNA()`;
- push 1, the levels dropped where rows are restricted on the matrix interface; at the formula door:
  step 1.2's `subset` rows;
- push 1, `dbartsSpec()` and a data object's levels, three mutations (a declaration leaving the levels
  the object had; a list that states no basis clearing them; a numeric declaration leaving them): step
  1.2's `dbartsSpec()` tests, each direction;
- push 1, the levels left on the stored block as an attribute: step 1.2's `attributes()`;
- push 1, a code basis's levels not recorded (its forest read as numeric): step 1.2's slot at the
  formula door, and step 1.4's table for a forest written as code;
- push 1, the count of the levels unchecked at one door, three mutations (the sampler's `initialize`;
  `dbarts()` handed a data object; `dbartsSpec()` handed one), and the validity's call left out: step
  1.2's edited object at that door;
- push 1, `validObject()` called on the data object where a sampler is made: step 1.2's object with the
  slot removed;
- push 1, a fit storing the block for a forest of levels; storing a factor with the levels of the value
  handed over and not those recorded; packaged from a data-door object, `bart(data)`, without its
  levels: step 1.2's tests of the fit, the one under a `subset` among them, and step 1.4's prediction
  by name on that fit;
- push 1, a value's emptied level kept as a column of zeros at one door, the matrix interface: step
  1.3's `identical()` to the code spelling;
- push 1, the levels dropped before `subset`, every level of the data kept: the same test;
- push 1, a level dropped for having no row of positive weight: step 1.3's weight-zero row;
- push 1, the zero-column check left out where `dbartsSpec()` installs; at the data door: step 1.3;
- push 1, the rows' function called from the sampler's `initialize`; from a reload's re-creation: step
  1.5's `copy()`, and its reload, after a swap that empties a level and after one to a column of zeros;
- push 1, a swap matching levels by position: step 1.4's reversed levels and lost level;
- push 1, a level matched without regard to case; with white space trimmed; by its first letters: step
  1.4's three near names;
- push 1, a level the forest lacks appended as a new column: step 1.4's refusal;
- push 1, a level declared and on no row refused at a swap: step 1.4's swap of the factor of creation;
- push 1, the recorded levels overwritten by those an accepted swap declares: step 1.4's slot after
  the swap that declares one more;
- push 1, a one-sided formula at a swap conformed by no table; a bare value at `predict` likewise: step
  1.4's table, its formula cells and its bare cells;
- push 1, `predict` left on the base build's expansion: step 1.4's `predict` rows;
- push 1, `predict` refusing a column of zeros for a numeric forest: step 1.4's cell for it;
- push 1, a factor's levels put on the wrong columns the same way where a fit is packaged and where
  `predict` reads it (the stored factor built with its levels in sorted order): step 1.4's training
  rows, whose factor is not in sorted order;
- push 1, indicator columns, the sampler's own block among them, accepted on a factor forest: step
  1.4's table and its test of the block handed back;
- push 1, the data object written before the table is consulted: step 1.4's twin after a refusal;
- push 1, the forest with no basis taken to be forest 1: step 1.4's last test;
- push 1, `copy()` handing on a data object without the levels: step 1.5;
- push 1, a state's coefficients matched to the recipient's level names at an install; the donor's
  levels handed to the recipient: step 1.5's state test;
- push 1, the check of `$setData` left out: step 1.7;
- push 2, a fixture's factor given other levels than the block's column order (`factor(z, levels = 1:0)`):
  the edited file's own pinned draws and pins of `data@bases`;
- push 2, the size check left out at one door, four mutations (the data door; where `dbartsSpec()`
  installs; a value in a list; a basis written as code): step 2.6's loop at that door;
- push 2, the size judged before the rows: step 2.6's `arm` under a `subset` that leaves two levels, as
  code and as a value;
- push 2, a numeric column no kept row has dropped as a level is, where it is refused: step 2.6's
  `cbind(1 - z, z)` through the data door under a `subset` that leaves one arm. Push 1 has no test
  that tells this build from a right one;
- push 2, two columns let through when they are 0s and 1s with one 1 a row: step 2.6's loop for
  `cbind(1 - z, z)`, and bartCause as it stood;
- push 2, the text for `1 +` chosen by a first column of ones: step 2.6's `cbind(1, dose)` as a value;
- push 2, a value's text showing a `forest()` call, or not naming the factor for indicator columns:
  step 2.6's loop, the texts at the value doors;
- push 2, a held basis of three levels let through to the hold's own text, or to the engine: step 2.6's
  loop with `amplitude = fixed()`;
- push 2, a one-column matrix refused as several columns: step 2.6's "what stays";
- push 2, the swap's width left free; the prediction's: step 2.6's "after creation";
- push 2, the size unchecked in the plain function at one door, three mutations (the sampler's
  `initialize`; `dbarts()` handed a data object; `dbartsSpec()` handed one), and the validity's call
  left out: step 2.6's edited objects at that door;
- push 2, the plain function judging the rows, a level or a column with none: step 2.6's `copy()`
  after each swap;
- push 2, the message's example written with another forest's predictors or with a level that is the
  first: step 2.6's run of the example.

Not a hot-path change: nothing a sweep runs is touched.

## Budget and stops

Planned ~4350 lines, by step:

| step | R | tests | the rest |
|---|---|---|---|
| 1.0 the trace, again | | | its counts, in the landing note |
| 1.1 one reader of a class | 100 | tinytest 250 | |
| 1.2 the record, its check and the fit | 125 | tinytest 180 | |
| 1.3 the rows at creation | 55 | tinytest 85 | |
| 1.4 after creation | 120 | tinytest 240 | |
| 1.5 the same sampler again | | tinytest 85 | |
| 1.6 respell | | tinytest 270 | benchmarks 5 |
| 1.7 `$setData` on a sampler of one forest | 10 | tinytest 25 | |
| 1.8 help and records | | | help 185; design note, public-surface, architecture, TODO, indexes 230 |
| 2.1 the blocks written by hand | | tinytest 380 | |
| 2.2 benchmarks and help | | | benchmarks 10; help 25 |
| 2.3 the baseline scenario | | | benchmarks 30; MANIFEST row and repins 15 |
| 2.4 the size at creation | 130 | | |
| 2.5 after creation | 50 | | |
| 2.6 the new file | | tinytest 425 | |
| 2.7 rewrite and remove | | tinytest 950 | |
| 2.8 help and records | | | help 210; design notes, TODO, indexes 160 |

R ~590, tests ~2890, the rest ~870; push 1 ~1965 and push 2 ~2385. Plans here run 1.5 to 2 times low:
the sd unit slice was planned at ~2340 and stood at ~3320 before its review, its tests at 1.7 times
their figure and its R code at its figure, and the four slices before it ran 1.8 to 2.3 times theirs.
Assumed here: the tests run to about twice their figure and everything else to about 1.3 times, which is
~7700 lines, inside a forecast of 6600 to 8700. An overrun in tests is expected and is not a reason to
stop; one in R code means the design moved. bartCause's part, ~15 lines, is outside these counts.

Stop and report when:
- push 1's diff passes ~3900 lines or its R code ~615; push 2's passes ~4750 or its R code ~270; the
  two together ~8700 or ~885 of R code: twice a diff's figure, 1.5 times the R code's;
- step 1.0's trace gives a count that sizes a step, the creations of a refused size, the blocks of two
  indicator columns, the swaps or the predictions the rule refuses, more than a fifth over Context's:
  before anything is built;
- anything under src/ or `inst/include/dbarts/dbarts.h` needs a change, or the bridge must be told a
  kind;
- a recorded scenario other than `bart2twoforest` moves in any of the three compares, or a snapshot
  file moves;
- a pair row that must be `identical()` is not, or one marked NOT identical is;
- a class must be tested anywhere but the one function, or a door needs a reading of its own;
- the levels cannot be kept off the stored block, or `data@bases[[f]]` must change what it holds;
- a check of a level's or a column's rows must stand where `copy()`, a reload or a sampler made again
  from its own data runs it;
- anything of a fit but the `bases` entry of a forest of levels must change, or one of its readers
  cannot go through the one function;
- a test fails in a file step 1.0's trace does not name for that push;
- more swaps or predictions meet push 1's rule than step 1.0's trace counts, by half again, or a
  creation of the suite meets one of its four changed sequences;
- bartCause unedited fails on push 1's tip, or fails on push 2's anywhere but `bcf()`'s creation.

Expected verdicts: tests/cpp unchanged; tinytest clean; the three compares 55, 15 and 11 identical on
push 1, and on push 2 54 with `bart2twoforest` moved, then 55 against the merged file; the four snapshot
files unchanged; all 42 exact-gate commands passing; the pair script's rows identical but 16, 18, 26 and
38; mutation-battery rows m45 to m47 killed at both pushes; the consumers' suites clean, bartCause as it
stood failing push 2 at `bcf()` alone.

## NEWS

No new item: forests, `dbartsData(bases = )` and `$setForestBasis` are new in 1.0-0 and nothing released
changes. One line of the forests item is corrected, in step 2.2.

## What a consumer sees

Read in each package's source at the tips Context names; no consumer suite was run.

- bartCause. Push 1 changes nothing for it: `bcf()`'s block through the data door is two numeric
  columns before and after, the fit a two-level factor has, and its own check refuses a treatment the
  kept rows leave with one value before dbarts sees a column of zeros. With step B.1's line it hands a
  factor, the same fit to the bit on push 1's tip. Without the line push 2 refuses `bcf()` at creation
  with the text for two numeric columns. Its coefficients keep their order, the levels 0 and 1, and the
  names it gives them itself. It refuses a treatment with a missing value itself, so the rows a
  missing value is judged on change nothing for it; it fits through a sampler and reads no fit's
  `bases`. The sd unit's edit to the same files goes first.
- stan4bart, treatSens and bairrtt see nothing: one forest, no basis, no swap, no `bases`; stan4bart
  refuses `forests` by name in its `bart_args` and hands `dbarts()` a list of one forest in its
  multivariate fitter. Their data objects gain a slot that is `NULL`, and one saved without it reads as
  before.
- A caller of dbarts. Push 1: `dbartsData(bases = )` takes a factor, a character vector and a logical
  vector as `forest(basis = )` does, and no longer takes a matrix of logicals, a column of zeros or a
  logical with one value; a factor handed over with a level no kept row has is created without it; a
  basis handed over with a missing value on a row the fit drops is taken, as the same basis written as
  code is; a factor with NA as a level is refused, with the remedy; `$setForestBasis` and
  `predict(bases = )` take the class a forest was created with and match a factor by its level names,
  so a factor forest no longer takes indicator columns, its own stored block among them; a fit's
  `bases` holds a factor for a forest of levels, which `predict(bases = )` takes back, and code that
  multiplied that entry into the coefficients takes `model.matrix(~ 0 + b)` of it or asks
  `extract(type = "forest", contribution = TRUE)`; `$setData` on a sampler of one forest refuses data
  that carry bases. Push 2: a basis of three or more levels or of several numeric columns is refused at
  every door, with the form to write.

## What this leaves for the multiplier slice

- The prior itself. After this slice a factor and a number differ in what they are called, in what a
  swap and a prediction take, and in nothing the sampler draws. The multiplier slice gives a number its
  own coefficient variance and default; it reads the kind from the one function this slice adds, and
  the bridge will need it at every construction, which this slice does not hand over.
- The creations on a 0/1 number and the 29 on a constant column of test-bcf-family.R: left as numbers,
  and the multiplier slice moves their priors. Which of them pin a literal draw is its plan's to count;
  a constant column is refused by it where no sd is stated.
- The held shapes: forest-defaults-by-kind's two interim refusals by position and dec-B353's of one
  numeric column stand. Which shapes may be held, and at what value, is the multiplier slice's
  (dec-B354).
- The sd unit's record of a forest created without a basis (`hasBasis`, and
  [`ForestSpec::createdWithoutBasis`](../../src/bartcore/combiner.hpp) in the engine) says nothing the
  data object does not once push 1 lands, no forest gaining or losing a basis afterwards. It stays,
  this slice changing nothing the bridge reads; the multiplier slice drops it or keeps it.
- A sampler saved before this slice has no record and reads as numbers wherever it has a basis. No
  release made one, and no refusal is written for one. From push 2 one whose block has two columns
  reloads and runs, the reload not being checked, and a `copy()` of it is refused for its size.
- The reader's entry for the kind, the printed line ("factor basis, 2 levels") and the help's advice
  that a treatment belongs in a factor (dec-B263): with the multiplier slice, where they become true of
  the prior.
- What the multiplier slice may assume: every forest has one size, one numeric column or two levels,
  and keeps it for life; every block in the tests, the benchmarks and bartCause that stands for a
  factor is a factor; no accepted swap or prediction changes a forest's kind or a factor's columns;
  a fit's `bases` entry for a forest of levels is the factor, with the levels of the record; nothing a
  1.0 user can write has more than one size for a forest.

## What is kept for later, and what is not built

- Kept: written-surface's reading and naming of a basis of several terms, to the point of the refusal.
  It protects the later meaning of several numeric columns (TODO `several-column-basis-meaning`) and a
  size for each column (dec-B272): the names a 1.0 message shows are the names those forms will use.
- Kept: the engine's coefficient for each column of any block, which is what it has; a two-level factor
  is two columns there. Not narrowed, so that a forest for each level (TODO `factor-basis-per-level`),
  which is R's to build from two-level forests, and either meaning of several columns need nothing
  taken back.
- Kept: `basis.levels` as a list of level vectors of any length, and a data object's `bases` as
  matrices of any width; only creation refuses a size.
- Not built: any meaning for three or more levels or for several columns; a width rule for a held
  forest beyond "none changes"; a refusal by content. Each later meaning is an addition: every call it
  would serve is refused by name in 1.0-0.

## What waits on what

- On the sd unit: this plan is written to its final tip, 37b59ea0. It shares R/model.R, R/spec.R and
  R/dbarts.R with this slice, and one test of its file is respelled here (step 1.6). Step 1.0 reads
  the tree push 1 is built on; nothing else of it is waited on.
- On written-surface and forest-defaults-by-kind, landed: a code basis is built through R's model frame
  on the rows kept and drops its empty levels, so step 1.1's one function is called from
  [`basisColumns`](../../R/forestBasis.R) too, step 1.2 records a code basis's levels from what
  [`buildCodeBasis`](../../R/forestBasis.R) returns, and step 2.4 is called where its block is put in
  place. A forest with no basis may stand anywhere in a list and in a data object's `bases`, so "none"
  is never "forest 1"; `<f>` is named by label; three of forest-defaults-by-kind's interim texts are
  replaced or retired by push 2 and two wait for the multiplier slice. dec-B284 has every basis read
  when `forest()` is called; it changes neither the class nor the size of what is read, and
  `$setForestBasis` is outside it.
- On bartCause: step B.1, between the pushes and after the sd unit's edit there. Push 2 is not pushed
  before it.
- Order: the sd unit, this slice, the multiplier slice. This slice and the multiplier slice are not one:
  that one is posterior-changing for every numeric multiplier, and the respelling of push 2 can be
  proved bit for bit only on a tip where the kind is recorded and the prior has not moved. This slice
  and the sd unit are not one either: that one changes engine arithmetic.

## Out of scope, and where it goes

- The multiplier slice, with everything listed above; `updateBasisScale`, after the merge to main
  (dec-B276).
- A forest for each level of a factor of three or more (TODO `factor-basis-per-level`), and what several
  numeric columns mean (TODO `several-column-basis-meaning`): after the release.
- A data frame as a basis, and a factor beside numbers in one basis: refused, as today; additions.
- `$setModel` and `$setData` on a sampler of several forests: refused whole today, so no model or data
  of another kind can reach one. The one way bases reached a sampler made without them, `$setData` on
  a sampler of one forest, is closed in step 1.7. The control-migration arc decides what the two take.
- A check of a sampler rebuilt at a reload: none is written.
- To TODO as new entries: a level that is to appear only later in a run cannot be declared, a level no
  kept row has being dropped at creation; `predict(bases = )` checks no name on a numeric column; the
  sd unit's record of a forest created without a basis has no case left from R; a sampler has no
  method that returns a forest's basis or its levels, which are read off its data object; and, found
  by the probes and no part of this slice, `bartBT()` handed a data object and no response stops on
  its own signature ("argument "y.train" is missing"), whatever the object carries.

## Open calls

None. A swap never changes what kind a forest is (dec-B452), and at predict, as at a swap, numbers are
refused for a forest created on a factor (dec-B453). What the blind critique raised beyond the rulings
is settled, and is in "Calls made in planning" with its alternative and its cost.

## Calls made in planning

- A matrix of characters, and one of logicals, of one column is refused at every door: a vector is
  levels, a matrix is numbers or nothing. Run on the build of 37b59ea0: the character matrix is levels
  at the forest doors, as code and as a value, and refused at the data door; the logical matrix is
  refused at the forest doors and taken as a 0/1 number at the data door. The alternative reads each
  as its vector, as `lm` reads both. The cost: 6 cells of the first table turn from created to refused
  for the character matrix and the data door's one for the logical; nothing hands either over (no
  test, benchmark, help example or consumer).
- A missing value in a basis is judged on the rows the fit keeps, at every door. The alternative
  judges it on every row handed over, as a value at a forest door is judged today: the data door would
  then refuse fits it creates today, typically a treatment missing where the response is, and one
  column would be two answers by whether it was written as code. The cost: 13 cells are created where
  they are refused today, and a caller is not told that a basis held missing values on rows the fit
  dropped. `lm` judges after `subset`.
- NA as a level is refused at every door, code included, where a kept row is at it. The alternative
  takes it as a level, as `lm` does, and needs a rule for matching NA by name at a swap and at
  predict. After push 2 a factor has two levels, so the case is one named level against NA; and a
  refusal can be lifted after the release, where an accepted meaning cannot be withdrawn. The cost: a
  value fit created today is refused, 3 cells, with the remedy in the text. An NA level with no kept
  row at it is not that case and goes with the levels no row has (dec-B308), as code has it today:
  refusing it too would refuse `addNA()` of a complete factor, which is created as code today.
- A fit stores the factor for a forest of levels, not its indicator block. Under dec-B452 and
  dec-B453 a factor forest refuses the block, so a fit that stored it could not hand its own `bases`
  back to `predict`, which the tests' first oracle of a prediction does. The alternatives: keep the
  block and say in the help that it cannot be handed back, or store the factor beside the block. The
  cost: three readers in R/generics.R, a shared helper of the tests and three test files that multiply
  the entry into the coefficients go through one function, and the pins of a factor entry's shape,
  values and names are rewritten in three more; no benchmark, help example or consumer reads the
  entry. The data object keeps the block, the engine's form, so the caller of a swap keeps the factor
  it handed over, which the help says.
- A state installed from a sampler whose factor has the same levels in another order goes by position.
  A state carries no data (dec-B231, dec-B418), and the levels are data. The alternative refuses an
  install where the two samplers' recorded levels differ, which needs the levels in the state. The
  cost: for a caller who moves a state between two such samplers the two levels' coefficients exchange
  with no message; the help says so in one sentence.
- A factor matrix is refused at every door, at a swap and at predict, of one column or
  of several, as a matrix of characters or of logicals is: a vector is levels, a matrix is numbers or
  nothing. Run on the base build: of one column it is levels at the forest doors, as code and as a
  value, at a swap and at predict, and refused at the data door as not numeric. The alternative reads
  it as its vector. The cost: the forest doors, a swap and a prediction turn from taken to refused for
  it; nothing hands one over (no test, benchmark, help example or consumer). An array of one dimension
  (`dim(f) <- n`), a factor or a character or logical one, is its vector at every door, as in a model
  frame, as on the base build a factor or character one was levels; there a logical one was refused as
  a logical matrix. The alternative refuses all three with the matrices; the cost of the call is none
  found.
- `dbartsSpec()` replaces a data object's blocks and levels whole when a forest of its list declares a
  basis, as the base build replaces the blocks (run on both). The alternative replaces forest by
  forest and keeps what the object carried for a forest the list gives no basis. The cost of the
  alternative is code for cases no test or consumer has, the grammar refusing most partial lists; the
  cost of this is that a data object of three forests under a list of two becomes two forests, and a
  list that moves the forest with no basis drops the block the object had at that place, with no
  message.
- `$setData` on a sampler of one forest refuses a data object that carries bases, in push 1. It is a
  defect of the tip and not of the rule: the run is the plain model's and `copy()` then fails. The
  alternative leaves it to the control-migration arc, with a TODO entry. The cost: one check, one text
  and one test.
- The record is checked by a plain function with four callers, and the rows are no part of it. No door
  that makes a sampler runs the class's validity, and `validObject()` there would refuse an object
  saved without a slot; a check in `initialize` runs on every `copy()`, so only what holds for a
  sampler's life stands there. A reload's re-creation is left unchecked: a check there would be a
  third copy of the call for an object edited in a file.
- A 0/1 number is not touched. It is one numeric column at all seven doors, and has been since the
  basis was built, as dec-B260's own text says ("a 0/1 number another"). No tip between this slice and
  the multiplier slice fits a model nobody asked for, because this slice changes the fit of no number.
- A factor is matched at a swap and at predict by the levels its rows carry. dec-B289 rules the match
  by name and dec-B308 refuses a "basis swap that carries the dropped level"; read as the levels the
  value declares, the factor a forest was created from would be refused at its first swap. `predict()` for
  `lm` takes a declared level no row has and refuses a row at a new one (run).
- A row of weight zero, or masked out, counts as a row for a level. dec-B308 says the rows the fit
  keeps; the rows in the likelihood are another set, the one a scale is read from, and a level whose
  rows are out of the likelihood for a while is the case dec-B346 takes at a swap.
- A value reaches the data object as it was given, and is expanded there. The forest doors expand a
  value before the data object cuts its rows today, which is why a value's emptied level stays a column
  of zeros where code drops it, and why `dbartsSpec()`, which cuts nothing, refuses what the list door
  keeps: one object read two ways by where it came in.
- A refusal the tip already words keeps its words, and the doors that lack it join it: a column of all
  zeros, a factor of one level, a logical with one value, a missing value. Only what is refused for
  the first time gets a new text, and a new text names its forest and advises what runs at its door.
- The kind is stored once, on the data, as the levels. The data door makes a data object with no
  control, so the record cannot live on the control alone; with levels matched by name neither a
  forest's kind nor its levels can change after creation, so the control's copy for a code basis
  (`basisTerms`) and its `hasBasis` cannot come to disagree with it.
- The levels are a slot and not an attribute of the block: a cut by rows drops an attribute, and every
  pin of `data@bases` against a bare matrix would see one.
- The rows are judged where the data object has them, and again where `dbartsSpec()` installs a
  declaration, the two places a block becomes final. A sampler made again from a data object whose swap
  emptied a column is not a creation and is not checked, nor is a copy or a reload: dec-B346 puts such
  a column in a sampler's data for good.
- The data door's column names are left as they are. It keeps a partly named block's names where a
  value at a forest door is named whole or not at all; with one numeric column or two levels the two
  cannot differ, so push 2 leaves the difference no case and nothing is decided for it.
- Both refusals of dec-B281 and dec-B282 land here, together, in the second push, at every door at
  once. The forest doors know a basis's class today and could have refused sooner; the data door cannot
  tell a factor's indicator columns from two numbers until push 1, and bartCause hands exactly that
  block through it. A refusal at some doors before the others would be one object read two ways by
  where it came in, the fault this slice exists to remove.
- bartCause moves between the pushes and not "the same day": the day a refusal lands is too late for a
  consumer whose only spelling it refuses, and push 1's tip takes both spellings.
- Push 2 is two commits and one gate battery, and its baseline a third commit after review. The
  respelling must be proved on a tree that still accepts both spellings, so it cannot be one commit
  with the refusal; a third push would buy a second battery for a commit that is tests and one help
  line.
- The size is refused in R and not in the engine. The engine keeps a coefficient for each column of any
  block; narrowing it would have to be undone for either later meaning.
- A basis of several terms is still read and named before it is refused, and the names are tested off
  the builder. dec-B272 fixes the names for whatever several columns come to mean; a refusal made
  before the columns exist could not show them and would leave the naming rule untested until then.
- The messages write the example in the caller's own words where they are code, and say it in words
  otherwise: a value door has no code of the caller's, a `forest()` call is not the data door's
  spelling, and `<arg>` is a list at two doors. For one term of several columns, a spline, there is no
  form to show that rebuilds its knots at new rows from one column of it, so the text does not suggest
  a subscript.
- `1 + dose` has its own sentence: its two columns are not two things the caller listed, and "a forest
  for each column" would tell them to multiply a forest by a constant.
- A block of three indicator columns in a test is not respelled as a factor of three levels in the
  first commit to be refused in the second; its test goes in the second.
- The baseline's scenario is given a factor of two levels, not the two-forest form the message shows:
  it is there as the one model of several forests the formula door reaches in that file, K = 2.
- The composition matrix gains no cell: it states which families and features compose, and a basis's
  kind is no row or column of it. Its runner is edited in one place, the `bart()` door of its `bcf`
  row, so that the row is one model at every door.
- The pair script is run at landing and not tracked: its old side needs the base build. What stays in
  the suite is step 1.1's loop across doors, step 1.4's identities and step 2.6's loop.

## Recheck against bartcore

On 2026-10-11, against e88b78eb, the tip of the sd unit slice: the plans and notes of what landed since
this plan was written were read ([written-surface.md](written-surface.md),
[forest-defaults-by-kind.md](forest-defaults-by-kind.md),
[small-rulings-1008-r.md](small-rulings-1008-r.md), [forest-sd-unit.md](forest-sd-unit.md) with its
design note, [forest-multiplier.md](forest-multiplier.md),
[feature-matrix.md](../design/feature-matrix.md) and its runner), then the code, the ledger for every
ruling cited here and every entry on a forest, a basis, a level, a swap or a prediction's bases through
dec-B447, and the runs below on a shipped build of that tip. The text above is written to that code
throughout.

What changed:

- What the plan waited on is in the tip: written-surface's third push, forest-defaults-by-kind and the
  sd unit. Every call is in the landed spelling, and a value is what `do.call()` hands over: a variable
  named in `forest()` is code.
- A basis is read in R/forestBasis.R now. Step 1.1's "one place" was
  [`expandForestBasis`](../../R/model.R) alone and is one function called from five; the plan's
  "block that arrives expanded keeps its levels" is a value handed on as given; steps 1.2, 2.4 and 2.6
  name the builders by their names.
- dec-B289 ruled the plan's one assumption, and it is worked into the rule and step 1.4. dec-B308
  replaced the refusal of a value's empty level with a drop in every spelling: the droplevels text
  goes with its rows of "The rule" and "Refused forms", step 1.3 is rewritten, and a fourth changed
  sequence is stated, since the tip keeps a value's emptied level as a column of zeros (dec-A175 and
  dec-B380 describe it as dropped; run, it is not). dec-B346, ruled on 2026-10-08 and built, takes a
  column of zeros at a swap: the plan's call is the tip's behavior and its test is no longer one that
  fails today.
- The sd unit records `hasBasis` and has the engine remember a forest created without a basis, because
  a swap can give such a forest one. Push 1 leaves that no case from R; the record stays, one of its
  tests is respelled, and the multiplier slice is told.
- Counts. 1412 models of several forests in 67 test files where the plan had 776 in 53 at its tip and
  1186 on written-surface's build; 222 forests of a refused size in 11 files where it had
  243 in 12; 91 two-column indicator blocks in 14 files where it had 98 in 10; 57 swaps reaching a
  sampler in 12 files, 27 of them against push 1's rule, where it had 29 in 7 and 20; 85 predictions
  with `bases` in 6 files, 29 against the rule, where it had 64 and 14.
- Names: the gaussian baseline is `equivalence-308c45be.rds`; the exact-gates workflow runs 42 commands;
  the refusal of a held numeric column is dec-B353's; the consumers' tips; the baseline is re-recorded
  in a commit of its own with its repins, as the last one was.
- The composition matrix's runner already multiplies its `bcf` samplers by a factor; one door of the
  row is left on a 0/1 integer, step 2.2's.
- The texts for a column of zeros and for an empty level, which the plan reworded, stay the tip's.
- The budget: ~3840 planned for ~3200, with a forecast and stops where there was an upper figure of
  5600.

The rulings as they stand. dec-B260 is the rule, and nothing since revises it. dec-B254 stands and
names `$setModel`. dec-B281 and dec-B282 stand and are push 2; dec-B280 as dec-B282 leaves it is built
in the sd unit; dec-B272 stands and is why a basis of several terms is still named. dec-B289 and
dec-B308 are built here. dec-B452 and dec-B453 are built here. dec-B346 is built, and revises dec-C14. dec-B353 and dec-B354 concern a held
coefficient and stay the multiplier slice's; dec-B354's "more than two levels or several columns held
stay refused" agrees with push 2. dec-B284, dec-B351, dec-A175 and dec-B380 describe the code path this
plan stands on. dec-B316 has a label matched exactly, which is how `<f>` is named. dec-A185 refuses
`bases` at `predict` with no new rows, so the table at predict is for new rows only. dec-B290, dec-B291,
dec-B293, dec-B303, dec-B307 and dec-B362 touch a forest's default, its scale or a re-derived response
scale, and nothing here. No entry after 2026-10-08 touches a forest's kind, a factor's levels, the data
door, a swap or a prediction's bases.

Runs, each one short process on the shipped build of e88b78eb: the 36 objects at the seven doors; the
rows probe (a level, a column and a logical emptied by `subset` and by a missing response, the code
shapes of several columns, the two forms the size refusals show); 8 kinds of forest by 17 kinds of value
at `$setForestBasis`, a forest with none given a basis at either position, and the block after a swap of
reversed levels; `predict(bases = )` on a factor fit and on a numeric fit, and the records on a data
object, a control and a fit; a level at weight zero, a copy, a reload, `new()` from a sampler's three
and from an edited data object, `$setData` and `$setModel`; `predict()` for `lm` on a declared level
and on reordered levels; and the 72 test files that name a basis, each traced in its own process with
the functions that read a basis replaced in the loaded namespace. Nothing was compiled or installed.

Not verified: nothing of the slice is built, so no mutant and no new test was run, and the 252 cells'
outcome under the rule is read off the rule and the tip, not measured. The prototype of the rule that
this plan was first measured on (227 of 252 cells with the same draws, the suite's files stopping at
their first refusal, bartCause's three-way identity over 11 seeded fits) predates written-surface's
third push and was not rebuilt. The tinytest suite was not run whole: its count at this tip is the sd
unit's gate run's, and the trace covers the files that name a basis, where one assertion of
test-predict-blend.R fails under the trace and passes without it. Whether the fourth sequence moves the
posterior of the fit or only its draws was not established; it is classed posterior-changing. No
reference build, bitwise compare, exact gate, sanitizer, x86 run or `R CMD check` was made, and the
composition matrix's runner was read and not run. The consumers were read at the tips named and not
built or run.

## Revised from the blind critique

On 2026-10-11 a blind critique read this plan and ran its own probes on the shipped build of 37b59ea0.
Each finding was checked against the source or by running its probe again before it was acted on. All
held, with the details at the end. What changed, by finding:

- B1, the hand-edit refusals rested on a validity no creation door runs. The count of the levels and
  the size are one plain function, called from the sampler's `initialize`, from `dbarts()` and
  `dbartsSpec()` handed a data object, and from the validity; the rows' check stands at a creation
  only, never where `copy()` or a reload runs it; a reload is said to be unchecked; each test of an
  edited object names its door, and the mutant on the validity is one for each caller ("The rule",
  steps 1.2, 1.3, 1.5, 2.4, 2.6).
- B2, the rule did not say on which rows a missing value is judged. It is the rows the fit keeps, at
  every door; the 28 cells are in Context, 13 of them created where refused; step 1.1's loop has the
  objects under `subset` and under a missing response.
- S1, the counts were of the sd unit slice before its last rebase and review round. Step 1.0 runs the trace again
  first; the stops are relative to its counts; Context says which commit each figure is of, and the
  figures a search gives are those of 37b59ea0.
- S2, step 1.1 against "Refused forms". The one function answers the kind or names the class, and each
  caller words its refusal; the code door's text for a matrix of one column is given.
- S3, the last example of man/forest.Rd fits a model push 2 refuses: named in step 2.8.
- S4, `$setData` on a sampler of one forest took bases in silence: refused in push 1, step 1.7, with
  its test and mutant.
- S5, a column of zeros at predict on a numeric forest: taken, as today; in the table and step 1.4.
- S6, tests that need numbers moved to a forest push 2 refuses: they move to a numeric forest of one
  column in step 1.6, and step 2.1 respells a numeric swap or prediction with the forest it respells.
- S7, the texts: each names the forest, advises what runs at its door and builds no code from
  `<arg>`; none advises `as.numeric()` of a factor; a forest's levels are named where a factor is
  asked for; the caller of `cbind(1 - z, z)` is told to hand the factor over; a level is in single
  quotes and the unknown level's text has a remedy; a forest with none has one text at both doors.
- S8, a factor with NA as a level: refused at every door where a kept row is at it, with one text;
  dropped where none is.
- S9, verification: mutation-battery rows m45 to m47 are run at both pushes, and push 1 accounts for
  its assertions file by file as push 2 does.
- The minor items. The first table's row for a data frame is spelled out and the NA level added, 37
  objects; the benchmarks' bases on a number are four, not two; a value of several numeric columns
  with one emptied is under "Refused where accepted"; the NEWS line is corrected in step 2.2; push 2's
  three commits go up in one push; the consumers' tips and stan4bart's `forests` list are as read;
  step 1.5's state test fails on a wrong build; step 1.4 has an identity on the training rows.
- The twelve missing mutants are in Verification with their checks: 1 in push 2, with its case at the
  data door in step 2.6, and 11 in push 2; the rest in push 1.
- The critique's calls. A matrix of logicals of one column is refused as the character matrix is; a
  missing value is judged on the kept rows; NA as a level is refused; a fit stores the factor for a
  forest of levels, and the caller of a swap keeps its own; a state goes by position, said in the
  help; a declared level with no row that the forest lacks is taken, as planned. Each that no ruling
  names is in "Calls made in planning".

What a fit's `bases` changing touches, by search of 37b59ea0: [`packageBartResults`](../../R/bart.R),
which writes it; [`predict.bart`](../../R/generics.R), [`extractForest`](../../R/generics.R) and
[`resolveForestBases`](../../R/generics.R), which read it; 21 lines of 8 test files and
[`recombine`](../../inst/common/recombine.R); the Value of man/bartBT.Rd. No benchmark, vignette,
help example or consumer (bartCause on dbarts-1.0, stan4bart on bartcore) reads one.

Where the critique was off, none of it changing a finding. test-update-scale-rows.R did not arrive
after the slice's rebase: the rebased slice has it, and the tree also gained
test-variance-forest-mixture.R, which names no basis. The first table's 36 objects and 84 refusals
were right and its row for a data frame short of one object. The "10 cells" of its first call are
every cell refused where created; the two matrices of one column are 7 of them. `bartBT(data)` fails
on its signature for any data object and takes one with a response beside it.

Run for this revision, each one short process on the shipped build of 37b59ea0: the critique's seven
probes again; this plan's own probes of the first table, the rows, a swap, a prediction, weights and
the data door's takers, each giving the log it gave before; the 28 cells of a missing value on a
dropped row and the NA level's doors; the spellings the texts advise at a code door; `predict` with a
missing `bases` value under each na.action, and against the stored training fit; `copy()`, a reload
and a sampler made again after a swap that empties a level and after one to a column of zeros; a
value that is not finite on a row `subset` drops, which goes as a missing one does.

Not verified at this revision: the trace was not run again, so every traced count is the earlier
commit's; nothing of the slice is built, so no mutant, no new test and no text was run as written, and
the cells' outcomes under the rule are read off the rule and the tip; `bartBT()` over a data object
was run for one block and no further; the 21 lines that read a fit's `bases` were found by search and
not by the trace, so a read through a variable is not among them; the budget's figures are estimates
and its multipliers are the earlier slices'; no suite, reference build, compare, exact gate, consumer
suite or pair script was run.
