# forest-kind-by-class: what a forest multiplies is decided by the class of its basis, at every door, fixed at creation, and is one numeric column or two levels

Status: PLANNED (dec-B260; dec-B254 for the kind as structure; dec-B281 and dec-B282 for the one size).
Follows the sd-unit slice (forest-sd-unit.md), [forest-defaults-by-kind.md](forest-defaults-by-kind.md)
and push 3 of [written-surface.md](written-surface.md), none of which has landed. Amended 2026-10-07:
the second push also refuses a basis of three or more levels and a basis of several numeric columns, at
every door.

agent: two pushes. Push 1 (the rule): opus implementer for the R code, sonnet for the help once the code is
fixed, opus reviewer told to refute. Push 2 (one size a forest), in two commits: sonnet implementer for
the first, which respells; opus for the R code of the second, which refuses, and sonnet for its tests
and help; opus reviewer for both. The reason for opus on push 1: it moves no arithmetic, and every slip
in it is a forest read as one kind at one door and another at the next, or a coefficient multiplying
another level's rows after a swap, with no message. The tip has both faults today (Context), and the
three reviews of written-surface's second push each found a fault of that shape one step to the side of
the last. On the refusals: a door left out is a model 1.0-0 was ruled not to fit, accepted in silence.
rng: stated per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- Push 1. NEUTRAL for every fit that both builds accept, the three sequences below apart. Measured on a
  prototype of the rule: 36 kinds of object at 7 creation doors, 252 cells, of which 227 have the same
  width, columns and seeded draws on both builds and 3 are the first sequence below; the suite, with the
  new refusals passed by, 16065 results and no seeded draw moved; 11 seeded `bcf()` fits of bartCause.
  POSTERIOR-CHANGING for three sequences, none in a test, a benchmark, a help example or a consumer
  (counted by running the suite on a build that logs every basis, and searched elsewhere):
  1. a logical vector handed to `dbartsData(bases = )`: one 0/1 numeric column before, two levels after,
     which is what the same vector is at every other door. The oracle is the base build's fit of it
     through a `forests` list: identical draws.
  2. `$setForestBasis` on a forest created with a factor or a character vector, given one whose levels
     stand in another order or lack one of them: matched by position before, by name after. The oracle is
     the base build given the same vector as `factor(v, levels = <the levels at creation>)`.
  3. `predict(bases = )` with such a vector: the same, with the same oracle.
- Push 2, first commit. No draw moves: each respelled call is the same block of 0s and 1s. One recorded
  scenario is rewritten, because the second commit refuses its model:
  ["bart2twoforest"](../../benchmarks/R/equivalence.R) of the gaussian equivalence baseline multiplies
  a forest by a factor of three levels. It is given two, recorded afresh on this commit's reference
  build and merged into a copy of the baseline; the other 54 scenarios are bitwise.
- Push 2, second commit. No draw of a model both builds accept moves. Refused where accepted: a factor,
  character or logical basis of three or more levels; a basis of several numeric columns, however
  written; a swap or a prediction that gives a forest either. On push 3's build the suite creates 243
  such forests in 12 files (Context); each is rewritten onto a shape that stays or goes with the test of
  a shape that does not.
window: pre-release, after the sd unit and before the multiplier slice, which makes a factor and a block of
numbers two different priors: until it lands a factor and its indicator columns handed over as numbers
are one fit, so this is the only tip on which every such block in the tests, the benchmarks and bartCause
can be respelled and the respelling proved bit for bit. Serial with any other work in
[`expandForestBasis`](../../R/model.R), [`validateForestBases`](../../R/data.R),
[`setForestBasis`](../../R/dbarts.R), [`resolveForestBases`](../../R/generics.R) or the multi-forest block
of [`resolveSamplerSpec`](../../R/spec.R), and in the file push 3 adds for a basis's reading. No engine,
bridge or header file changes. bartCause's edit lands between the two pushes.
budget: ~3200 lines changed, upper figure 5600. Push 1 ~1100 (upper 2000): R ~280, tinytest ~520, of
which ~400 in one new file, help ~130, records ~170. Push 2 ~2100 (upper 3600): its first commit ~500
(tinytest ~440, benchmarks ~40, help ~20); its second ~1600 (R ~170; tinytest ~1050, of which ~300 in
one new file and ~750 rewritten or removed in 12, test-forest-basis-terms.R most of it; benchmarks ~30;
help ~200; design notes, TODO and the two indexes ~150). The design estimated 900 and planned for 1700
before the refusals; written-surface's two landed pushes ran at 1.8 and 2.7 times their plans.

## Assumed, pending the maintainer

One point of what a user sees is settled here by the design's critic and not by a ruling. It is planned
as written and flagged in the notes for the coordinator.

- A factor handed to `$setForestBasis` or to `predict(bases = )` is matched to the forest's levels by
  name, and a level the forest was not created with is refused. The tip matches by position (Context
  shows what that does); the design had "any number of levels" accepted. If the ruling is "by position,
  documented": step 1.4's alignment is dropped, its texts about a level go, and the second and third
  changed sequences go with them.

## Ruled after this plan was amended

To be worked into the steps at the recheck before push 1 is built; where a step says otherwise, this
section stands.

- dec-B289. Levels are matched by name at a swap and at predict: the assumption above is the ruling.
- dec-B308. A factor basis drops a level no kept row has, without a message, in every spelling and at
  every door that takes a factor; the refusal "... drop it with droplevels()" for a value goes, with its
  text, its pin and the rows of "The rule" and "Refused forms" that carry it. The size of the basis is
  judged after the drop. A numeric column that is zero on every kept row is not covered by the ruling
  and stays as planned.

      forests = list(forest(), forest(basis = arm)), subset = arm != "c"   # created, two levels

## Goal

A forest multiplies nothing, the levels of a factor, or a numeric column, and which of the three follows
from the class of what the caller handed over and from nothing else: a factor, a character vector or a
logical vector is levels; a number is numbers, whatever its values. One object is one kind at every door,
`dbartsData(bases = )` included. The kind and a factor's levels are recorded when the sampler is created
and hold for its life: `$setForestBasis` and `predict(bases = )` take the class the forest was created
with, a factor is matched to the recorded levels by name, and a forest created without a basis takes
none. In 1.0-0 what a forest multiplies has one size: a basis is one numeric column, or a factor, a
character vector or a logical vector of two levels. Three or more levels, and several numeric columns
however they are written, are refused by name at every door, with the form that works: a forest for each
level but the first, or for each column.

## Context

Measured at the tip (4f2e79e6; written-surface pushes 1 and 2 landed, push 3 and forest-defaults-by-kind
not) on the shipped build, R 4.6.1, unless a line says push 3's build (4eeaf03d, read off 2026-10-07):
150 rows; x1, x2; z a 0/1 double with 66 zeros and 84 ones, and the same column as an integer, a logical,
a factor, an ordered factor and a character vector; g a character vector of three levels; dose, age. The
spellings are the tip's. "Draws" is the sum of 10 seeded train draws of a two-forest sampler, 15 trees;
the same number is the same fit.

- The suite. At the tip 231 files, of which 4 exit off a reference build and 3 ask for more than two
  threads; the other 224 give 16065 results, none failing, 102 seconds in one process. On push 3's
  build 233 files, 230 without the three, 17121 results, none failing.
- The kind at creation, today. The forest doors are: a `forest()` term of a formula at `dbarts()` and at
  `bart()`; a `forests` list on a formula, its basis written as code or handed over as a value; the same
  list on the matrix interface and in `dbartsSpec()`. They read one object alike in every case tried. The
  data door, `dbartsData(bases = )`, reads it otherwise.

  | handed over | the forest doors | the data door |
  |---|---|---|
  | a factor; an ordered factor; `I()` of one | levels: 2 columns (66, 84), draws 2103.634934 | refused: "'bases' must be numeric or logical" |
  | a character vector | levels, the same fit | refused, the same |
  | a logical vector; `z == 1`; `I()` of one | levels, the same fit | one numeric 0/1 column: draws 2104.118263, another model |
  | a 0/1 integer; a 0/1 double; `I(z)`; `matrix(z)`; `as.numeric()` of a logical; a time series | numbers, one column, draws 2104.118263 | the same |
  | two values that are not 0 and 1 (-0.5, 0.5; 1, 2); three values; dose | numbers, one column | the same |
  | `cbind(1 - z, z)`; `model.matrix(~ g - 1)` | numbers, 2 and 3 columns; draws 2103.634934, the factor's own | the same |
  | `model.matrix(~ g)[, -1]`, a level omitted; `cbind(dose, age)` | numbers, 2 columns | the same |
  | a logical matrix | refused: "a 'basis' must be numeric, a factor, or a character vector" | numbers, 2 columns |
  | a logical matrix of one column | refused, the same | numbers, one column |
  | a character matrix of one column | levels | refused |
  | a data frame of one column or two | refused, two texts by door | refused |
  | a Date, a list, a complex | refused | refused |
  | a factor with a level no row has | refused: "... drop it with droplevels()" | refused, as not numeric |
  | a factor of one level | refused | refused |
  | a number with a missing value | refused: "a 'basis' cannot be NA" | refused: "'bases' values must all be finite" |
  | `rep(1, n)` | numbers, one constant column | the same |
  | `rep(0, n)` | refused: "a 'basis' column of all zeros contributes nothing to a forest" | accepted: a forest multiplied by zero |

  So a 0/1 number is one numeric column at every door today, and is not the factor's fit; dec-B260
  changes nothing for it. What the tip decides by something other than class is the data door alone.
- The same on push 3's build. The data door is as above. At the forest doors a basis written as code is
  the columns of its model frame: `arm`, a factor of three levels, 3 columns (`arma`, `armb`, `armc`);
  `dose + age` 2; `ns(dose, 3)` 3; `poly(dose, 2)` 2; `1 + dose` 2 (`(Intercept)`, `dose`); and one
  column each for `dose`, `0 + dose`, `scale(age)`, `poly(dose, 1)` and `dose:age`. A matrix written in
  place in a `forests` list, `forest(basis = cbind(dose, age))`, is read as code and refused by the
  grammar ("the columns of a basis are separated by '+'"); a matrix is a value only when a variable
  holds it. `forest(x1, basis = arm == "b") + forest(x1, basis = arm == "c")` beside a forest with no
  basis is created: three forests, widths 0, 2, 2, labels `arm == "b"` and `arm == "c"`. So is
  `forest(x1, basis = dose) + forest(x1, basis = age)`. A code basis travels through the data object as
  one column of row numbers and is built on the rows kept afterwards; its record for new rows holds the
  levels kept.
- Rows. Under a `subset` that leaves one level of g with no row, and under a missing response on those
  rows:

  | basis | door | result |
  |---|---|---|
  | a factor column named in a formula term | formula | refused, the droplevels text |
  | `factor(g)`, or the character column, in a formula term | formula | 2 columns, the level dropped; draws 1680.610731 |
  | the factor column, or the character column, written in a `forests` list; the factor handed over as a value; on the matrix interface; its indicator columns at the data door | list, matrix, data | 3 columns, one of them all zeros, accepted; draws 1669.160511 |
  | a level emptied by a missing response | formula and list | 3 columns, one all zeros, accepted |
  | a number left all zero by `subset`, handed over in a list | list | accepted: one column of zeros |
  | the same written as code in a formula term | formula | refused |
  | a logical left all `TRUE`, handed over in a list; through the data door | list; data | 2 columns, one all zeros; one constant column |

  Push 3 of written-surface drops an emptied level for a basis written as code, at both doors (read off
  its build: `arm` under a subset that leaves two levels has 2 columns). A value is this slice's.
- After creation. `$setForestBasis` takes any class on any forest: 8 kinds of forest by 16 kinds of value
  were tried, and every pair is accepted except a factor of one level, a column of zeros, a logical
  matrix, and one numeric column on a held forest (dec-A171). A forest created with no basis takes one,
  of any kind; a factor forest takes dose, or `cbind(dose, age)`; a numeric forest takes a factor; the
  width follows the value. On push 3's build the same: a forest created on dose took two columns, then
  one, then a factor of three levels.
- A factor's levels are matched by position. A held two-level character forest (control, treated):
  mean absolute term 0.000 on the control rows and 0.129 on the treated; after a swap to the same vector
  as a factor with levels (treated, control), 0.133 and 0.000. A drawn three-level forest with
  coefficients -0.785, -0.158, 1.436: swapped to the vector with one level redrawn into another it has 2
  columns and the rows of the third level take the second coefficient; swapped back, the third
  coefficient is 1. A vector left with one value is refused.
- `predict(bases = )` expands what it is given by its own levels and checks the width ("'bases' gives
  forest 2 2 columns; its amplitudes take 1"). On a fit created with `factor(z)`: the same factor with
  its levels reversed is accepted and predicts 1.267 where the aligned one predicts 1.433; a two-level
  character vector unrelated to z is accepted; `cbind(dose, age)` is accepted (mean prediction 26.7). On
  a fit created with a 0/1 number a factor of one level is accepted as a column of ones.
- What is recorded. A data object has `bases`, bare numeric matrices, and no record of levels; the
  forests' configuration has no kind. A fit carries `bases` and, for a formula term, the term with its
  levels ([`packageBartResults`](../../R/bart.R), [`replayForestBasis`](../../R/model.R)).
- What the suite hands over, counted at the tip by running it on a build that logs every basis. 776
  models of several forests are created in 53 files. Their bases: a factor 281 times; a 0/1 integer or
  double vector 403 times in 15 files (239 in test-forest-predictors.R, 127 in test-formula-terms.R);
  another numeric vector 28; a numeric matrix at a forest door 32, of which 7 are a factor's indicator
  columns; through the data door 139, of which 91 are indicator columns, 13 indicator columns times a
  constant, 33 a constant column and 2 other. Of the 98 blocks of indicator columns, 61 are in
  test-bcf-family.R and 14 in test-bcf-creation.R. No character vector, logical vector or ordered factor
  reaches a creation.
- What the suite creates on push 3's build, every creation traced. 1186 models of several forests.
  Per forest with a basis: two indicator columns 434; one numeric column 374; one 0/1 column 341; one
  constant column 55; three indicator columns 123; two other numeric columns 74; two columns with one
  nonzero entry a row that are not 0s and 1s 25; three numeric columns 19; four 2. The last five rows,
  243 forests, are what dec-B281 and dec-B282 refuse whatever class they came in as: 160 in
  test-forest-basis-terms.R (86 of three levels; 74 of two, three or four numeric columns), 50 in
  test-bcf-family.R (26 three-column blocks and 24 scaled two-column blocks through the data door), 12
  in test-forest-capture.R, 5 in test-forest-labels.R, 4 each in test-forest-basis-r5.R and
  test-multiforest-leaf-prior-writer.R, 3 in test-formula-terms.R, 2 in test-sampler-residuals.R, one
  each in test-argument-surface.R, test-bcf-creation.R and test-predict-blend.R. After creation, values
  handed over at the forest doors: swaps with two or three numeric columns 26, in 7 files; predictions
  with two to four numeric columns 123 and with a factor of three levels 15, nearly all in
  test-forest-basis-terms.R and test-predict-blend.R. Many of the two-column swaps and predictions are
  indicator columns given to a factor forest, which push 1 already respells.
- `$setForestBasis` is named on 53 lines of 11 test files at the tip and reaches a sampler 29 times in 7.
  Under push 1's rule 20 of those meet a refusal: 13 put numbers on a factor forest
  (test-forest-basis-r5.R 4, test-forest-arguments.R 4, test-bcf.R 2, and one each in
  test-bcf-mutation-pins.R, test-multi-forest-seam.R and test-multiforest-leaf-prior-writer.R), two of
  which are today's tests of dec-A171's refusal, and 7 give a basis to the forest with none
  (test-forest-basis-r5.R 4, and one each in test-bcf-r5-surface.R, test-forest-arguments.R and
  test-multiforest-leaf-prior-writer.R). Three more that are refused today for their length would meet
  the rule first. `predict` is given `bases` 64 times in 6 files, always numbers; 14 of them, all in
  test-predict-blend.R, hand indicator columns to a factor forest.
- The prototype against the suite, at the tip. As push 1's rule stands, 6 files stop at their first
  refused call (test-bcf.R, test-forest-arguments.R, test-forest-basis-r5.R, test-multi-forest-seam.R,
  test-multiforest-leaf-prior-writer.R, test-predict-blend.R: 2895 lines, 545 results at the tip). With
  the refusals logged and passed by, 16065 results run and 8 fail: one pin of the droplevels text, and 7
  pins of `data@bases` against a bare matrix that fail only because the prototype keeps the levels as an
  attribute of the block, which the plan does not.
- Elsewhere. benchmarks/R: every basis at creation is a factor of two levels, but for the `bcf` row of
  `composition-matrix.R`, a 0/1 number, and for the scenario
  ["bart2twoforest"](../../benchmarks/R/equivalence.R) of `equivalence.R`, a factor of three levels,
  which is one of the 55 recorded scenarios of the gaussian baseline. `bcf-equivalence.R` swaps
  `cbind(1 - z2, z2)` onto a factor forest once, in the scenario `set_treatment`. Help: the example of
  man/bart.Rd multiplies by a 0/1 integer and calls it the Bayesian causal forest shape; on push 3's
  build man/forest.Rd and man/bart.Rd describe `basis = a + b` as two columns of one forest, with
  `1 + dose`, `poly()` and `ns()` beside it. The vignettes fit no model of several forests.
- Consumers. bartCause's `bcf()` (dbarts-1.0, 6c1bff9) hands `cbind(1 - z, z)` to `dbartsData(bases = )`
  under its own `subset`, having refused a treatment that the kept rows leave with one value; its
  treatment forest takes `amplitude = fixed()` when `update.b = FALSE`, which on two columns as the
  second forest holds (0, 1). Two of its test files build the same sampler by hand from the same block
  (tests/testthat/test-14-bcf.R in three places, test-03-responseFit.R in one). Run for this plan: 11
  seeded `bcf()` fits (the default, a binary and a logistic response, `subset`, a logical treatment,
  weights with an offset, stated sds, each coefficient held and both, another coefficient variance) are
  identical between bartCause as it stands on the tip, bartCause with
  `basis <- factor(as.integer(z), levels = 0:1)` on the prototype, and bartCause as it stands on the
  prototype; its two bcf test files pass on all three, 156 and 47 expectations. stan4bart (bartcore,
  a9d081b), treatSens (dbarts-1.0, aecec71) and bairrtt (main, 3f57f61) declare one forest with
  `forest(n.trees = )` and no basis (searched again 2026-10-07).

## The rule

At creation, at every door. The bracket is the push that makes the row so.

| handed over | kind | coefficients |
|---|---|---|
| nothing | none | one, on the forest itself |
| a factor, ordered or not; a character vector; a logical vector; `I()` of any of them: of two levels [1] | levels | one per level, in the factor's own order of levels, the sorted values of a character vector, `FALSE` then `TRUE` |
| an integer or double vector, or a matrix of one column, whatever its values; `I()` of one; one numeric term of one column [1] | numbers | one |
| a factor, a character vector or a logical vector of three or more levels | levels until push 2; then refused by name [2] | |
| a numeric matrix of several columns; several numeric terms; one term of several columns, a spline or a polynomial among them; `1 +` beside a column | numbers until push 2; then refused by name [2] | |
| a logical or character matrix; a data frame; anything else | refused by name [1] | |

- The class is looked at in one place, [`expandForestBasis`](../../R/model.R), which every creation
  door reaches. The data door hands it what it is given.
- The content is never looked at to decide a kind: a 0/1 number is a number, and a factor's indicator
  columns handed over as numbers are several numeric columns, refused as such.
- A basis written as code is the class its code gives in the fit's model frame, as push 3 builds it;
  the same vector handed over as a value is that class too.
- Empty at creation. A level no kept row has, and a numeric column that is zero on every kept row, are
  refused for a value, once, where the data object has its rows: after `subset` and the na.action, at
  every door. Code has dropped its empty levels by then (push 3).
- The size is judged after the rows: a factor of three levels written as code that `subset` leaves
  with two is a basis of two levels; handed over as a value it meets the empty-level text first, and
  with `droplevels()` is two levels.
- A hold changes nothing here: a basis of three levels or of several columns is refused drawn or held,
  in the same words.
- The record. A factor's levels are kept with the data, per forest; a forest's kind is read from the
  data by one function and from nowhere else.

After creation, at `$setForestBasis` and at `predict(bases = )`, by one function:

| forest created as | takes | refuses |
|---|---|---|
| none | nothing | any basis |
| levels | a factor, a character vector or a logical vector, matched to the forest's levels by name; a level with no row is a column of zeros | numbers, indicator columns included; a level the forest does not have |
| numbers | a numeric vector or a matrix of one column; a column of zeros | a factor, a character vector, a logical vector; from push 2, several columns |

- A levels forest has the same columns for life: its width cannot change, and coefficient j multiplies
  level j's rows after any swap.
- From push 2 no forest changes width, held or drawn: a numeric forest is one column. Between the
  pushes a drawn numeric forest takes another width as the tip does, and a held forest refuses one by
  forest-defaults-by-kind's interim text.
- A basis written as code and rebuilt at new rows from the fit's own record is untouched.
- Nothing of this moves a draw of a sampler that is swapped to a block the tip would have built the same
  way.

Accepted where refused: a factor, an ordered factor and a character vector at the data door; a factor of
one value at a swap and at predict, where its level is one of the forest's (every new row treated); a
column of zeros on a swap of a numeric forest. Refused where accepted: a logical matrix, a column of
zeros and an empty level at the data door; an emptied level or column of a value under `subset` or a
missing response, at the list and the matrix doors; a character matrix of one column; every cell of the
second table's last column; and, with push 2, three or more levels and several numeric columns.

Before and after, per door.

| door | before | after push 1 | after push 2 |
|---|---|---|---|
| a `forest()` term at `bart()` and `dbarts()`; a `forests` list on a formula, on a matrix, in `dbartsSpec()` | by class | by class; the same fits | one numeric column or two levels; anything larger refused |
| `dbartsData(bases = )` | numbers, a logical vector among them; a factor refused | by class | the same sizes; two numeric columns, a factor's indicators among them, refused |
| `$setForestBasis` | any class, any forest, levels by position | the table | the table; no forest changes width |
| `predict(bases = )` | any class at the fit's width, levels by position | the table | the table |
| `predict` from a formula term's stored record; `fitted`; `extract` | | unchanged | unchanged |
| `$getLeafPrior`, `$setLeafPrior`, `print`, `show` | say nothing of a kind | unchanged | unchanged; the reader's entry and the printed block are the multiplier slice's |
| `copy()`, a reload, `new("dbartsSampler", control, model, data)` | carry the block | carry the block and its levels: the same kind | a data object edited to a refused size is refused |
| a state stored or installed | carries no basis | unchanged | unchanged |
| the engine | a coefficient for each column of any block | unchanged | unchanged: nothing is narrowed there |

## Refused forms, with their texts

Base R's style, as in written-surface. `<arg>` is `basis` at a forest door and at `$setForestBasis`,
`bases` at the data door and at `predict`. `<f>` is `forest 2`, or `forest 2 ("dose")` where the forest
has a label. The push that adds each is in brackets.

    [1] '<arg>' must be a factor, a character or logical vector, or numeric; a logical matrix is none of these: write it as numbers (<arg> + 0), or give one logical column            (a character matrix: "give one character column")
    [1] '<arg>' must be a factor, a character or logical vector, or numeric, not a data frame; give one of its columns            (a Date, a list, a complex: named by its class)
    [1] <f> has a basis level with no observations: "c"; drop it with droplevels()
    [1] <f> has a basis column of all zeros, which contributes nothing; drop it
    [1] <f> multiplies nothing: it was created without a basis, and what a forest multiplies is fixed when the sampler is created; make a new sampler to give it one
    [1] <f> multiplies a factor, so '<arg>' must be a factor, a character vector or a logical vector; numbers are another model, an indicator matrix included: write factor(), or make a new sampler to multiply numbers
    [1] <f> multiplies numbers, so '<arg>' must be a numeric vector, not a logical vector; write as.numeric(<arg>) for a 0/1 column, or make a new sampler to multiply a factor            (a factor, a character vector)
    [1] <f> has no level "s"; its levels are "p", "q", fixed when the sampler was created
    [1] 'basis.levels' gives forest 2 3 levels and its basis has 2 columns
    [2] <f> has a basis of 3 levels ("a", "b", "c"); in this version a factor, a character or a logical basis has two levels. Give each level but the first a forest of its own, as forest(x1, basis = arm == "b") + forest(x1, basis = arm == "c")
    [2] <f> has a basis of 2 numeric columns (dose, age); in this version a basis is one numeric column, or a factor, a character or a logical vector of two levels. Give each column a forest of its own, as forest(x1, basis = dose) + forest(x1, basis = age)
    [2] <f> has a basis of 2 numeric columns ((Intercept), dose): '1 +' in a basis adds a constant column, and in this version a basis is one numeric column. Write the column alone, as basis = dose
    [2] <f> multiplies one numeric column, and '<arg>' has 2; in this version a basis is one numeric column, or a factor, a character or a logical vector of two levels

In the first two [2] texts the example is written with the caller's own words where they are code: the
forest's predictors as written, the factor's name, each term that is one column. Where the basis is a
value, or one term of several columns, the example is the fixed one above. The last [2] text is the
swap's and the prediction's.

The ninth is the data object's validity, for an object edited by hand. Replaced by push 1: "'bases' must
be numeric or logical", "a 'basis' must be numeric, a factor, or a character vector", "a 'basis' factor
level with no observations contributes nothing to a forest" and "a 'basis' column of all zeros
contributes nothing to a forest". Replaced by push 2: forest-defaults-by-kind's interim texts for a held
basis of three or more columns and for a held forest's width at a swap; "'bases' gives forest 2 2
columns; its amplitudes take 1" where the forest is numeric; push 3's refusal of a swap whose columns
stand in another order, which has no case left. Kept as they are: "a 'basis' cannot be NA", "a 'basis'
factor must have at least two levels" at creation, the texts about a basis's length, the
at-least-two-forests text, dec-A171's refusal of a held single numeric column at creation,
forest-defaults-by-kind's two interim texts by position, and at `predict` the text for a forest given a
basis it has none of. A text of push 3 for a basis written as code (a logical matrix as a term, a
logical with one value, a factor beside other terms, an operator it does not take) comes first where
both apply.

## Constraints

- Every fit that both builds accept, the three sequences of the `rng:` line apart, draws what it draws
  on the base build. The gaussian baseline is re-recorded for the one scenario push 2 rewrites and for
  nothing else; no snapshot file is regenerated.
- No change to src/, to the flat C header, to the stored state or to what the bridge reads. The engine
  keeps a coefficient for each column of whatever block it is handed, and is not told what 1.0-0
  refuses.
- The class is read in one function at creation and the table applied by one function afterwards; the
  size is refused by one function; no door has a copy of any of the three.
- `data@bases[[f]]` stays a bare numeric matrix, with the names push 3 gives it and no other attribute.
- No law changes: until the multiplier slice a factor of two levels and a 0/1 number keep the fits they
  have.
- What push 3 fixed for later stays fixed (dec-B272): how the columns of a basis of several terms are
  named, and every refusal of a spelling that could later state a size for each column. A basis of
  several terms is still read and named, and then refused.
- The reader, `extract` and `print` are untouched. The unit and the default of an sd are untouched.
- Nothing of `updateBasisScale` rides along.
- Push 1 adds no test on a factor of three or more levels or on several numeric columns, the rows of
  its loop apart that push 2 turns.
- Each push leaves the help saying what the code does. Base R calls stay within DESCRIPTION's R floor.

## Pushes

1. The rule. Changes the three sequences; the design note lands here. After it the data door takes a
   factor, which is what lets bartCause hand one over.
2. One size a forest. Two commits, gated together.
   - First: every call in the tests, the benchmarks and the help that hands a two-level factor's
     indicator columns over as numbers is written as the factor it stands for, bit for bit; the one
     baseline scenario on three levels is rewritten and recorded.
   - Second: a basis of three or more levels, and of several numeric columns, is refused at every door
     (dec-B281, dec-B282), and the tests of those shapes are rewritten or go.

bartCause's line lands on its own branch after push 1 and before push 2: its `cbind(1 - z, z)` through
the data door is two numeric columns, which push 2's second commit refuses. No tip breaks it: push 1's
tip takes both spellings, and push 2 is not pushed until bartCause's suite has passed on it with the
factor.

Push 1 alone leaves a coherent tip: its own respelling is the 20 swaps and 14 predictions it refuses,
which it carries. The two commits of push 2 are cut where the review changes: the first is read against
"the same draws" and nothing else, the second against "which door did it miss" and "what coverage went
with a deleted test".

## Steps

"Fails today" is what the base build does where the test expects otherwise. New names are for the
reader's sake. Calls are in push 3's spelling; the fixture is Context's, with `arm` a factor of three
levels beside g. A value is always held in a variable: push 3 reads a call written in place as code.

### Push 1: the rule

1.1 One reader of a class. [`expandForestBasis`](../../R/model.R) is the only place a value's class is
    read at creation, and returns the block with, for levels, the levels it was expanded from. It
    refuses the classes of "The rule" with the first two texts and no longer looks for an empty level or
    column (step 1.3). [`validateForestBases`](../../R/data.R) hands it any entry that is not already a
    block it expanded, so the data door has no rule of its own; a block that arrives expanded keeps its
    levels and is not read again as numbers. Tests, a new file test-forest-kind.R, built as one loop:
    each object of Context's first table is taken through the seven creation doors, and the test asserts
    against the FIRST door, so that no door can be right alone. Per object: the kind, the width, the
    levels, `data@bases[[2]]` `identical()` across doors, and the 10 seeded draws `identical()` across
    the sampler doors; a refused object is refused at every door with one text but for `<arg>` (fails
    today at the data door for a factor, an ordered factor, a character vector, a logical vector, a
    logical matrix and `rep(0, n)`). The loop's table has one column for what push 2 expects, and the
    objects of three levels and of several columns sit in it as created here. Pinned beside it, so that
    a reading by content fails: a 0/1 integer, a 0/1 double, `I(z)` and `matrix(z)` are numbers of one
    column and do NOT draw what `factor(z)` draws; `cbind(1 - z, z)` held in a variable is numbers of
    two columns. The same column written as code and handed over as a value is one kind and one fit.
1.2 The record. A slot `basis.levels` on [`dbartsData`](../../R/A_class.R), beside `response.levels`:
    `NULL`, or a list with one entry per forest, `NULL` for a forest with no basis or with numbers and
    the levels for one with a factor. It is written wherever `bases` is: where
    [`dbartsData`](../../R/data.R) finishes its rows, with the levels carried across every restriction
    of rows ([`restrictBasesToRows`](../../R/data.R), [`alignForestBasisToSubset`](../../R/model.R)) and
    taken off the block before it is stored; and where a declaration's bases are put on a data object
    ([`resolveSamplerSpec`](../../R/spec.R)), replacing what the object had. A basis written as code
    reaches the data object as push 3's column of row numbers and is built afterwards: its levels are
    recorded when the block is, from the levels push 3's record keeps, and are those its stored `terms`
    rebuild at new rows. The slot is read through an accessor that answers `NULL` for an object saved
    without it, as [`dataRowNames`](../../R/data.R) does. One function, `forestBasisKind(data, f)`,
    answers "none", "levels" or "numeric" and is the only reader of a kind after creation. The class's
    validity refuses recorded levels whose count is not the block's width.
    [`packageBartResults`](../../R/bart.R) puts the levels on a fit beside `bases`. Tests: the slot at
    each of the seven doors for a factor, a character and a logical vector, a number and no basis; the
    same under `subset` at the formula, list, matrix and data doors, and after a dropped missing
    response (a record lost with the rows makes the forest numeric there); `dbartsSpec()` over a data
    object that carried a numeric block for that forest records the declaration's levels;
    `attributes(data@bases[[2]])` holds the dimensions and names and nothing else; an object with the
    slot removed creates, runs and reads as numbers; an object with levels edited to another count is
    refused with the ninth text.
1.3 Empty at creation. One check where the data object has its final rows, and the same check where
    `dbartsSpec()` installs a declaration: a level with no row, or a column of all zeros, is refused
    with the third and fourth texts. Tests: a factor value with a level never used; one whose level
    `subset` empties; one emptied by a missing response; a number left all zero by `subset`;
    `rep(0, n)`: each at the list, matrix, data and `dbartsSpec` doors, one text (fails today: all but
    the first and fifth are created at the list and matrix doors, and the fifth at the data door); and
    beside the value the same column written as code under the same `subset`, which push 3 creates with
    the level dropped, pinned so that the two readings stay told apart.
1.4 After creation. One function (`conformBasis`: a value, the forest's kind, its levels, its width)
    returns the block for the engine or refuses by "The rule";
    [`setForestBasis`](../../R/dbarts.R) and [`resolveForestBases`](../../R/generics.R) both call it,
    the sampler reading kind and levels from its data object and a fit from what it carries. In
    `$setForestBasis` nothing is stored before it returns. A one-sided formula is evaluated first, as
    today. In this push it adds no rule about a numeric forest's width: a drawn one takes what the tip
    takes, and a held forest's width is forest-defaults-by-kind's interim refusal. Tests, every factor
    in them of two levels:
    - The table: each kind of forest (none; a factor, a character and a logical forest; a 0/1 number;
      dose) by each kind of value at `$setForestBasis` and at `predict`, each cell accepted or refused
      with its text (fails today: all but four kinds of value are accepted everywhere). After every
      refused swap the data object, `getForestAmplitudes()` and the next 5 sweeps are `identical()` to
      an untouched twin's.
    - Which model, by name. A held two-level character forest swapped to the same vector with its
      levels reversed keeps its term on the treated rows (fails today: it moves to the control rows). A
      drawn two-level forest swapped to a vector left with one value keeps 2 columns, one of them zero,
      and its coefficients; swapped back, the two coefficients are `identical()` to those before (fails
      today: refused). A level the forest lacks is refused.
    - `predict`: a factor with its levels reversed predicts `identical()` to the aligned one (fails
      today: 1.267 against 1.433); every row at one level, written `rep("1", n)`, predicts what
      `factor(rep(1, n), levels = 0:1)` does today; a renamed level is refused (fails today: taken);
      numbers on a factor forest and a factor on a numeric forest are refused (fail today: taken where
      the width fits).
    - A column of zeros swapped onto a numeric forest is accepted and the run goes on (fails today:
      refused).
    - With the forest with no basis second in a list, as forest-defaults-by-kind accepts it, that
      forest refuses a basis and the first takes one: the kind is not a position.
1.5 The same sampler again. Tests: after `copy()`, after `saveRDS`, `readRDS` and a first use, and on
    `new("dbartsSampler", control, model, data)` from the sampler's own three, the kind and the levels
    are those of the original, the by-name swap of step 1.4 gives the same block, and a refused swap is
    refused; a state taken from a sampler created on the same factor with its levels in another order
    installs, and the recipient's levels are its own.
1.6 Respell what the rule refuses. The 20 swaps and the 14 predictions of "Context", recounted on the
    landed tip: a numeric block swapped onto, or predicted for, a factor forest becomes the factor
    (`cbind(1 - z2, z2)` becomes `factor(z2)`); a basis given to the forest with none becomes a pin of
    the fifth text where the test is about that, and otherwise the forest is created with a basis. The
    three swaps refused today for their length are respelled so that the length is still what is tested.
    The pin of
    ["factor level with no observations contributes nothing"](../../inst/tinytest/test-forest-basis-r5.R)
    takes the third text. `bcf-equivalence.R`'s swap becomes `factor(z2)`; its scenario is unchanged to
    the bit. Run the suite first and repair what it shows.
1.7 Help and records. man/forest.Rd, the `basis` item as push 3 leaves it: which class is which kind;
    that a 0/1 number is one column with one coefficient and not two levels, and `factor(z)` is how
    levels are written; that a value's empty level or column is refused and code's is dropped.
    man/dbartsData.Rd, `bases`: the same classes as `forest(basis = )`.
    [`dbartsSampler$setForestBasis`](../../man/dbartsSampler-class.Rd) and the method's docstring: the
    table, and that levels are matched by name. man/bartBT.Rd, the `bases` item of `predict`: the
    table. docs/design/forest-kind-by-class.md with its index row: the rule, the record, the three
    changed sequences with their oracles, what a sampler saved before reads as.
    docs/design/public-surface.md and docs/architecture.md where they describe a basis and the data
    object's slots. TODO: `forest-prior-args` names this push landed.
1.8 Mutations (Verification): apply each, install, run the named test, record the failing count, revert,
    `touch` the file.

### Between the pushes: bartCause

B.1 Its own commit on dbarts-1.0, after push 1 has landed on bartcore. R/bcf.R:
    `basis <- cbind(1 - z, z)` becomes `basis <- factor(as.integer(z), levels = 0:1)`; `as.integer`
    because `bcf()` takes a logical treatment, which `factor(z, levels = 0:1)` would turn into missing
    values. The comment above it loses its sentence about the column order of `cbind`, the order now
    being that of the levels 0, 1. Its hand-built samplers hand the factor over in the same commit, so
    that they stay the comparator of `bcf()` under the multiplier slice: tests/testthat/test-14-bcf.R in
    its three places (the builder's `bases = `, the aligned block, and the block already cut to the kept
    rows, which keeps its refusal for its length) and test-03-responseFit.R in its one. Its pin of the
    second column of the aligned block against the kept treatment stands: a factor's block is still
    stored as its indicator columns. No other file there reads a basis.
B.2 Gate: bartCause's suite whole, on push 1's tip, before and after the commit, none failing; the 11
    seeded `bcf()` fits of Context identical before and after. Not on any earlier tip: before push 1
    the data door refuses a factor.

### Push 2: one size a forest

First commit, the blocks written by hand.

2.1 Tests. Each creation from a two-level factor's indicator columns handed over as numbers (at the
    tip 98 blocks, 91 through the data door and 7 at a forest door, in 10 files, some of them of three
    columns; recounted on the landed tip) is read at its fixture. Where a block of two columns stands
    for a factor (the treatment blocks of test-bcf-family.R, test-bcf-creation.R, test-bcf-loglik.R,
    test-bcf-reporting.R, test-forest-basis-subset.R and the others the log names) the fixture hands
    the factor over, and a pin of `data@bases` against the block is kept, against the block. Not
    touched in this commit: a block of three columns, a block times a constant and a constant column;
    and the creations on a 0/1 number, which are one numeric column before and after. Every edited
    file's results are those of the base build, assertion for assertion.
2.2 Benchmarks and help. `composition-matrix.R`'s `bcf` row is written `basis = factor(z)`; man/bart.Rd's
    example of the Bayesian causal forest shape is written with `factor(z)`, and the line of its
    section "Formula Terms" that glosses a numeric z says one coefficient, not a treatment. These two
    change the model the text fits, which nothing pins.
2.3 The baseline scenario. ["bart2twoforest"](../../benchmarks/R/equivalence.R) draws its factor from
    two values where it drew from three, with nothing else of its text changed, and its comments say
    two levels; its response moves with the factor, the scenario's data being drawn after it from one
    seed of its own, and no other scenario's data does. Record it alone on this commit's
    reference build (`EQUIVALENCE_SCENARIOS=bart2twoforest`, `EQUIVALENCE_CORES=1`), merge it into a
    copy of `equivalence-1b7d730c.rds` in its scenario order, and name the file after this commit. The
    MANIFEST row says what moved and why (the scenario's own data and model, not the sampler: the other
    54 are bitwise on the same build), and names as oracle what gates a two-level factor forest, the
    BCF exact gates in `quick` on this build.

Second commit, the refusal.

2.4 One function refuses a size (`refuseBasisSize`: the built block, its levels or none, the forest's
    label, the basis's text and its terms where it is code, the forest's own predictors as written).
    More than two levels: the first [2] text. More than one numeric column: the second, or the third
    where the first column is the constant `1 +` adds. It is called wherever a creation's block is
    final: beside step 1.3's check where the data object has its rows and where `dbartsSpec()` installs
    a declaration, after that check, so that a value's empty level is met first; where push 3 builds a
    code basis on the rows kept, which has dropped its empty levels; and in the data class's validity
    beside step 1.2's, for an object edited by hand. It is called before anything is read of a hold, so
    forest-defaults-by-kind's interim refusal of a held basis of three or more columns can no longer be
    reached and goes.
2.5 After creation. `conformBasis` refuses a numeric value of several columns for a numeric forest, with
    the last [2] text, at `$setForestBasis` and at `predict`; with it no forest changes width, and
    forest-defaults-by-kind's interim refusal of a held forest's width at a swap goes, as does push 3's
    refusal of a swap whose columns stand in another order.
    [`validateForestSd`](../../R/model.R)'s text for a longer `sd` is reread: it no longer speaks of
    columns since forest-sd-unit.
2.6 Tests, a new file test-forest-basis-one-size.R.
    - The loop. Each refused object: a factor, an ordered factor and a character vector of three
      levels; `dose + age`; `ns(dose, 3)`; `poly(dose, 2)`; `1 + dose`; and, held in variables, a
      matrix of two numeric columns, `cbind(1 - z, z)`, `model.matrix(~ g - 1)` and
      `model.matrix(~ g)[, -1]`. Each at every door it can be written at: code as a formula term at
      `dbarts` and at `bart`, in a list on a formula and in `dbartsSpec`; a value in a list on a formula
      and on a matrix, in `dbartsSpec` and at `dbartsData(bases = )`. Each drawn and with
      `amplitude = fixed()`: one text, the basis's, asserted against the first door's but for `<arg>`
      and `<f>` (fail today: created). The message's example, taken from the text and run as written
      where it is the caller's own words, creates a sampler.
    - What stays. A factor, a character and a logical vector of two levels; `arm == "b"`; one numeric
      column as a vector, as a matrix of one column, and as `poly(dose, 1)`, `scale(age)`, `dose:age`,
      `I(dose * age)`, `0 + dose` and `dose - 1`: each created at every door, with the block step 1.1
      pins.
    - After the rows. `arm` written as code under a `subset` that leaves two of its levels is created
      with those two. The same factor handed over as a value under that `subset` is refused with step
      1.3's text, and with `droplevels()` is created. Under a `subset` that keeps all three it is
      refused with the first [2] text at both.
    - The forms the messages show. `y ~ forest(x1 + x2) + forest(x1, basis = arm == "b") +
      forest(x1, basis = arm == "c")` is three forests with those labels; `forest(x1, basis = dose) +
      forest(x1, basis = age)` beside a forest with no basis is three forests.
    - After creation. A numeric forest swapped to a matrix of two columns is refused, the sampler
      `identical()` to an untouched twin's afterwards and over the next 5 sweeps; a matrix of one column
      is taken. `predict(bases = )` with two columns for a numeric forest is refused in the same words.
      forest-defaults-by-kind's three swaps of a held forest to another width take the text.
    - Edited by hand. A data object whose `bases` entry is set to two numeric columns, and one whose
      recorded levels are set to three with a block to match, are refused at `dbarts()` and at
      `new("dbartsSampler", control, model, data)`.
    - The names, for later (dec-B272). The shapes push 3 pins against
      `names(coef(lm(y ~ 0 + <basis>)))` are read from the function that builds a code basis's columns,
      called on the model frame with no sampler, so the names of several columns stay pinned though no
      model has them; the refusal's text quotes those names.
2.7 Rewrite and remove, by what "Context" traced on push 3's build, recounted on the landed tip; run
    the suite first and repair what it shows.
    - test-forest-basis-terms.R, 160 creations and nearly all the predictions. Its fixture on g, three
      levels, moves to a column of two. Where a block tests the grammar or the names of several terms
      it reads the builder, as step 2.6's last test does, and fits nothing. Where it tests what a basis
      of several columns does in a fit or at new rows (the two-column rows of its blocks on columns,
      on rows and on new rows; the changed texts `dose + age` and `1 + dose`) it becomes that shape's
      refusal or goes, and the one-column rows stand. `scale()`, `poly(dose, 1)` and `I()` at new rows
      stand.
    - test-bcf-family.R, 50 creations through the data door: its arms on three groups and on a block
      times a constant go. They pin the tip's row norm on shapes that no longer exist; what they pin of
      a two-level block is held by the two-column arms the first commit made factors. Its 54 constant
      columns are one numeric column and stay until the multiplier slice.
    - test-forest-capture.R 12 and test-forest-labels.R 5 creations, with 10 swaps: a basis of
      `dose + age` used only as some basis becomes `I(dose + age)` or one of its columns, the label and
      the capture under test unchanged; the tests of a swap's column order go with the refusal.
    - test-forest-basis-r5.R 4, test-multiforest-leaf-prior-writer.R 4, test-formula-terms.R 3,
      test-sampler-residuals.R 2, and one each in test-argument-surface.R, test-bcf-creation.R and
      test-predict-blend.R: a factor of three levels becomes one of two where the test is not about the
      third, and the test goes where it is. test-predict-blend.R's predictions with several columns the
      same.
    - test-forest-kind.R: the loop's column for push 2 is switched on.
    - test-forest-defaults.R: the held table's row for three or more columns takes the basis's text,
      and `cbind(1 - z, z)` through the data door, held second, turns from created to refused.
2.8 Help and records. man/forest.Rd: the `basis` item and the grammar's table say what a basis is in
    this version, one numeric column or two levels, give the two refusals with the form each shows, and
    lose the rows that describe several columns, a constant beside a column, splines and polynomials of
    a higher degree as accepted; `sd` and `amplitude` lose what they say of several columns. man/bart.Rd,
    "Formula Terms": the line for `basis = a + b` becomes the two-forest form. man/dbartsData.Rd,
    `bases`; [`dbartsSampler$setForestBasis`](../../man/dbartsSampler-class.Rd) with its docstring;
    man/bartBT.Rd's `bases`: one column or two levels. docs/design/forest-kind-by-class.md: the two
    refusals, where the check sits, what is read and named and then refused, and what was kept for
    later with the addition each protects. docs/design/written-surface.md and public-surface.md where
    they say a basis of several columns is accepted. TODO: `factor-basis-refusal` and
    `several-column-basis-refusal` close; `factor-basis-per-level` and `several-column-basis-meaning`
    stay.
2.9 Mutations (Verification).

## Verification

Every push, against the slice's own library (`R CMD INSTALL -l <lib> .`, `R_LIBS=<lib>` on every call;
check `dbarts:::buildInfo()$mode` and that the install postdates the source), run in series, at most two
cores (`MAKEFLAGS=-j2`, `EQUIVALENCE_CORES=2`):

- `cd tests/cpp && make && ./test_bartcore`: unchanged and passing (nothing under src/ moves).
- The full tinytest suite on the shipped build, in one process, counted per file: no failure, no file
  stopping. Push 1: at least the base build's count plus the new file's assertions. Push 2: the landing
  note gives, file by file, the assertions removed with a refused shape and those written in their
  place; a file whose count falls is named with what it no longer covers.
- On a reference build (`--preclean --configure-args=--enable-reference-build`): the four
  `test-reproducibility-*.R` files pass unchanged, and the three compares are bitwise, every scenario
  reporting identical draws, counted per scenario with no `max |z|` line: 55 against the gaussian
  baseline, 15 against `bcf-equivalence-1b7d730c.rds` (its respelled swap among them), 11 against
  `multinomial-equivalence-80b1c8d4.rds`. In push 1 the gaussian baseline is
  `equivalence-1b7d730c.rds` and nothing is re-recorded. In push 2 the first commit is compared with
  it first, 54 identical and `bart2twoforest` not, and then with the merged file, 55; the second commit
  with the merged file, 55.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick`, unchanged: a data object and a fit
  carry another record, and push 1's class is posterior-changing.
- The pair script (below), old side on the base build, new side on the slice's.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD build` with every vignette rebuilt and `R CMD check --as-cran` on
  a tarball from a clean copy.
- The consumers, each suite whole against a private install, none failing. Push 1: bartCause on
  dbarts-1.0 as it stands, and again with step B.1's edit (1412 expectations at its last run). Push 2:
  bartCause with the edit; and once as it stood before it, which must fail at `bcf()` with the second
  [2] text and nowhere else, the proof that the edit is what the push needs. stan4bart on bartcore
  (582), treatSens on dbarts-1.0 (306) and bairrtt on main (207) at both pushes, unedited.

The pair script. Each row is one model fitted on the base build as written there and on the slice's as
written here, the same seed: sampler fits compare the train draws, sigma and the coefficients, `bart`
fits `yhat.train`, sigma and `predict` at 40 new rows. `identical()` on every row but those marked.
Rows 01 to 30 are push 1's against the tip before it; rows 31 to 36 push 2's against push 1.

| | base build | slice |
|---|---|---|
| 01 to 06 | a formula term on a factor, an ordered factor, a character, a logical, a 0/1 number, dose | the same text |
| 07 to 10 | a list on a formula: a factor value, a logical value, `cbind(1 - z, z)`, two numeric columns | the same text |
| 11 to 13 | the matrix interface; `dbartsSpec()`; indicator columns through the data door | the same text |
| 14 | `dbartsData(bases = list(NULL, cbind(1 - z, z)))` | `bases = list(NULL, factor(z))` |
| 15 | a logical vector in a `forests` list | the same vector through the data door: the first changed sequence against its oracle |
| 16 | a logical vector through the data door | the same text: NOT identical |
| 17 to 20 | probit; logistic; three forests of the three kinds; `subset`, weights and an offset | the same text |
| 21, 22 | a swap to `cbind(1 - z2, z2)` on a factor forest, 25 sweeps either side; a swap to `factor(z2)` | `factor(z2)`; the same text |
| 23 | a swap to `factor(v, levels = <those at creation>)` | the same vector with its levels in another order: the second sequence against its oracle |
| 24 | a swap to that vector in the other order | the same text: NOT identical |
| 25, 26 | `predict(bases = )` with indicator columns; with the aligned factor | with the factor; with its levels reversed |
| 27, 28 | a held factor through a list and its block through the data door; `bart` with a term and predict from the stored record | the factor at both; the same text |
| 29 | a copy after a swap, 10 sweeps | the same text |
| 30 | bartCause's `bcf()` as it stands, 11 settings | with step B.1's line |
| 31 to 33 | the three kinds that stay, at a formula term, a list and the data door | the same text |
| 34 | a model of three forests, `arm == "b"` and `arm == "c"` | the same text |
| 35 | `cbind(1 - z, z)` through the data door, drawn and held second | `factor(z)` there |
| 36 | a factor of three levels; `dose + age` | the same text: refused |

Mutations, each expected to fail the named test and no gate before it:

- push 1, a block of 0s and 1s with one 1 a row read as levels at the data door: step 1.1's pins beside
  the loop;
- push 1, a logical vector coerced to a number at the data door; at a forest door: step 1.1's loop;
- push 1, the data door keeping a rule of its own for a factor (another order of levels): step 1.1's
  `identical()` blocks;
- push 1, the levels dropped where rows are restricted on the matrix interface; at the formula door:
  step 1.2's `subset` rows;
- push 1, `dbartsSpec()` leaving the levels the data object had: step 1.2's row for it;
- push 1, the levels left on the stored block as an attribute: step 1.2's `attributes()`;
- push 1, a code basis's levels not recorded (its forest read as numeric): step 1.2's slot at the
  formula door, and step 1.4's table for a forest written as code;
- push 1, the empty check made before `subset`; left out where `dbartsSpec()` installs: step 1.3;
- push 1, a swap matching levels by position: step 1.4's reversed levels and lost level;
- push 1, a level the forest lacks appended as a new column: step 1.4's refusal;
- push 1, `predict` left on the base build's expansion: step 1.4's `predict` rows;
- push 1, indicator columns accepted on a factor forest: step 1.4's table;
- push 1, the data object written before the table is consulted: step 1.4's twin after a refusal;
- push 1, the forest with no basis taken to be forest 1: step 1.4's last test;
- push 1, `copy()` handing on a data object without the levels: step 1.5;
- push 2, a fixture's factor given other levels than the block's column order (`factor(z, levels = 1:0)`):
  the edited file's own pinned draws and pins of `data@bases`;
- push 2, the size check left out at one door, four mutations (the data door; where `dbartsSpec()`
  installs; a value in a list; a basis written as code): step 2.6's loop at that door;
- push 2, the size judged before `subset`: step 2.6's `arm` under a `subset` that leaves two levels;
  judged before the empty check: the same factor as a value, which must meet the droplevels text;
- push 2, two columns let through when they are 0s and 1s with one 1 a row: step 2.6's loop for
  `cbind(1 - z, z)`, and bartCause as it stood;
- push 2, a held basis of three levels let through to the hold's own text, or to the engine: step 2.6's
  loop with `amplitude = fixed()`;
- push 2, a one-column matrix refused as several columns: step 2.6's "what stays";
- push 2, the swap's width left free; the prediction's: step 2.6's "after creation";
- push 2, the validity's check left out: step 2.6's "edited by hand";
- push 2, the message's example written with another forest's predictors or with a level that is the
  first: step 2.6's run of the example.

Not a hot-path change: nothing a sweep runs is touched.

## NEWS

No new item: forests, `dbartsData(bases = )` and `$setForestBasis` are new in 1.0-0 and nothing released
changes.

## What this leaves for the multiplier slice

- The law itself. After this slice a factor and a number differ in what they are called, in what a swap
  and a prediction take, and in nothing the sampler draws. The law gives a number its own coefficient
  variance and default; it reads the kind from the one function this slice adds, and the bridge will
  need it at every construction, which this slice does not hand over.
- The creations on a 0/1 number and the 54 on a constant column of test-bcf-family.R: left as numbers,
  and the law moves their priors. Which of them pin a literal draw is the law's plan to count; a
  constant column is refused by it where no sd is stated.
- The held shapes: forest-defaults-by-kind's two interim refusals by position and dec-A171's of one
  numeric column stand. Which shapes may be held, and at what value, is the law's.
- A sampler saved before this slice has no record and reads as numbers wherever it has a basis. No
  release made one, and the law writes no refusal for one.
- The reader's entry for the kind, the printed line ("factor basis, 2 levels") and the help's advice
  that a treatment belongs in a factor (dec-B263): with the law, where they become true of the prior.
- What the law may assume: every forest has one size, one numeric column or two levels, and keeps it
  for life; every block in the tests, the benchmarks and bartCause that stands for a factor is a factor;
  no accepted swap or prediction changes a forest's kind or a factor's columns; nothing a 1.0 user can
  write has more than one size for a forest.

## What is kept for later, and what is not built

- Kept: push 3's reading and naming of a basis of several terms, to the point of the refusal. It
  protects the later meaning of several numeric columns (TODO `several-column-basis-meaning`) and a size
  for each column (dec-B272): the names a 1.0 message shows are the names those forms will use.
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

- On push 3: it builds a basis written as code through R's model frame and drops its empty levels, so
  step 1.1's "one place" is to be read against its builder, which must still hand a factor term to
  [`expandForestBasis`](../../R/model.R) and nothing else that decides a kind; step 1.2 records the
  levels of a code basis from what push 3 builds; step 2.4 is called where push 3 builds a code basis
  on the rows kept, and step 2.6's last test calls its builder; the texts name a forest by push 3's
  label; every call here is in its spelling, a value held in a variable. dec-B284, ruled on push 3's
  review, moves when a formula object's basis is read, to the call of `forest()`; it changes neither the
  class nor the size of what is read, and `$setForestBasis` is outside it.
- On forest-defaults-by-kind: a forest with no basis may stand anywhere in a list and in a data
  object's `bases`, so "none" is never "forest 1" (step 1.4's last test); `<f>` may be selected by label;
  its interim refusals of held shapes, three of whose texts push 2 replaces or retires and two of which
  wait for the multiplier slice.
- On the sd unit: nothing but the order of edits to shared files, and its reworded text for a longer
  `sd`.
- On bartCause: step B.1, between the pushes. Push 2 is not pushed before it.
- To recheck once they have landed: Context's two tables, by running the two probes again; the counts
  from the traced suite (1186 models, 243 forests of a refused size, the swaps and predictions), by
  running the trace again; that no character or logical basis has entered a test whose block push 2
  would then have to read; that no baseline scenario but `bart2twoforest`, no exact gate and no help
  example fits a refused size (searched at the tip and on push 3's build); the names of the functions
  cited here.
- Order: the sd unit, this slice, the multiplier slice. This slice and the law are not one: the law is
  posterior-changing for every numeric multiplier, and the respelling of push 2 can be proved bit for
  bit only on a tip where the kind is recorded and the law has not moved. This slice and the sd unit
  are not one either: that one changes engine arithmetic.

## Out of scope, and where it goes

- The multiplier slice, with everything listed above; `updateBasisScale`, after the merge to main
  (dec-B276).
- A forest for each level of a factor of three or more (TODO `factor-basis-per-level`), and what several
  numeric columns mean (TODO `several-column-basis-meaning`): after the release.
- A data frame as a basis, and a factor beside numbers in one basis: refused, as today; additions.
- `setModel` given a model whose forest multiplies another kind: there is no model record of a kind to
  differ in until the control-migration arc moves the forests to the model; that arc refuses it.
- To TODO as new entries: a level that is to appear only later in a run cannot be declared, a value's
  empty level being refused at creation; `predict(bases = )` checks no name on a numeric column.

## Calls made in planning

- A 0/1 number is not touched. The brief for this plan had it read as two levels today and become one
  numeric column; measured, it is one numeric column at all seven doors and has been since the basis was
  built, as dec-B260's own text says ("a 0/1 number another"). No tip between this slice and the law
  fits a model nobody asked for, because this slice changes the fit of no number.
- Levels are matched by name at a swap and at predict, and an unknown level is refused (the item under
  "Assumed"). The design accepted any number of levels by position; its critic measured a held contrast
  moving to the other arm. `predict.lm` and this package's own `predict` for a factor predictor refuse a
  level the fit never saw.
- `predict(bases = )` takes the swap's table, numbers on a factor forest refused. The design aligned a
  factor there and said nothing of numbers. Two tables would be a kind judged one way at
  `$setForestBasis` and another at `predict`; the cost is 14 calls in one test file. Lost: a prediction
  at a fractional mix of levels, which no test or consumer makes.
- The kind is stored once, on the data, as the levels; the forests' configuration gains nothing. The
  design put the levels on the data and a `kinds` entry in the configuration. With levels matched by
  name neither a forest's kind nor its levels can change after creation, so a second copy could only
  disagree with the first. The multiplier slice's bridge reads the data object already.
- The levels are a slot and not an attribute of the block, as the design has it: the prototype's
  attribute broke 7 pins of `data@bases` and none of them is this slice's to change.
- The empty check sits where the data object has its rows, and again where `dbartsSpec()` installs a
  declaration. The prototype had it at the first place only and `dbartsSpec()` then created a forest
  with a level no row has: one door reading otherwise, found by the probe of step 1.1's loop. A sampler
  made again from a data object whose swap emptied a column is not a creation and is not checked.
- A character matrix of one column is refused, where the tip reads it as levels at the forest doors:
  "a character vector" is the rule, and a matrix of two columns was a length error.
- A zero column is accepted on a swap of a numeric forest, as an empty level is of a factor: a sampler
  that redraws an indicator inside a larger sampler must not stop on a sweep that empties it.
- Both refusals of dec-B281 and dec-B282 land here, together, in the second push, at every door at
  once. The forest doors know a basis's class today and could have refused sooner; the data door cannot
  tell a factor's indicator columns from two numbers until push 1, and bartCause hands exactly that
  block through it. A refusal at some doors three slices before the others would be one object read
  two ways by where it came in, the fault this slice exists to remove.
- bartCause moves between the pushes and not "the same day": the day a refusal lands is too late for a
  consumer whose only spelling it refuses, and push 1's tip takes both spellings.
- Push 2 is two commits and one gate battery. The respelling must be proved on a tree that still
  accepts both spellings, so it cannot be one commit with the refusal; a third push would buy a second
  battery for a commit that is tests and one help line.
- The size is refused in R and not in the engine. The engine keeps a coefficient for each column of any
  block; narrowing it would have to be undone for either later meaning, and a host other than R may
  mean something else by two columns.
- The size is judged on the rows the fit keeps, after an emptied level has been dropped (code) or
  refused (a value): a factor that `subset` leaves with two levels is the two-level model the caller
  would get by subsetting first.
- A basis of several terms is still read and named before it is refused, and the names are tested off
  the builder. dec-B272 fixes the names for whatever several columns come to mean; a refusal made
  before the columns exist could not show them and would leave the naming rule untested until then.
- The messages write the example in the caller's own words where they are code, and a fixed example
  otherwise. For one term of several columns, a spline, there is no form to show that rebuilds its
  knots at new rows from one column of it, so the text shows the fixed example and does not suggest a
  subscript.
- `1 + dose` has its own sentence: its two columns are not two things the caller listed, and "a forest
  for each column" would tell them to multiply a forest by a constant.
- A block of three indicator columns in a test is not respelled as a factor of three levels in the
  first commit to be refused in the second; its test goes in the second.
- The baseline's scenario is given a factor of two levels, not the two-forest form the message shows:
  it is there as the one model of several forests the formula door reaches in that file, K = 2.
- The pair script is run at landing and not tracked: its old side needs the base build. What stays in
  the suite is step 1.1's loop across doors, step 1.4's identities and step 2.6's loop.
- Measured for this plan on a prototype of push 1 (R only, 141 changed lines, the levels as an
  attribute, no data slot, no check at `dbartsSpec()`, forests named by position): the 252 cells, the
  suite's 6 stopping files and 8 failures, the 37 refusals passed by, and bartCause's three-way
  identity. Read off push 3's build, with nothing of this slice in it: what each door does with the
  refused sizes, and the traced counts. Not built: the slot, the validity, the labels in texts, the
  size check, push 3's code path.
