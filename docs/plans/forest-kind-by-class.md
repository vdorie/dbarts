# forest-kind-by-class: what a forest multiplies is decided by the class of its basis, at every door, and fixed at creation

Status: PLANNED (dec-B260; dec-B254 for the kind as structure). Follows the sd-unit slice
(forest-sd-unit.md), [forest-defaults-by-kind.md](forest-defaults-by-kind.md) and push 3 of
[written-surface.md](written-surface.md), none of which has landed.

agent: two pushes. Push 1 (the rule): opus implementer for the R code, sonnet for the help once the code is
fixed, opus reviewer told to refute. Push 2 (the blocks written by hand): sonnet implementer, opus
reviewer. The reason for opus on push 1: it moves no arithmetic, and every slip in it is a forest read as
one kind at one door and another at the next, or a coefficient multiplying another level's rows after a
swap, with no message. The tip has both faults today (Context), and the three reviews of written-surface's
second push each found a fault of that shape one step to the side of the last.
rng: stated per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- NEUTRAL for every fit that both builds accept, the three sequences below apart. Measured on a prototype
  of the rule: 36 kinds of object at 7 creation doors, 252 cells, of which 227 have the same width, columns
  and seeded draws on both builds and 3 are the first sequence below; the suite, with the new refusals
  passed by, 16065 results and no seeded draw moved; 11 seeded `bcf()` fits of bartCause.
- POSTERIOR-CHANGING for three sequences, none in a test, a benchmark, a help example or a consumer
  (counted by running the suite on a build that logs every basis, and searched elsewhere):
  1. a logical vector handed to `dbartsData(bases = )`: one 0/1 numeric column before, two levels after,
     which is what the same vector is at every other door. The oracle is the base build's fit of it
     through a `forests` list: identical draws.
  2. `$setForestBasis` on a forest created with a factor or a character vector, given one whose levels
     stand in another order or lack one of them: matched by position before, by name after. The oracle is
     the base build given the same vector as `factor(v, levels = <the levels at creation>)`.
  3. `predict(bases = )` with such a vector: the same, with the same oracle.
- Refused where accepted, and accepted where refused: listed under "The rule". No draw is involved.
- Push 2 moves no draw: each respelled call is the same block of 0s and 1s.
window: pre-release, after the sd unit and before the multiplier law, which makes a factor and a block of
numbers two different priors: until it lands a factor and its indicator columns handed over as numbers
are one fit, so this is the only tip on which every such block in the tests, the benchmarks and bartCause
can be respelled and the respelling proved bit for bit. Serial with any other work in
[`expandForestBasis`](../../R/model.R), [`validateForestBases`](../../R/data.R),
[`setForestBasis`](../../R/dbarts.R), [`resolveForestBases`](../../R/generics.R) or the multi-forest block
of [`resolveSamplerSpec`](../../R/spec.R). No engine, bridge or header file changes.
budget: ~1650 lines changed (R ~300; tinytest ~900, of which ~550 in one new file and ~350 changed in
about sixteen; benchmarks ~5; help ~150; design note, architecture, public-surface, TODO and the two
indexes ~200; the rest slack for texts), upper figure 3000 (R 600, tinytest 1700, help 300, records 350,
benchmarks 50). By push: 1 ~1150 (upper 2100), 2 ~500 (upper 900). The design estimated 900 and planned
for 1700; written-surface's two landed pushes ran at 1.8 and 2.7 times their plans.

## Assumed, pending the maintainer

One point of what a user sees is settled here by the design's critic and not by a ruling. It is planned
as written and flagged in the notes for the coordinator.

- A factor handed to `$setForestBasis` or to `predict(bases = )` is matched to the forest's levels by
  name, and a level the forest was not created with is refused. The tip matches by position (Context
  shows what that does); the design had "any number of levels" accepted. If the ruling is "by position,
  documented": step 1.4's alignment is dropped, its texts about a level go, and the second and third
  changed sequences go with them.

## Goal

A forest multiplies nothing, the levels of a factor, or numeric columns, and which of the three follows
from the class of what the caller handed over and from nothing else: a factor, a character vector or a
logical vector is levels; a number is numbers, whatever its values. One object is one kind at every door,
`dbartsData(bases = )` included. The kind and a factor's levels are recorded when the sampler is created
and hold for its life: `$setForestBasis` and `predict(bases = )` take the class the forest was created
with, a factor is matched to the recorded levels by name, and a forest created without a basis takes
none.

## Context

Measured at the tip (4f2e79e6; written-surface pushes 1 and 2 landed, push 3 and forest-defaults-by-kind
not) on the shipped build, R 4.6.1: 150 rows; x1, x2; z a 0/1 double with 66 zeros and 84 ones, and the
same column as an integer, a logical, a factor, an ordered factor and a character vector; g a character
vector of three levels; dose, age. The spellings are the tip's. "Draws" is the sum of 10 seeded train
draws of a two-forest sampler, 15 trees; the same number is the same fit.

- The suite. 231 files, of which 4 exit off a reference build and 3 ask for more than two threads; the
  other 224 give 16065 results, none failing, 102 seconds in one process.
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

  Push 3 of written-surface drops an emptied level for a basis written as code, at both doors. A value
  is this slice's.
- After creation. `$setForestBasis` takes any class on any forest: 8 kinds of forest by 16 kinds of value
  were tried, and every pair is accepted except a factor of one level, a column of zeros, a logical
  matrix, and one numeric column on a held forest (dec-A171). A forest created with no basis takes one,
  of any kind; a factor forest takes dose, or `cbind(dose, age)`; a numeric forest takes a factor; the
  width follows the value.
- A factor's levels are matched by position. A held two-level character forest (control, treated):
  mean absolute term 0.000 on the control rows and 0.129 on the treated; after a swap to the same vector
  as a factor with levels (treated, control), 0.133 and 0.000. A drawn three-level forest with
  coefficients -0.785, -0.158, 1.436: swapped to the vector with one level redrawn into another it has 2
  columns and the rows of the third level take the second coefficient; swapped back, the third
  coefficient is 1. A vector left with one value is refused.
- `predict(bases = )` expands what it is given by its own levels and checks the width. On a fit created
  with `factor(z)`: the same factor with its levels reversed is accepted and predicts 1.267 where the
  aligned one predicts 1.433; a two-level character vector unrelated to z is accepted; `cbind(dose, age)`
  is accepted (mean prediction 26.7). On a three-level fit a vector with one level renamed is accepted
  and takes the old level's coefficient. On a fit created with `cbind(dose, age)` a factor, a logical
  and `cbind(age, dose)` are accepted. On a fit created with a 0/1 number a factor of one level is
  accepted as a column of ones.
- What is recorded. A data object has `bases`, bare numeric matrices, and no record of levels; the
  forests' configuration has no kind. A fit carries `bases` and, for a formula term, the term with its
  levels ([`packageBartResults`](../../R/bart.R), [`replayForestBasis`](../../R/model.R)).
- What the suite hands over, counted by running it on a build that logs every basis. 776 models of
  several forests are created in 53 files. Their bases: a factor 281 times; a 0/1 integer or double
  vector 403 times in 15 files (239 in test-forest-predictors.R, 127 in test-formula-terms.R); another
  numeric vector 28; a numeric matrix at a forest door 32, of which 7 are a factor's indicator columns;
  through the data door 139, of which 91 are indicator columns, 13 indicator columns times a constant,
  33 a constant column and 2 other. Of the 98 blocks of indicator columns, 61 are in test-bcf-family.R
  and 14 in test-bcf-creation.R. No character vector, logical vector or ordered factor reaches a
  creation.
- `$setForestBasis` is named on 53 lines of 11 test files and reaches a sampler 29 times in 7. Under the
  rule 20 of those meet a refusal: 13 put numbers on a factor forest (test-forest-basis-r5.R 4,
  test-forest-arguments.R 4, test-bcf.R 2, and one each in test-bcf-mutation-pins.R,
  test-multi-forest-seam.R and test-multiforest-leaf-prior-writer.R), two of which are today's tests of
  dec-A171's refusal, and 7 give a basis to the forest with none (test-forest-basis-r5.R 4, and one each
  in test-bcf-r5-surface.R, test-forest-arguments.R and test-multiforest-leaf-prior-writer.R). Three
  more that are refused today for their length would meet the rule first. `predict` is given `bases` 64
  times in 6 files, always numbers; 14 of them, all in test-predict-blend.R, hand indicator columns to a
  factor forest.
- The prototype against the suite. As the rule stands, 6 files stop at their first refused call
  (test-bcf.R, test-forest-arguments.R, test-forest-basis-r5.R, test-multi-forest-seam.R,
  test-multiforest-leaf-prior-writer.R, test-predict-blend.R: 2895 lines, 545 results at the tip). With
  the refusals logged and passed by, 16065 results run and 8 fail: one pin of the droplevels text, and 7
  pins of `data@bases` against a bare matrix that fail only because the prototype keeps the levels as an
  attribute of the block, which the plan does not.
- Elsewhere. benchmarks/R: every basis at creation is a factor, but for the `bcf` row of
  `composition-matrix.R`, a 0/1 number; `bcf-equivalence.R` swaps `cbind(1 - z2, z2)` onto a factor
  forest once, in the scenario `set_treatment`. Help: the example of man/bart.Rd multiplies by a 0/1
  integer and calls it the Bayesian causal forest shape. The vignettes fit no model of several forests.
- Consumers. bartCause's `bcf()` (dbarts-1.0, 6c1bff9) hands `cbind(1 - z, z)` to `dbartsData(bases = )`
  under its own `subset`, having refused a treatment that the kept rows leave with one value; its
  treatment forest takes `amplitude = fixed()` when `update.b = FALSE`, which on two columns holds (0,
  1) and is not the one-column case dec-A171 refuses. Two of its test files build the same sampler by
  hand from the same block. Run for this plan: 11 seeded `bcf()` fits (the default, a binary and a
  logistic response, `subset`, a logical treatment, weights with an offset, stated sds, each coefficient
  held and both, another coefficient variance) are identical between bartCause as it stands on the tip,
  bartCause with `basis <- factor(as.integer(z), levels = 0:1)` on the prototype, and bartCause as it
  stands on the prototype; its two bcf test files pass on all three, 156 and 47 expectations. stan4bart
  (bartcore, a9d081b), treatSens (dbarts-1.0, aecec71) and bairrtt (main, 3f57f61) declare one forest
  and no basis.

## The rule

At creation, at every door:

| handed over | kind | coefficients |
|---|---|---|
| nothing | none | one, on the forest itself |
| a factor, ordered or not; a character vector; a logical vector; `I()` of any of them | levels | one per level, in the factor's own order of levels, the sorted values of a character vector, `FALSE` then `TRUE` |
| an integer or double vector or matrix, whatever its values; `I()` of one | numbers | one per column |
| a logical or character matrix; a data frame; anything else | refused by name | |

- The class is looked at in one place, [`expandForestBasis`](../../R/model.R), which every creation
  door reaches. The data door hands it what it is given.
- The content is never looked at to decide a kind: a 0/1 number, a factor's indicator columns and
  dummies with a level left out are numbers with every row, with some rows and with none at one value.
- A basis written as code is the class its code gives in the fit's model frame, as push 3 builds it;
  the same vector handed over as a value is that class too.
- Empty at creation. A level no kept row has, and a numeric column that is zero on every kept row, are
  refused for a value, once, where the data object has its rows: after `subset` and the na.action, at
  every door. Code has dropped its empty levels by then (push 3).
- The record. A factor's levels are kept with the data, per forest; a forest's kind is read from the
  data by one function and from nowhere else.

After creation, at `$setForestBasis` and at `predict(bases = )`, by one function:

| forest created as | takes | refuses |
|---|---|---|
| none | nothing | any basis |
| levels | a factor, a character vector or a logical vector, matched to the forest's levels by name; a level with no row is a column of zeros | numbers, indicator columns included; a level the forest does not have |
| numbers | a numeric vector or matrix; a column of zeros | a factor, a character vector, a logical vector |
| either, coefficient held | the same, at the width it was created with | another width |

- A levels forest has the same columns for life: its width cannot change, and coefficient j multiplies
  level j's rows after any swap.
- A numeric forest whose coefficient is drawn takes another width, as today, the added coefficients
  entering at 1; the multiplier law narrows that for a forest at its default sd.
- At `predict` a numeric block has the fit's width, as today.
- A basis written as code and rebuilt at new rows from the fit's own record is untouched.
- Nothing of this moves a draw of a sampler that is swapped to a block the tip would have built the same
  way.

Accepted where refused: a factor, an ordered factor and a character vector at the data door; a factor of
one value at a swap and at predict, where its level is one of the forest's (every new row treated); a
column of zeros on a swap of a numeric forest. Refused where accepted: a logical matrix, a column of
zeros and an empty level at the data door; an emptied level or column of a value under `subset` or a
missing response, at the list and the matrix doors; a character matrix of one column; every cell of the
table's last column.

Before and after, per door.

| door | before | after |
|---|---|---|
| a `forest()` term at `bart()` and `dbarts()`; a `forests` list on a formula, on a matrix, in `dbartsSpec()` | by class | by class; the same fits |
| `dbartsData(bases = )` | numbers, a logical vector among them; a factor refused | by class |
| `$setForestBasis` | any class, any forest, levels by position | the table |
| `predict(bases = )` | any class at the fit's width, levels by position | the table |
| `predict` from a formula term's stored record; `fitted`; `extract` | | unchanged |
| `$getLeafPrior`, `$setLeafPrior`, `print`, `show` | say nothing of a kind | unchanged; the reader's entry and the printed block are the multiplier law's |
| `copy()`, a reload, `new("dbartsSampler", control, model, data)` | carry the block | carry the block and its levels: the same kind |
| a state stored or installed | carries no basis | unchanged |

## Refused forms, with their texts

Base R's style, as in written-surface. `<arg>` is `basis` at a forest door and at `$setForestBasis`,
`bases` at the data door and at `predict`. `<f>` is `forest 2`, or `forest 2 ("dose")` where push 3 gives
the forest a label. The push that adds each is in brackets.

    [1] '<arg>' must be a factor, a character or logical vector, or numeric; a logical matrix is none of these: write it as numbers (<arg> + 0), or give one logical column            (a character matrix: "give one character column")
    [1] '<arg>' must be a factor, a character or logical vector, or numeric, not a data frame; give one of its columns, or its numeric columns as a matrix            (a Date, a list, a complex: named by its class)
    [1] <f> has a basis level with no observations: "c"; drop it with droplevels()
    [1] <f> has a basis column of all zeros (column 2), which contributes nothing; drop it
    [1] <f> multiplies nothing: it was created without a basis, and what a forest multiplies is fixed when the sampler is created; make a new sampler to give it one
    [1] <f> multiplies a factor, so '<arg>' must be a factor, a character vector or a logical vector; numbers are another model, an indicator matrix included: write factor(), or make a new sampler to multiply numbers
    [1] <f> multiplies numbers, so '<arg>' must be a numeric vector or matrix, not a logical vector; write as.numeric(<arg>) for a 0/1 column, or make a new sampler to multiply a factor            (a factor, a character vector)
    [1] <f> has no level "s"; its levels are "p", "q", "r", fixed when the sampler was created
    [1] $setForestBasis cannot change the width of the basis of <f> (2 to 3): its coefficient is held (amplitude = fixed()), and the held value is defined for that width only; make a new sampler
    [1] 'basis.levels' gives forest 2 3 levels and its basis has 2 columns

The last is the data object's validity, for an object edited by hand. Replaced: "'bases' must be numeric
or logical", "a 'basis' must be numeric, a factor, or a character vector", "a 'basis' factor level with
no observations contributes nothing to a forest" and "a 'basis' column of all zeros contributes nothing
to a forest". Kept as they are: "a 'basis' cannot be NA", "a 'basis' factor must have at least two
levels" at creation, the texts about a basis's length, the at-least-two-forests text, dec-A171's refusal
of a held single numeric column at creation, and at `predict` the texts for a forest given a basis it has
none of and for a numeric block of another width. A text of push 3 for a basis written as code (a logical matrix
as a term, a logical with one value) comes first where both apply.

## Constraints

- Every fit that both builds accept, the three sequences of the `rng:` line apart, draws what it draws
  on the base build. No baseline is re-recorded and no snapshot file regenerated.
- No change to src/, to the flat C header, to the stored state or to what the bridge reads.
- The class is read in one function at creation and the table applied by one function afterwards; no
  door has a copy of either.
- `data@bases[[f]]` stays a bare numeric matrix, with the names push 3 gives it and no other attribute.
- No law changes: a block of indicator columns handed over as numbers and the factor it came from are
  one fit after this slice as before it.
- The reader, `extract` and `print` are untouched. The unit and the default of an sd are untouched.
- Nothing of `updateBasisScale` rides along, and no width rule for a forest at its default sd.
- Each push leaves the help saying what the code does. Base R calls stay within DESCRIPTION's R floor.

## Pushes

1. The rule. Changes the three sequences; the design note lands here.
2. The blocks written by hand: every call in the tests, the benchmarks and the help that hands a
   factor's indicator columns over as numbers, or swaps them onto a factor forest, is read and either
   written as the factor it stands for or left as the numeric forest it is. No draw moves. bartCause's
   line lands any time from push 1 until the multiplier law.

Push 1 alone leaves a coherent tip: its own respelling is the 20 swaps and 14 predictions it refuses,
which it carries. Push 2 is what the multiplier law needs done before it.

## Steps

"Fails today" is what the base build does where the test expects otherwise. New names are for the
reader's sake. Calls are in push 3's spelling; the fixture is Context's.

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
    logical matrix and `rep(0, n)`). Pinned beside it, so that a reading by content fails: a 0/1
    integer, a 0/1 double, `I(z)` and `matrix(z)` are numbers of one column and do NOT draw what
    `factor(z)` draws; `cbind(1 - z, z)` and `model.matrix(~ g - 1)` are numbers and DO, with a comment
    that the multiplier law turns that assertion around; `model.matrix(~ g)[, -1]` under a `subset`
    that leaves no row at the omitted level is numbers. The same column written as code and handed over
    as a value is one kind and one fit.
1.2 The record. A slot `basis.levels` on [`dbartsData`](../../R/A_class.R), beside `response.levels`:
    `NULL`, or a list with one entry per forest, `NULL` for a forest with no basis or with numbers and
    the levels for one with a factor. It is written wherever `bases` is: where
    [`dbartsData`](../../R/data.R) finishes its rows, with the levels carried across every restriction
    of rows ([`restrictBasesToRows`](../../R/data.R), [`alignForestBasisToSubset`](../../R/model.R)) and
    taken off the block before it is stored; and where `dbartsSpec()` puts a declaration's bases on a
    data object ([`resolveSamplerSpec`](../../R/spec.R)), replacing what the object had. It is read
    through an accessor that answers `NULL` for an object saved without the slot, as
    [`dataRowNames`](../../R/data.R) does. One function, `forestBasisKind(data, f)`, answers "none",
    "levels" or "numeric" and is the only reader of a kind after creation. The class's validity refuses
    recorded levels whose count is not the block's width.
    [`packageBartResults`](../../R/bart.R) puts the levels on a fit beside `bases`. Tests: the slot at
    each of the seven doors for a factor, a character and a logical vector, a number and no basis; the
    same under `subset` at the formula, list, matrix and data doors, and after a dropped missing
    response (a record lost with the rows makes the forest numeric there); `dbartsSpec()` over a data
    object that carried a numeric block for that forest records the declaration's levels;
    `attributes(data@bases[[2]])` holds the dimensions and names and nothing else; an object with the
    slot removed creates, runs and reads as numbers; an object with levels edited to another count is
    refused with the last text.
1.3 Empty at creation. One check where the data object has its final rows, and the same check where
    `dbartsSpec()` installs a declaration: a level with no row, or a column of all zeros, is refused
    with the third and fourth texts. Tests: a factor value with a level never used; one whose level
    `subset` empties; one emptied by a missing response; a number left all zero by `subset`;
    `rep(0, n)`; indicator columns with one all zero: each at the list, matrix, data and `dbartsSpec`
    doors, one text (fails today: all but the first and fifth are created at the list and matrix doors,
    and the fifth at the data door); and beside the value the same column written as code under the
    same `subset`, which push 3 creates with the level dropped, pinned so that the two readings stay
    told apart.
1.4 After creation. One function (`conformBasis`: a value, the forest's kind, its levels, its width,
    whether its coefficient is held) returns the block for the engine or refuses by "The rule";
    [`setForestBasis`](../../R/dbarts.R) and [`resolveForestBases`](../../R/generics.R) both call it,
    the sampler reading kind and levels from its data object and a fit from what it carries. In
    `$setForestBasis` nothing is stored before it returns. A one-sided formula is evaluated first, as
    today. dec-A171's refusal of one numeric column on a held forest can no longer be reached at a swap
    (a held forest keeps its width and a factor forest takes no number) and goes from there; at
    creation it stays. Tests:
    - The table: 8 kinds of forest by 16 kinds of value at `$setForestBasis`, 4 by 13 at `predict`,
      each cell accepted or refused with its text (fails today: all but four kinds of value are
      accepted everywhere). After every refused swap the data object, `getForestAmplitudes()` and the
      next 5 sweeps are `identical()` to an untouched twin's.
    - Which model, by name. A held two-level character forest swapped to the same vector with its
      levels reversed keeps its term on the treated rows (fails today: it moves to the control rows). A
      three-level forest swapped to a vector that has lost its middle level keeps 3 columns, the middle
      one zero, and its coefficients; swapped back, the three coefficients are `identical()` to those
      before (fail today: 2 columns, and the third restarts at 1). A vector left with one value is
      accepted (fails today: refused). A level the forest lacks is refused.
    - `predict`: a factor with its levels reversed predicts `identical()` to the aligned one (fails
      today: 1.267 against 1.433); every row at one level, written `rep("1", n)`, predicts what
      `factor(rep(1, n), levels = 0:1)` does today; a renamed level is refused (fails today: taken);
      numbers on a factor forest and a factor on a numeric forest are refused (fail today: taken where
      the width fits).
    - A column of zeros swapped onto a numeric forest is accepted and the run goes on (fails today:
      refused). A held forest of two numeric columns refuses three.
    - With the forest with no basis second in a list, as forest-defaults-by-kind accepts it, that
      forest refuses a basis and the first takes one: the kind is not a position.
1.5 The same sampler again. Tests: after `copy()`, after `saveRDS`, `readRDS` and a first use, and on
    `new("dbartsSampler", control, model, data)` from the sampler's own three, the kind and the levels
    are those of the original, the by-name swap of step 1.4 gives the same block, and a refused swap is
    refused; a state taken from a sampler created on the same factor with its levels in another order
    installs, and the recipient's levels are its own.
1.6 Respell what the rule refuses. The 20 swaps and the 14 predictions of "Context" (the two tests of
    dec-A171's refusal at a swap become tests of the held width): a numeric block
    swapped onto, or predicted for, a factor forest becomes the factor (`cbind(1 - z2, z2)` becomes
    `factor(z2)`); a basis given to the forest with none becomes a pin of the fifth text where the test
    is about that, and otherwise the forest is created with a basis. The three swaps refused today for
    their length are respelled so that the length is still what is tested. The pin of
    ["factor level with no observations contributes nothing"](../../inst/tinytest/test-forest-basis-r5.R)
    takes the third text. `bcf-equivalence.R`'s swap becomes `factor(z2)`; its scenario is unchanged to
    the bit. Run the suite first and repair what it shows.
1.7 Help and records. man/forest.Rd, the `basis` item as push 3 leaves it: which class is which kind, in
    the three lines of "The rule"; that a 0/1 number is one column with one coefficient and not two
    levels, and `factor(z)` is how levels are written; that a value's empty level or column is refused
    and code's is dropped. man/dbartsData.Rd, `bases`: the same classes as `forest(basis = )`.
    [`dbartsSampler$setForestBasis`](../../man/dbartsSampler-class.Rd) and the method's docstring: the
    table, and that levels are matched by name. man/bartBT.Rd, the `bases` item of `predict`: the
    table. docs/design/forest-kind-by-class.md with its index row: the rule, the record, the three
    changed sequences with their oracles, what a sampler saved before reads as.
    docs/design/public-surface.md and docs/architecture.md where they describe a basis and the data
    object's slots. TODO: `forest-prior-args` names this slice landed.
1.8 Mutations (Verification): apply each, install, run the named test, record the failing count, revert,
    `touch` the file.

### Push 2: the blocks written by hand

2.1 Tests. Each of the 98 creations from a factor's indicator columns (91 through the data door, 7 at a
    forest door, in 10 files) is read at its fixture. Where the block stands for a factor (the
    treatment and group blocks of test-bcf-family.R, test-bcf-creation.R, test-bcf-loglik.R,
    test-bcf-reporting.R, test-forest-basis-subset.R and the others the log names) the fixture hands
    the factor over, and a pin of `data@bases` against the block is kept, against the block. Left as
    numbers, each with a comment that it is meant so: the 13 blocks times a constant and the 33
    constant columns of test-bcf-family.R, which test the tip's law and which the multiplier law
    rewrites; and the 403 creations on a 0/1 number, which are one numeric column before and after.
    Every edited file's results are those of the base build, assertion for assertion.
2.2 Benchmarks and help. `composition-matrix.R`'s `bcf` row is written `basis = factor(z)`; man/bart.Rd's
    example of the Bayesian causal forest shape is written with `factor(z)`, and the line of its
    section "Formula Terms" that glosses a numeric z says one coefficient, not a treatment. These two
    change the model the text fits, which nothing pins.
2.3 bartCause (its own commit on dbarts-1.0, any day from push 1 to the multiplier law). R/bcf.R:
    `basis <- cbind(1 - z, z)` becomes `basis <- factor(as.integer(z), levels = 0:1)`; `as.integer`
    because `bcf()` takes a logical treatment, which `factor(z, levels = 0:1)` would turn into missing
    values. The comment above it loses its sentence about the column order of `cbind`, the order now
    being that of the levels 0, 1. Its two hand-built samplers (tests/testthat/test-14-bcf.R, three
    blocks; test-03-responseFit.R, one) are written with `factor()` in the same commit, so that they
    stay the comparator of `bcf()` under the multiplier law. No other file there reads a basis.

## Verification

Every push, against the slice's own library (`R CMD INSTALL -l <lib> .`, `R_LIBS=<lib>` on every call;
check `dbarts:::buildInfo()$mode` and that the install postdates the source), run in series, at most two
cores (`MAKEFLAGS=-j2`, `EQUIVALENCE_CORES=2`):

- `cd tests/cpp && make && ./test_bartcore`: unchanged and passing (nothing under src/ moves).
- The full tinytest suite on the shipped build, in one process, counted per file: no failure, no file
  stopping, and at least the base build's count plus the new file's assertions; the landing note gives
  the figures.
- On a reference build (`--preclean --configure-args=--enable-reference-build`): the four
  `test-reproducibility-*.R` files pass unchanged, and the three compares are bitwise, every scenario
  reporting identical draws, counted per scenario with no `max |z|` line: 55 against
  `equivalence-1b7d730c.rds`, 15 against `bcf-equivalence-1b7d730c.rds` (its respelled swap among them),
  11 against `multinomial-equivalence-80b1c8d4.rds`. Nothing is re-recorded.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick`, unchanged: a data object and a fit
  carry another record, and push 1's class is posterior-changing.
- The pair script (below), old side on the base build, new side on the slice's.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD build` with every vignette rebuilt and `R CMD check --as-cran` on
  a tarball from a clean copy.
- The consumers, each suite whole against a private install of the slice, none failing: bartCause on
  dbarts-1.0 as it stands and again with step 2.3's edit (1412 expectations at its last run), stan4bart
  on bartcore (582), treatSens on dbarts-1.0 (306), bairrtt on main (207).

The pair script. Each row is one model fitted on the base build as written there and on the slice's as
written here, the same seed: sampler fits compare the train draws, sigma and the coefficients, `bart`
fits `yhat.train`, sigma and `predict` at 40 new rows. `identical()` on every row but those marked.

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
| 30 | bartCause's `bcf()` as it stands, 11 settings | with step 2.3's line |

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
- push 1, the empty check made before `subset`; left out where `dbartsSpec()` installs: step 1.3;
- push 1, a swap matching levels by position: step 1.4's reversed levels and lost level;
- push 1, a level the forest lacks appended as a new column: step 1.4's refusal;
- push 1, `predict` left on the base build's expansion: step 1.4's `predict` rows;
- push 1, indicator columns accepted on a factor forest: step 1.4's table;
- push 1, the data object written before the table is consulted: step 1.4's twin after a refusal;
- push 1, the forest with no basis taken to be forest 1: step 1.4's last test;
- push 1, the held-width refusal dropped: step 1.4;
- push 1, `copy()` handing on a data object without the levels: step 1.5;
- push 2, a fixture's factor given other levels than the block's column order (`factor(z, levels = 1:0)`):
  the edited file's own pinned draws and pins of `data@bases`.

Not a hot-path change: nothing a sweep runs is touched.

## NEWS

No new item: forests, `dbartsData(bases = )` and `$setForestBasis` are new in 1.0-0 and nothing released
changes.

## What this leaves for the multiplier law

- The law itself. After this slice a factor and a block of numbers differ in what they are called, in
  what a swap and a prediction take, and in nothing the sampler draws. The law gives numbers their own
  coefficient variance and default; it reads the kind from the one function this slice adds, and the
  bridge will need it at every construction, which this slice does not hand over.
- The 403 creations on a 0/1 number in 15 test files, the 13 scaled blocks and the 33 constant columns:
  left as numbers, and the law moves their priors. Which of them pin a literal draw is the law's plan to
  count; the constant columns are refused by it where no sd is stated.
- The width of a numeric forest on a swap: free here where the coefficient is drawn. The law fixes it
  for a forest at its default sd, and its text ends "state one with $setLeafPrior".
- The held shapes: here a held forest keeps its width and dec-A171's interim refusal stands. Which
  shapes may be held, and at what value, is the law's.
- A sampler saved before this slice has no record and reads as numbers wherever it has a basis. No
  release made one. The law decides whether such a sampler is refused or run as numbers.
- The reader's entry for the kind, the printed line ("factor basis, 2 levels") and the help's advice
  that a treatment belongs in a factor (dec-B263): with the law, where they become true of the prior.
- What the law may assume: every block in the tests, the benchmarks and bartCause that stands for a
  factor is a factor; no accepted swap or prediction changes a forest's kind or a factor's columns; a
  levels forest never changes width.

## What waits on what

- On push 3: it builds a basis written as code through R's model frame and drops its empty levels, so
  step 1.1's "one place" is to be read against its builder, which must still hand a factor term to
  [`expandForestBasis`](../../R/model.R) and nothing else that decides a kind; step 1.2 records the
  levels of a code basis from what push 3 builds, and they must be the levels its stored `terms` would
  rebuild at new rows; step 1.4's swap keeps push 3's rule for the names of a numeric block; the texts
  name a forest by push 3's label; every call here is in its spelling.
- On forest-defaults-by-kind: a forest with no basis may stand anywhere in a list and in a data
  object's `bases`, so "none" is never "forest 1" (step 1.4's last test); `<f>` may be selected by label.
- On the sd unit: nothing but the order of edits to shared files.
- To recheck once they have landed: Context's two tables, by running the two probes again (the list
  door under `subset` changes with push 3); the counts from the logging build (776, 403, 98, 29, 64);
  that no character or logical basis has entered a test whose block push 2 would then have to read; the
  names of the functions cited here.
- Order: the sd unit, this slice, the multiplier law. This slice and the law are not one: the law is
  posterior-changing for every numeric multiplier, and the respelling of push 2 can be proved bit for
  bit only on a tip where the kind is recorded and the law has not moved. This slice and the sd unit
  are not one either: that one changes engine arithmetic.

## Out of scope, and where it goes

- The multiplier law, with everything listed above; `updateBasisScale`, after the merge to main
  (dec-B276).
- A data frame as a basis, and a factor beside numbers in one basis: refused, as today; additions.
- `setModel` given a model whose forest multiplies another kind: there is no model record of a kind to
  differ in until the control-migration arc moves the forests to the model; that arc refuses it.
- To TODO as new entries: a level that is to appear only later in a run cannot be declared, a value's
  empty level being refused at creation; `predict(bases = )` checks no names on a numeric block, which
  `cbind(age, dose)` for a forest created on dose and age passes (push 3 fixes the swap, not
  `predict`).

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
  disagree with the first. The multiplier law's bridge reads the data object already.
- The levels are a slot and not an attribute of the block, as the design has it: the prototype's
  attribute broke 7 pins of `data@bases` and none of them is this slice's to change.
- The empty check sits where the data object has its rows, and again where `dbartsSpec()` installs a
  declaration. The prototype had it at the first place only and `dbartsSpec()` then created a forest
  with a level no row has: one door reading otherwise, found by the probe of step 1.1's loop. A sampler
  made again from a data object whose swap emptied a column is not a creation and is not checked.
- A character matrix of one column is refused, where the tip reads it as levels at the forest doors:
  "a character vector" is the rule, and a matrix of two columns was a length error.
- A held forest keeps its width, of either kind, as the design's table has it. For levels it follows
  from the rule; for numbers it is one more refusal, and without it a held block of two swapped to
  three would take a held value nobody defined.
- A zero column is accepted on a swap of a numeric forest, as an empty level is of a factor: a sampler
  that redraws an indicator inside a larger sampler must not stop on a sweep that empties it.
- Two pushes. The rule's own respelling is 34 calls it refuses; the 98 blocks of push 2 are accepted
  either way and are respelled for the law's sake. Cutting there lets the rule be reviewed against
  model risk and the respelling against nothing but "the same draws".
- Left as numbers in push 2: the scaled and constant blocks. They test the row norm and the anchor of
  the tip's law, which the law removes with their tests.
- The pair script is run at landing and not tracked: its old side needs the base build. What stays in
  the suite is step 1.1's loop across doors and step 1.4's identities.
- Measured for this plan on a prototype (R only, 141 changed lines, the levels as an attribute, no data
  slot, no check at `dbartsSpec()`, forests named by position): the 252 cells, the suite's 6 stopping
  files and 8 failures, the 37 refusals passed by, and bartCause's three-way identity. Not built: the
  slot, the validity, the labels in texts, push 3's code path.
