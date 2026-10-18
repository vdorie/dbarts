# factor-levels-across-setdata: setData reads a factor by label against the sampler's levels

Status: LANDED 2026-10-17 (92ae2714..6beef48d; dec-B454, dec-A213; [Landing note](#landing-note); follow-up 3eece920..7318bd9b, dec-A215); option 5 (C
below) ruled by the maintainer (dec-B454); written to the tree of
f514c1d8 (origin/bartcore); revised from its blind critique the same day ([Revised from the blind
critique](#revised-from-the-blind-critique)), which kept the recommendation. Closed TODO `factor-levels-across-setdata`. Rulings it stands on: dec-A115 (declared
levels), dec-A85 and dec-B155 (a level matched by label, an unseen one refused), dec-B279, dec-B286,
dec-B359 and dec-B360 (every route that takes new predictors codes a factor column by label and refuses
numbers for it), dec-B231 (setData changes derived data only unambiguously), dec-B195, dec-B418 and
dec-A210 (a state carries no data), dec-B449 and dec-B450 (a state installs no grid), dec-B230; dec-B346
by analogy only. It settles what dec-A209 scoped out. The maintainer's call
([For the maintainer](#for-the-maintainer)) is ruled: option 5, dec-B454.

agent: tier "Changes draws" ([Process by risk](README.md#process-by-risk)): the diff is R only, but four
call sequences fit another model after it. The blind critique ran (scratch/flsc-crit.md). One opus
implementer for the R code and tests (a slip here reads one level as another with no message, the defect
being fixed); one opus review that runs the mutants under [Tests](#tests) and checks the fix round.
Nothing under src/ changes, so no new numerics are reachable and no sanitizer run is owed.
rng: NEUTRAL, bit for bit, for every `setData` whose columns carry the sampler's own level tables and
names, and for every fit that never calls `setData`. POSTERIOR-CHANGING for four sequences, none in a
benchmark, a snapshot, a help example or a consumer: (1) `setData` handed a factor whose shared levels sit
at other positions (a level dropped before the last, the levels reordered, an ordered factor included):
its rows now carry their own level's code, where today they carry the code of the level at their
position; (2) a `copy()` or a reload after `setData` handed a factor declaring other levels: it now has the
sampler's level count and continues as the original does, where today it is refused or runs on another
count; (3) a label the sampler lacks, on a training or test row, wherever the count fits (renamed levels,
or one level replaced by another): refused, where today it is read by position with no message; (4)
predictor columns handed in another order under their names: read by name, where today they are read by
position.
window: R surface, file-disjoint from the engine slices; after forest-kind-by-class push 2 or before it,
rebase only (both touch [R/dbarts.R](../../R/dbarts.R) and the sampler Rd). Before TODO
`state-frame-prior`'s record slices, which it leaves less to do.
budget: planned ~530 lines (R ~100, tests ~360, help, NEWS and notes ~70); upper figure 1000. Stops:
[Budget and stops](#budget-and-stops).

## Goal

`$setData(newData)` reads each predictor column of `newData` (training and test) by name and each factor
column by label against the level table the sampler was created with, as `setPredictor`,
`setTestPredictor` and `predict` already read one, and stores the data object re-coded onto that table.
A factor declaring fewer levels, the same levels in another order, or more levels that no row has, goes
in with each row at its own level; a row at a label the sampler lacks is refused by name. The engine's
level table and the record a re-creation reads are then always the same, so a copy and a reload continue
as the original does and every install of the original's stored state goes in.

## Context

Rulings, in the register's words unless quoted as the maintainer's:

- dec-A115, the maintainer on 2026-10-05: "Declared levels." A factor's split grid and level count come
  from its declared levels, so "an unobserved interior level stays splittable and predictable", and the
  grid does not "move whenever a subset missed a level".
- dec-A85: "A factor level in predict's new data that training never saw is refused by name ... as base
  R's predict.lm and randomForest refuse it." dec-B155, the maintainer on 2026-10-01 ("Use option 1"):
  the indicator route "recodes a test or newdata factor against them by label, so a test factor carrying
  only some of the levels predicts as the full-level one would".
- dec-B279 and dec-B286, the maintainer on 2026-10-07: "Refuse numbers for a factor column." The joint
  update "takes a factor or labels, matched to the column's levels by name". dec-B359, the maintainer on
  2026-10-08: "matrices should only work when matrices were given as an input"; its rule for
  `setPredictor` with no column: a data frame "coded by label as the column form is; a numeric matrix is
  taken only where every predictor column is numeric ... and refused by name on a sampler with a factor
  column". dec-B360, the maintainer the same day: "A. Apply the same rule to setTestPredictor and
  predict"; the rule: "every route that takes new predictors, setPredictor, setTestPredictor and predict,
  takes a numeric matrix only where every predictor column is numeric, judged by the column types, and
  refuses it by name on a fit with a factor column". Neither names `setData`; this plan extends dec-B360
  to it.
- dec-B346, the maintainer on 2026-10-08: "Take both on a swap", of a basis, where a factor level no row
  holds keeps its part "drawn under its prior until a row returns". For a predictor the analogy is the
  level kept in the table; a predictor level with no rows empties leaves instead, which the swap merges
  (measured below).
- dec-B231, the maintainer on 2026-10-05: "You don't want to update derived data in an ambiguous case. It
  needs to be explicit, because a caller might not understand that setData is changing a different part
  of the sampler." Its rule re-derives at a plain `setData` the values derived from the rows (split
  points, standardization, response range); a level table is not derived from the rows (below).
- dec-B195, the maintainer on 2026-10-02: "To the degree possible, state != model", and "state (and
  scratch) describe the current literal state of an MCMC sampler as it walks the posterior, as defined by
  a model and data." dec-A210 reads dec-B231 and dec-B418 as "a state carrying no data" and has a state
  between samplers whose levels differ in order go by position; it also notes "a refusal can be lifted
  after the release where an accepted meaning cannot be withdrawn". dec-B449: "a live setState never
  writes the sampler's cut points"; TODO `state-frame-prior` moves the grid off the state.
- dec-B230: a state from other data is "neither supported nor prevented".
- dec-A209 scoped this out: "whether setData gives the live sampler the new table or a re-creation takes
  the table from the state (about 30 engine lines, a draw change for such copies)".

Code as it stands:

- A data object made from a data frame keeps each factor as one column of 0-based codes into its own
  level table, `attr(x, "factor.levels")` ([`makeCategoricalModelMatrix`](../../R/utility.R)).
- A creation reads each factor's level count from that table
  ([`readDeclaredCategoryCounts`](../../src/R_interface_bartcore.cpp)); an ordered factor's grid is its
  count less one. The count is fixed for the engine's life.
- [`bartcore_setData`](../../src/R_interface_bartcore.cpp) refuses a change of column kind ("the kind
  fixes the grid, so a replacement that changes it is a different design rather than new values for this
  one"), refuses a sparse or mixed-sparse store, refuses a code at or past the engine's count, and
  installs the codes as given. The R side ([`bartcoreSamplerSetData`](../../R/bartcore.R)) already
  overwrites the handed object's `n.cuts` and `sigma` with the sampler's, and stores the result as the
  record. So the live engine reads the new rows by position in the new table against its old count, and
  the record declares the new table; columns are taken by position, whatever their names.
- A copy, a reload and a dead pointer's `setState` re-create from the record
  ([`recreatePointer`](../../R/dbarts.R)), so they take the new table's count, then install the stored
  state, whose trees hold category masks and ordered thresholds by code.
- `predict` re-codes a data frame against the record's table by label and matches columns by name
  ([`mapFactorColumnsToTrainingLevels`](../../R/utility.R)); `setPredictor` by column codes labels the
  same way ([`codeCategoricalColumnUpdate`](../../R/bartcore.R)), and a whole-frame `setPredictor`
  matches columns by name too ([`codePredictorFrame`](../../R/bartcore.R)). Both read an ordered factor
  in another order by label. `setData` is the one route that does neither.
- The `setData` docstring promises "Saved draws are rewritten so that they predict what they did"; the
  positional reading breaks that promise.
- retired: [`heldAfterAdoption`](../../R/dbarts.R) (gone since state-frame-prior slice 1) refreshes a re-created sampler's `state` field where its
  ordered column's grid differs from the recorded one; the `setData` sequence below is its only
  documented reach.
- No consumer calls `setData` (rechecked by the critique): stan4bart (branch bartcore) re-creates its
  samplers from data it built and calls `setState`; bartCause (dbarts-1.0) calls neither; treatSens
  (master and dbarts-1.0) only `setTestPredictor`; bairrtt calls `setPredictor` on a numeric column and
  the joint update. Nothing in [dbarts.h](../../inst/include/dbarts/dbarts.h) swaps data. Outside the
  tests, dbarts calls `setData` only in the composition matrix, on subsets of a sampler's own data
  object. Of the 40 tinytest files that call `setData`, two hand it another level table (the critique's
  logging shim): test-bartcore.R (a level 11) and test-state-grid-install.R (4 of 6 levels); a third,
  test-data-code-channel.R, hands it a numeric matrix carrying the sampler's own tables.

### Measured: the tip beside 0.9-34

Scripts in scratch/fls of this worktree (p2.R, sequences; p3.R to p8.R, corners; crit/, the critique's
probes rerun), on a build of f514c1d8 (scratch/libs/fls) and on 0.9-34 (scratch/libs/cran). 200 rows, a
numeric column `a` held fixed across data sets and a factor `f` of levels a to f, the response `a` plus 0,
2, 4, 0, 2, 4 by level, every data set on the first one's response range; 20 trees, seeded before
creation, 30 sweeps of burn-in and 10 kept draws, one chain; then `setData` with the levels in the first
column. "Draws as the original" means one sweep from the same R seed gives the original's training fits
to rounding (at most 1.8e-15, the fits re-summed) or identically. "By position" means the kept draws'
prediction at each label after the swap is the one they gave the level at that position before it: 3.83
to 4.04 away from the same label, and as close to the position as a swap of the same levels gives (0 on
the tip; 0.073 on 0.9-34 with NA, where rows dropped for NA move the cut points of `a`). Where the two
readings agree the row says so. `setState(stored)` hands back the state `setData` stored.

| Levels handed over | 0.9-34 (both kinds, with or without NA) | tip, unordered (with or without NA) | tip, ordered | tip, ordered with NA on it |
|---|---|---|---|---|
| a..f, the same | taken; copy and reload draw as the original; `setState` NULL everywhere | taken; copy and reload draw as the original; `setState` TRUE everywhere | the same | the same |
| a..d, last two dropped | refused: "maximum number of cuts must be of length equal to 5" (a factor is indicator columns) | taken (label = position); copy and reload refused, "state is not consistent with this sampler"; live `setState` TRUE | taken; copy and reload made on 3 points, draws differ by up to 0.14; their `setState(stored)` refused except a reload's first call | as unordered: copy and reload refused |
| a, c..f, b dropped | refused (length 6) | taken, by position; copy and reload refused | taken, by position; copy and reload made, draws differ by up to 0.38; `setState` as above | taken, by position; copy and reload refused |
| a, c..g, b replaced by g | taken, by position; copy and reload as the original | taken, by position (g read as f, no message); copy and reload as the original | the same | the same |
| a..h, g and h on no row | taken (0.9-34 drops a level no row has); copy and reload draw as the original | taken; copy and reload made on 8 levels, draws differ by up to 0.30 (0.50 with NA) | taken; copy and reload made and draw as the original; their `setState(stored)` refused except a reload's first call | the same |
| a..h, rows at g and h | refused (length 9) | refused: "categorical predictor values must be existing category codes" | refused: "ordered factor predictor values must be existing level codes" | the same |
| f..a, reordered | taken, by position; copy and reload as the original | taken, by position; copy and reload as the original | the same | the same |
| A..F, renamed | taken, by position (label = position) | taken, by position; copy and reload as the original | the same | the same |

On every tip row: `predict` reads a new frame by label against the new table, so it agrees with the live
sampler's positional reading, and refuses a label the new table lacks ("test data factor 'f' has levels
not present in the training data") although the engine still holds that level; a reload predicts as the
live sampler. After the two "dropped" rows a later `setData` back to the six levels, with rows at all six,
is taken by the live sampler and refused by its copy (the copy has fewer levels, or cannot be made). The
critique's twin probe (crit/q1c.R, seeded before creation) confirms the reading: the sampler handed the
b-dropped levels draws bit for bit as a twin handed the same rows relabeled by position, and 0.35 away
from one relabeled by label.

An example of "by position" (p8.R): with b dropped, the kept draws predict 2.63 at level c after the swap,
the value they learned for b (4.66 before it). The live chain re-learns fast in this fixture: one sweep
later its mean fit at the c rows is 4.48 against a truth of 4.51. So the lasting harm is to the kept
draws, which no longer "predict what they did", and to copies and reloads; inside a larger sampler, a
chain pushed off its posterior at every swap whose table moves.

Corners on the tip:

- A data object made from a matrix of codes declared categorical (`attr(x, "varTypes")`) with no level
  table, handed to a sampler made from a data frame: taken, the record then has no level table, and
  `copy()` is refused (p3.R). A sampler made from such a matrix, handed codes 0 to 3 of its 0 to 5:
  `copy()` refused, the re-creation inferring the count from the codes present (p4.R).
- The test set of a data object with other levels is coded against that object's own training table, so
  it reads as its training rows do (p3.R).
- A character column rebuilt from a draw that missed b (`as.factor` gives sorted levels with rows only):
  taken, its c rows coded 1, as b (crit/q3.R; 0.9-34 refuses it for its column count).
- Columns handed in another order under their names, `y ~ g + f` to a `y ~ f + g` sampler with yes/no
  labels on both: taken by position, the record's names become "g f", and the fit at f = yes / no is 0.03 /
  2.98 against a truth of 3 / 0 (crit/q4.R). A whole-frame `setPredictor` matches them by name.
- A sparseFactor column reaches the bridge and is refused there, "requires a dense predictor matrix"
  (crit/q5.R).
- A predictor factor with an NA level is created; `predict` at a row of that level is refused as a level
  training never saw (p5.R): the by-label helper matches `as.character`, which loses an NA level.
- Rules on a level whose rows leave (crit/q6b.R): the swap merges a leaf it empties, under every option
  and today with identical tables. 84 nodes before; 82 after a swap to the same data less its b rows
  (and 82 after the by-label swap to new rows with no b); 78 less its c rows; 84 after new rows at all
  six. The kept draws still predict at b after the by-label swap (2.54).

Emulating the recommended rule on the tip (p2.R `recode`: the new factor re-coded by label onto the
sampler's table before `dbartsData`; it is the rule's level reading, without its column-name match and
its matrix cases), on all four kinds: every taken row draws as the original in copy and reload (at most
8.9e-16), `setState(stored)` is TRUE on the live sampler, the copy and a reload as first and later call,
the kept draws predict by label, `predict` at a level the new rows lack is taken (2.58 at e, the value
learned for e), a later `setData` back to the six levels is taken by the copy, and the rows at g and h,
the replaced level and the renamed levels are refused. The reordered ordered factor is taken by label
and its copy draws as the original.

### What base R does

Measured with lm (scratch/fls/lm.R). `predict.lm` reads a new factor by label against the fit's levels:
fewer levels, another order (an ordered factor in reverse order predicts 0.03 2.02 4.01 0.97, as the
original) and an extra declared level no row has all predict as the full factor would; a new label on a
row is refused ("factor f has new levels e"), renamed levels too. `update(fit, data = )` refits from
scratch and takes the new data's levels: that is a new model, the analogue of a new sampler, not of a
swap.

## What a user means

`setData` gives the same model new observations: inside a larger Gibbs sampler that imputes a predictor or
moves rows, once per sweep. A factor there is a variable whose possible values were fixed when the model
was stated; a row's value is its label. The level table a data object carries is not a statement the
caller makes about the model: every factor carries one, and it shrinks or reorders whenever a factor is
rebuilt from characters (`factor(x)` on a draw that missed a rare level drops it; a character column's
levels are sorted). A user handing such data means "these rows, at these levels", never "level c is now
level b" or "the variable now has four values". So the correct model keeps the sampler's levels, reads
each row by label, keeps a level no row has in the table (dec-A115) so that it keeps its side in every
rule that survives the swap and can still be predicted, and refuses a label the model never had, as
`predict` does. A rule whose one side held only that level's rows is merged at the swap, as it is today
with identical tables; that is the same under every option. A different set of values is another model,
which is a new sampler.

dec-B231's rule re-derives at a plain `setData` the values derived from the rows, and uses values a data
object carries. A level table is neither: it is not computed from the rows, and taking the one every
factor carries as a deliberate record would make every `setData` change the sampler's levels, the
ambiguity the ruling forbids. The precedent is in `setData` already: it refuses a change of a column's
kind as "a different design rather than new values for this one", and it overwrites the handed object's
`n.cuts` and `sigma` with the sampler's. The level table is the same kind of thing, the design the
sampler was made with. Where dec-B231 moves a numeric split onto a re-derived grid by value, so that a
threshold keeps its meaning, the factor's analogue is reading by label, which this plan does by
re-coding the rows and not the trees.

## Options

**A. setData gives the live sampler the new data's table, rows by position (the TODO's first).** As a
model: no, while rows are read by position; a level dropped before the last moves every later one, and
one replaced reads as another. Draws: the live sampler's change on every swap whose count differs (the
tree prior counts levels, and rules on a vanished level are repaired). Copy, reload and state: agree with
the live sampler. Engine ~200 lines: a per-column count through `SamplerBase::setData` (internal virtual,
`--preclean`), the store's counts, codes and an ordered grid resized, a mask's word count and the
missing-value code moved where the count crosses 64, live trees with rules on a vanished code repaired
through the forced-update path, kept draws' masks trimmed. Bridge and R ~30. dbarts.h untouched.
Consumers: none calls `setData`. A dropped level leaves every rule and `predict`; if a later table brings
it back, it has no side in any rule.

**A2. The same, with the trees re-mapped by label.** As a model: yes, for a variable whose values really
changed. Draws: as A, and rules move with their labels. Engine ~300 to 450 lines beyond R's map: A's, plus
every live and kept category mask permuted, ordered thresholds carried through a map that must be
increasing (another order refused), and a rule for the side a level new to the table takes, which needs
a ruling. A level missing for one sweep is gone from every rule, as under A; the difference from C is that
when it returns it needs that ruling, where under C it kept its side.

**B. A re-creation takes the level count from the state (the TODO's second).** As a model: the live
sampler is left as it is, so rows stay read by position. Draws: copies and reloads after such a swap now
continue as the original. Engine ~30 lines (an adopting install sets a factor column's count from a
record), bridge and R ~40 (an unordered factor's count recorded in the state, which records nothing for
it today; an ordered one's is its recorded grid). It adds data to the state against dec-B195, dec-B449
and TODO `state-frame-prior`, and the state's format changes: stan4bart stores states in its fits, so a
fit saved before must still restore (an absent record read as today).

**C. setData reads columns by name and a factor by label against the sampler's table (recommended).** As
a model: yes, the one above. Draws: the four sequences under `rng`; bit for bit elsewhere. Copy, reload
and state: the record's table is the engine's, so all agree with the live sampler and the stored state
installs everywhere. Engine 0 lines; bridge 0; R ~100. dbarts.h untouched. Consumers: none. A categorical
column with no level table, which the tip takes, is refused on a sampler with a factor (dec-B360). NEWS
states it as a fix (the saved draws now "predict what they did") and the change from 0.9-34: other level
counts taken; reordered, renamed or replaced levels, which 0.9-34 read by position, read by label or
refused.

**D. setData refuses a factor whose level table differs from the sampler's.** As a model: safe, no wrong
reading reachable. Draws: none move; sequences that run today are refused. Copy, reload and state: never
reach a mismatch. R ~15 lines. The cost is the larger sampler's: a loop that rebuilds a factor each sweep
stops the first time a rare level is not drawn, and the caller must write `factor(x, levels = <the
sampler's levels>)` on every call, which is what C does for them; `setPredictor` and `predict` would
take what `setData` refuses.

## Recommendation

C. It is what a user means and what base R does with a factor's levels; it gives `setData` the rule every
other route already has (dec-A85, dec-B155, dec-B279, dec-B359, dec-B360), keeps dec-A115's declared
levels through a swap, and needs no data in the state (dec-B195, dec-A210). It is the cheapest correct
option: no engine, bridge or ABI change, and the TODO's draw change for copies comes with it. A and B leave
a level read as another; A2 is correct for a different variable but costs an engine slice and a ruling;
D is safe but refuses what the other routes take.

## The rule

At `$setData(newData)`, before anything is installed:

1. Columns are matched by name where both the sampler's predictors and `newData`'s carry unique names,
   as `codePredictorFrame` and `predict` match them, and put in the sampler's order; names that differ
   from the sampler's are refused, naming them (as built: read by position where the names differ, see [Built](#built)). Where either side has no names, by position, as
   today.
2. Each column the sampler records a level table `T` for (from `attr(sampler$data@x, "factor.levels")`)
   is judged by the table `newData` records for it, `N`, however the column is stored (a factor in the
   container, or a numeric matrix column carrying `factor.levels`), in the training and the test
   predictors:
   - `N` identical to `T`: nothing is touched, the codes go in as handed and the bridge's code check
     stays the refusal of a code outside the table, in its words.
   - `N` another table: each level of `N` is matched to `T` by label, an NA level to an NA level of `T`.
     A level of `N` that `T` lacks and that some training or test row has is refused, naming the column
     and the labels, in the words of `setPredictor` by column ("column 'f' has label 'g' not among its
     training levels"); one no row has is dropped. A missing value stays missing. The column's codes
     become each row's level's position in `T`, and its table `T`. An ordered factor is read the same way,
     whatever order `N` gives its levels: the sampler's order stands.
   - No `N`, a categorical column of bare codes: refused by name, in dec-B360's words, pointing to a data
     frame.
3. A column whose source is sparse (a negative entry of the container's map) is left to the bridge,
   whose refusal of a sparse store stays as it is; the rule does not read it.
4. The data object the sampler stores carries the sampler's tables and names (as built: true where columns
   are matched by name; a positional read keeps the new data's names, see [Built](#built)).

A column the sampler records no table for (a sampler made from a matrix of codes) is left as today
([Out of scope, and where it goes](#out-of-scope-and-where-it-goes)).

## Constraints

- `bartcore_setData`, the engine and dbarts.h have no diff; the bridge's code check stays the backstop.
- Identical tables and names are the code that runs today: the record is the handed object, bit for bit.
- A refusal leaves the sampler, its record and its `state` field as they were.
- The handed object is the caller's: R's copy semantics leave it untouched; only the stored record is
  re-coded.
- `heldAfterAdoption` stays; with C it is reached only by a record assigned by hand or an object saved by
  an earlier build of this branch ([Calls made in planning](#calls-made-in-planning)).
- No new argument, no warning. The by-label match and its refusal reuse the words of
  `codeCategoricalColumnUpdate`; the column match follows [`alignContainerFactorLevels`](../../R/utility.R)'s
  dense arm, but matches levels and not `as.character` of rows, so an NA level keeps its name.
- Twin tests seed before creation: chain seeds are drawn from R's stream when a sampler is created
  ([docs/architecture.md](../architecture.md), its RNG section); seeded only before `run`, a session's
  first sampler draws 0.83 away from later identical ones (crit/q2.R).

## Steps

1. The re-coding (~100 lines, [R/bartcore.R](../../R/bartcore.R), or [R/utility.R](../../R/utility.R)
   beside its siblings). A function (working name `alignDataFactorLevels(newData, x.train)`) applying
   [The rule](#the-rule) to `newData@x` and `newData@x.test` (a dense matrix, or the mixed container,
   whose sparse-mapped columns it skips) and returning the re-coded object. Called in
   [`bartcoreSamplerSetData`](../../R/bartcore.R) after `checkBasisRecord` and before the `.Call`, so
   every refusal precedes the swap; the record stored is its result.
2. Tests ([Tests](#tests); ~360 lines). A new `inst/tinytest/test-setdata-factor-levels.R`. Rewrite, not
   delete: in [test-state-grid-install.R](../../inst/tinytest/test-state-grid-install.R) the block "an
   ordered factor after a whole-data replacement declaring fewer levels" (about 105 lines) pins today's
   rule (a copy on 3 points, the original's state refused). Under C that sequence gives a copy on 5 points
   that takes the original's state; the block's refresh pins move to the same fixture with the record then
   assigned by hand (`sampler$data <- <the data handed to setData>`, `n.cuts` and `sigma` copied), which
   recreates the mismatch C closes. In [test-bartcore.R](../../inst/tinytest/test-bartcore.R) the refusal
   of a row at level 11 expects "existing category codes"; it becomes the by-label text.
   [test-data-code-channel.R](../../inst/tinytest/test-data-code-channel.R)'s matrix with the sampler's own
   tables and a code past the end keeps the bridge's text and is not edited.
3. Help, NEWS and notes (~70 lines). The `setData` docstring in [R/dbarts.R](../../R/dbarts.R): one
   sentence on the rule and its refusals. [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd),
   Saving: the sentences on `setData` to an ordered factor declaring fewer levels, and on the original's
   state refused by its copy, go; in their place, a copy and a reload hold the sampler's levels.
   [dbartsData.Rd](../../man/dbartsData.Rd): a swap reads the factor by label. NEWS: the entry "New
   predictors for a fit with a factor column come as a data frame, coded by label" gains `setData`, as a
   fix (the kept draws now predict what they did, and a copy continues as the original) and with the
   change from 0.9-34 (other level counts taken, columns matched by name, reordered levels read by label,
   renamed or replaced levels refused). [install-surface.md](../design/install-surface.md): the paragraph
   on a factor column's level table across `setData` is rewritten to the rule;
   [data-store.md](../design/data-store.md) one sentence. The ledger entry for the calls made, at landing.
   TODO: drop `factor-levels-across-setdata`; add the two lines under
   [Out of scope, and where it goes](#out-of-scope-and-where-it-goes).

## Tests

tinytest only; tests/cpp has nothing to pin. Each row: what it pins; the mutant it must kill. Fixtures:
p2.R's data sets, 20 trees, kept draws, every sampler seeded before creation; "draws as the original" is
`expect_equal` of five sweeps' fits from one seed (the fits are re-summed at a re-creation, 9e-16) and
`expect_identical` of the trees.

| Test | Pins | Mutant |
|---|---|---|
| same table neutral | the record is identical to the handed object; next draws identical to a twin's, both seeded before creation, given the object directly | the object re-coded where the tables are identical (non-identical attributes) |
| fewer, tail | unordered and ordered: taken; the record's table is the sampler's six; copy and reload made and draw as the original; `setState(stored)` TRUE on the live sampler, the copy, a reload's first and a later call | the record keeps the handed table (copy refused, unordered) |
| fewer, middle | b dropped: each row's code is its label's position in the sampler's table; kept draws predict at c what they predicted at c before | codes left as handed (by position) |
| character draw | a character column rebuilt without b (`as.factor`): its c rows coded as c | (pins the imputation shape) |
| reordered | unordered and ordered f..a: by label as above; copy as the original | tables compared by length, not identity |
| replaced | a, c..g with rows at g: refused naming g; sampler, record and field untouched | unknown labels caught only by the bridge's code-range check |
| more, no rows | g and h declared, no row: taken; record table six; copy as the original | an unused unknown label refused |
| unknown on a row | rows at g and h, and renamed A..F: refused naming the labels; sampler untouched | the refusal dropped (unknown labels become missing values) |
| unknown on a test row | a test set alone holding g: refused | the unknown-label check run on training rows only |
| test set | new training table with b dropped, test set declaring c and e only: test fits as by label | x.test not re-coded |
| columns by name | `y ~ g + f` to a `y ~ f + g` sampler: record in the sampler's order, fits as by name; names that differ refused (as built: read by position where the names differ, see [Built](#built)) | columns taken by position |
| back again | after the tail case, `setData` of all six levels with rows at each: taken by the live sampler and by its copy | (pins the closed TODO sequence) |
| NA level | a sampler made with an NA level, handed the levels reordered: the NA level's rows keep it, a missing value stays missing | match through `as.character` (NA level rows become missing) |
| missing values | the ordered NA fixture of p2.R, tail case: copy and reload made and draw as the original | (pins the TODO's missing-value reach) |
| codes with no table | a sampler made from a frame, handed a categorical code matrix with no table: refused in dec-B360's words, sampler untouched | the no-table check dropped |
| matrix with the tables | test-data-code-channel.R's matrix carrying the sampler's tables, a code past the end: the bridge's text, unchanged | every numeric matrix refused |
| sparse left alone | a sparseFactor column: the bridge's "requires a dense predictor matrix", unchanged | the helper reading a sparse-mapped column |
| code sampler untouched | a sampler made from a code matrix, handed codes: taken as today | the rule applied to a column with no recorded table |
| refresh by hand | the rewritten block of test-state-grid-install.R, the mismatch made by hand: its refresh pins as before | (keeps the refresh's mutants of the state-grid-install rounds killed) |

## Verification

The tier "Changes draws", nothing under src/:

- `R CMD INSTALL --preclean -l <lib> .` and the full tinytest suite with `R_LIBS=<lib>`.
- Reference build: the four `test-reproducibility-*.R` files pass unedited, and the three compares print
  identical draws for every scenario and no `max |z|` line: 55 against `equivalence-308c45be.rds`
  (`--bitwise --strict-coverage`), 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. None calls `setData` with another level table.
- The 42 commands of `.github/workflows/exact-gates.yaml` in quick, unchanged, as the class requires.
- `Rscript benchmarks/R/composition-matrix.R`, no disagreement (its `setData` probes hand a sampler its own
  data's subsets, same tables).
- `src/` and dbarts.h have no diff (`git diff --stat`).
- Lint chain, each on its own exit status: `lintr::lint_package()`, `air format --check .`,
  `tools/check-rc-codoc.R`, `tools/check-win-drift.R`, `tools/check-doc-freshness.R`; NEWS parses;
  `R CMD check --as-cran` from a clean tarball.
- stan4bart's suite against the slice's library, expected unchanged; bartCause's, which never calls
  `setData`, at the orchestrator's discretion.
- The x86-64 leg as the orchestrator's mechanics require; no engine diff, so its sanitizer half is not
  owed.

## Budget and stops

| Step | Lines |
|---|---|
| 1 R | ~100 |
| 2 tests | ~360 |
| 3 help, NEWS, notes, TODO | ~70 |

Estimates here have run 1.5 to 2 times low. Stop and report, building nothing further: if a bitwise
compare prints a `max |z|` line or a snapshot file moves; if any change under src/ or dbarts.h appears
needed; if a tinytest file outside step 2's two edited files fails (list it, do not edit it); if
stan4bart's suite moves; at 150 R lines or 1000 in all.

## Interactions

- forest-kind-by-class push 2: [R/dbarts.R](../../R/dbarts.R), the sampler Rd and NEWS; rebase only. Its
  `checkBasisRecord` call in `bartcoreSamplerSetData` comes first, as now.
- TODO `state-frame-prior`: its record of the grid on the data object is the numeric columns' business;
  after this slice a factor column's grid already follows the record, and the record the sampler's table.
- TODO `warm-start-salvage` and the warm start: a donor's state goes in by position, dec-A210 unchanged.

## Out of scope, and where it goes

- A sampler made from a matrix of codes records no table, and a re-creation infers its count from the
  codes present, so `copy()` is refused after `setData` of fewer codes (p4.R); and such a sampler handed a
  frame gains the frame's table in its record (test-bartcore.R, seen in the critique's shim), which a copy
  then reads as its count. New TODO line `code-matrix-category-count`: record the count such a sampler
  was made with, and decide what a frame handed to it means.
- `predict` and `setPredictor` match a factor through `as.character`, so a predictor's NA level is
  refused at `predict` (p5.R), and at `setPredictor` a missing value would read as that level. New TODO
  line `na-level-predictor-doors`, beside dec-A210's refusal of an NA level for a basis.
- A state moved between samplers whose levels differ: by position (dec-A210), unchanged.

## Calls made in planning

- An ordered factor whose levels come in another order is read by label, the sampler's order standing,
  as `predict`, `setPredictor` and base R's `predict.lm` read one (settled by the orchestrator,
  2026-10-17, from the critique). The alternative refuses it, which a later release could lift (dec-A210:
  "a refusal can be lifted after the release where an accepted meaning cannot be withdrawn"); its cost
  is a refusal only this route makes.
- Columns are matched by name where both sides carry names, and names that differ are refused (settled by
  the orchestrator, 2026-10-17, from the critique; as built: read by position where the names differ, see
  [Built](#built)), as a whole-frame `setPredictor` and `predict` match them; about 5 lines and a test. The alternative, a TODO line, leaves a correctly read factor in the
  wrong column (0.03 / 2.98 against 3 / 0 in crit/q4.R).
- A column is judged by the level table it carries, not by how it is stored: a numeric matrix column
  carrying the sampler's own table goes in as today; only a categorical column with no table is refused,
  extending dec-B360 to `setData`. The alternative refuses every numeric matrix, as the other routes'
  `refuseMatrixOfCodes` does, which refuses test-data-code-channel.R's matrix before the bridge names its
  bad code; the cost taken is one route reading a tabled matrix the others refuse.
- A declared level no row has and the sampler lacks is dropped, not refused, as `predict.lm` takes one.
  The alternative refuses it, at the cost of a loop whose factor declares a superset of levels.
- An NA level is matched by name. The alternative refuses a row at an NA level with dec-A210's remedy;
  the cost taken is none at this route, creation accepting such a factor.
- `heldAfterAdoption` is kept and its pins move to a mismatch made by hand. The alternative removes it
  (about 10 lines and five pins), reached then by nothing documented.
- The record stored is the re-coded object, not the handed one, where the tables or names differ. The
  alternative, the handed object, is what makes a copy differ today.

## Open calls

None. The maintainer ruled option 5 on 2026-10-17, "For decision 2, yes, use your recommendation"
(dec-B454); the plan is written for it (C).

## Built

Built on 92ae2714 (origin/bartcore), after the ruling (dec-B454). R only: `src/` and dbarts.h have no diff.

- Step 1: [`alignDataFactorLevels`](../../R/bartcore.R) with three helpers beside
  `codeCategoricalColumnUpdate`: `refuseNewDataLabels`, `recodePredictorColumn` and
  `permutePredictorColumns`. It is called in `bartcoreSamplerSetData` after `checkBasisRecord` and
  the column count check, before anything is stored, and the record is its result.
- Step 2: inst/tinytest/test-setdata-factor-levels.R, the block of test-state-grid-install.R
  rewritten to a record assigned by hand, and the level-11 refusal of test-bartcore.R moved to the
  by-label text.
- Step 3: the `setData` docstring, the `newData` item and the Saving paragraph of the sampler Rd,
  dbartsData.Rd, NEWS, install-surface.md, data-store.md and TODO, as the step lists.

Calls made building, each where the plan was silent or wrong:

- Names that differ are read by position, not refused. The rule said "names that differ ... are
  refused", but test-sampler-setData.R hands `y ~ x.new` to a `y ~ x` sampler, and
  `codePredictorFrame`, which the rule cites, matches by name only where the frame carries every one
  of the sampler's names and reads the rest by position. The refusal would have failed a file
  outside step 2, a stop. The orchestrator's call, columns matched by name where both sides carry
  them, stands where the names are the sampler's.
- The refusal of a label the sampler never had is one function, `refuseNewDataLabels`, called once
  per column with the training and test labels, with its own tests, so that it can be changed alone
  if the maintainer rules otherwise on new levels (the orchestrator's instruction). No acceptance of
  new levels is built.
- A sparse source, training or test, is left whole to the bridge, which refuses any sparse store;
  the rule's column-by-column skip gives the same refusal.
- A number in a tabled matrix that is no code of its own table becomes -1, so the bridge refuses it
  in its own words; read through the table it would have become a missing value.
- A test set is read through the training table of the same data, which dbartsData codes it against;
  dbartsData leaves no table on a test set, so a branch reading one was dead and was removed in the
  fix round.
- Columns put in the sampler's order take the sampler's design attributes (`term.labels`, `terms`,
  `drop`, `indicator.levels`), so the record carries the sampler's names (rule 4; as built: true where
  columns are matched by name; a positional read keeps the new data's names).
- The no-table refusal is worded for a data object: "the data hold factor 'f' as numbers with no
  level table; give dbartsData a data frame, whose factor columns are coded by label".
- The "same table neutral" test pins the record's identity only. Its mutant, re-coding where the
  tables are identical, is equivalent: the re-coded object is identical to the handed one.
- The "matrix with the tables" row is test-data-code-channel.R's pin, unedited; the new file adds a
  tabled matrix in another order, read by label, and a stray code refused by the bridge.

Gates on the tip: the full suite in one process, the three bitwise compares on the reference build (55, 15
and 11 scenarios identical, against equivalence-32418b6e.rds, the baseline current at 92ae2714), the four
snapshot files, the 38 commands of exact-gates.yaml's main job loop (its monotone job's other 4 not run:
they reach only the engine, which has no diff), the composition matrix, the lint chain and `R CMD check
--as-cran`. Of the plan's mutants and the build's own, all were killed but the re-coding of identical
tables, which is equivalent. The review's survivors were closed in the fix round with tests of two factor
columns with other tables swapped by name, a numeric and a factor column swapped by name, a tabled matrix
holding missing values and a sparse test set; the full suite, the lint chain and `R CMD check` were rerun
on its tip.

### Follow-up after landing

A review after landing found that columns put in the sampler's order took the sampler's `terms`, whose
predvars were fixed from the creation data, while the columns held the new data's basis: with
`y ~ poly(a, 2) + f` handed as `y ~ f + poly(a, 2)`, `predict` was up to 0.71 from an in-order twin, and
`getTrees(newdata = )` and a frame `setTestPredictor` were off too (`ns()` the same). The sampler's terms
now take the new data's predvars, matched by variable (`permutePredictorColumns`); the record, `predict`,
`getTrees(newdata = )`, a frame `setTestPredictor` and a seeded run are then a twin's, with new data of
another row count. `pdbart` over the factor agreed before and after; over a `poly()` variable it refuses
the column on any fit.

- Read by position, a name of the sampler's that sits in the data at another place is refused, naming
  the data's column and the sampler's name in its place (settled by the orchestrator, 2026-10-17, as more
  faithful to dec-B454's "each column by its name"): a sampler of (a, f, g) handed (a, g, h) read g as f
  with no message. Data whose names share none of the sampler's in another place, `y ~ x.new` for
  `y ~ x`, are still read by position.
- A number read in a factor's place names the data's column and the place, and tells a caller who
  handed a data frame to give a factor, keeping the data-frame advice for a matrix.
- The help says how to hand over a factor relabeled on purpose (`levels<-` to the sampler's labels, or a
  new sampler) and that a positional read gives the record the new names, which `predict` then needs.
- Tests: the `poly()` and `ns()` swap, the two refusals, and a swap on a sampler of three chains and two
  threads under case weights and on a probit sampler under a 0/1 mask.

## Revised from the blind critique

The critique (scratch/flsc-crit.md, probes in scratch/flsc/, rerun on this worktree's build in
scratch/fls/crit) kept the recommendation. What changed:

1. The rule judges a column by its table, not its storage; only a categorical column with no table is
   refused; test-data-code-channel.R is named and pinned.
2. Option 5 no longer claims to keep a rare level's rules: the swap merges an emptied leaf under every
   option (84 to 82 nodes); C keeps the level's side in the rules that survive and its prediction. dec-B346
   is cited as an analogy.
3. dec-B360 is cited; the matrix refusal extends it.
4. The replaced-level row is measured and added (tip and 0.9-34, all kinds), with rng item (3), a test and
   its mutant.
5. An ordered factor in another order is read by label (the orchestrator's call); the refusal, its test
   and its mutant are gone.
6. Columns are matched by name, refusing names that differ (the orchestrator's call), with a test (as built: read by position where the names differ, see [Built](#built)).
7. Sparse-mapped columns are left to the bridge, pinned.
8. Tests gain a test set whose training table differs, a test-only unknown label and a character column.
9. Every sampler in the probes and tests is seeded before creation; every p2.R number was rerun so.
10. The maintainer text gives option 3 its bridge and R lines and the state-format cost, and option 5 its
    matrix refusal.
11. The dec-B231 argument cites the kind refusal and the `n.cuts` and `sigma` overwrite.
12. The preamble's "by position" allows 0.9-34's 0.073 with NA; 0.9-34 with NA on an unordered factor was
    run and matches.
13. The emulation is said to be the level reading without the name match and the matrix cases (rule 4 is
    gone).
14. The code-matrix sampler handed a frame joins `code-matrix-category-count`.
15. The docstring's promise that saved draws predict what they did frames NEWS as a fix.
16. The harm sentence says the live chain recovers within a sweep here (4.48 against 4.51) and the
    lasting harm is to kept draws and copies.

## For the maintainer

Ruled 2026-10-17: option 5 (dec-B454). The text as shown:

A sampler can be handed new data in the middle of a run, which a larger model does every sweep, for
example when it imputes a predictor. When a predictor is a factor, the new data's factor carries its own
list of levels, which need not be the sampler's: a level no row has may have been dropped, which happens
when a draw misses a rare level and the factor is rebuilt; the levels may come in another order, since a
text column's levels are sorted; and a level the sampler never had may appear.

Today the running sampler keeps its own list but reads each new row by its position in the new list, not
by its name. If "b" is dropped, rows at "c" are read as "b": a stored draw that predicted 4.66 at "c"
before the swap predicts 2.63, what it learned for "b", after it. If "b" is replaced by a new level "g",
rows at "g" are taken as "f" with no message. The running chain re-learns quickly: one sweep later its fit
at the "c" rows was 4.48, against a truth of 4.51. The stored draws stay wrong, though. A copy or a saved
and reloaded sampler is rebuilt from the new list, so it either cannot be made or runs a different chain,
its fits 0.38 apart from the original's after one sweep in one probe. 0.9-34 refused a factor whose number
of levels with rows differed, and otherwise read levels by position, whether reordered, renamed or
replaced.

Under every option, when a level's rows vanish for one swap, tree rules whose one side held only that
level's rows are merged, as today: 2 of 84 tree nodes went in one run. The options differ in what happens
to the level itself.

The question is which list governs when new data brings a factor with other levels.

1. The new list replaces the sampler's, and rows are still read by position. Copies then match the
   original, but "c" is still read as "b". A dropped level is gone from every rule and from prediction,
   and if a later list brings it back it has no side in any rule. About 200 lines of engine code.
2. The same, with rules moved by name. Right for a variable whose values really changed. A level missing
   for one sweep is still gone from every rule, and when it returns a new ruling must say which side of
   each rule it takes. About 300 to 450 lines of engine code.
3. A copy or a reload takes the level count from the saved state. About 70 lines, 30 of them engine code,
   and a change to what a saved state holds, which fits the stan4bart package saved before must still
   read. Copies then match, but "c" is still read as "b", and a state then carries data, against the
   rulings that it carries none.
4. Refuse a factor whose list differs. About 15 lines, and nothing that runs today draws differently.
   But a loop that rebuilds its factor each sweep stops the first time a rare level is not drawn, the
   caller must restate the sampler's levels on every call, and the other ways of handing in predictors
   accept what this one refuses.
5. Read each new row by name against the sampler's own list, and each column by its name, as replacing
   predictors and prediction already do.
   - Fewer levels, another order and extra levels no row has are all taken.
   - A level with no rows stays in the list, keeps its side in every rule that survives, and can still be
     predicted: 2.58 at "e" in one probe, the value learned for it.
   - A level the sampler never had, on a row, is refused by name. That includes renamed or replaced
     levels, which 0.9-34 read by position.
   - A factor column given as bare numbers with no list of levels is refused, as replacing predictors
     refuses one; today it is taken.
   - About 100 lines of R and 360 of tests, no engine change, nothing for the consumer packages. Draws
     change only where a level used to be read as another, or a column taken for another.

Recommended: 5. It is what a user means ("b" is "b"), it is what base R's predict does with a factor's
levels, and it gives this method the rule the other ways in already have, at no engine cost.

## Landing note

Landed 2026-10-17 as 92ae2714..6beef48d on bartcore: the plan and its revision from the critique (3dc7da24,
a51abb71), the ruling's ledger entry (78840c0a), the build (20c2ac88), the help, NEWS and notes (46c08f31),
the fix round's tests and removal (70396cb1), and the fix round's text (6beef48d). One planner, one blind
critique, one implementer, one Sonnet review whose reviewer also checked the fix round. The tier is the
header's, "Changes draws", whose review the header names as opus; the review ran on Sonnet.

Built as planned: `alignDataFactorLevels` with `refuseNewDataLabels`, `recodePredictorColumn` and
`permutePredictorColumns`, called in `bartcoreSamplerSetData` before anything is stored; no change under
`src/`, in dbarts.h or to a consumer. Departures, the calls of [Built](#built) and of the ledger entry dec-A213:
names that differ are read by position, not refused, so that the plan's lines that said "refused" are marked as
built, and the record keeps the sampler's names only where columns are matched by name; a sparse source is left
to the bridge; a stray number in a tabled matrix becomes -1 for the bridge to refuse. The plan's Verification
still names equivalence-308c45be.rds and 42 exact-gate commands; the baseline current at 92ae2714 was
equivalence-32418b6e.rds, and the 42 is the whole file, the main job's loop being 38.

Review: LAND AFTER FIXES. The reviewer found the code correct on every probe (two factor columns with
different tables, ordered and unordered, with and without NA; fewer, reversed and extra levels; test rows; a
tabled matrix; a sparse source) and ran 16 mutants: seven killed, five survivors that were gaps in the tests
(r01, r02, r03b, r06, r07), two dead code (r04, r08) and two equivalent (r09, r16). The fix round added the tests
that kill the five, removed the dead branch that read a test set's own level table (`dbartsData` never leaves
one), corrected the sampler's help, install-surface.md, NEWS and the gate count, and marked dec-B454 built.
The check of the round: LAND. r01, r02 and r07 rerun and killed, a new mutant r18 (test rows read through the
sampler's table) killed, the builder's m05 rewritten as r17 on the remaining table comparison and killed, and
the removal safe, only a hand-set slot reaching the old behavior. It asked for the notes on rule 4 and on the
Built line about the record's names, which are made above.

Gates, on the tip: the full tinytest suite 27145 results, 0 failed; lint clean; `R CMD check --as-cran` one
NOTE (the Date field); the reference build's bitwise compares at 46c08f31 identical, 55 of 55 against
equivalence-32418b6e, 15 of 15 against bcf-equivalence-1b7d730c and 11 of 11 against
multinomial-equivalence-80b1c8d4, the fix round having changed R only in the deleted branch and a comment; the
exact gates 38 of 38 for the main job's loop. Not run: the monotone job's 4, which reach only the engine; the
x86-64 leg and the sanitizers, nothing under `src/` having changed; the consumers' suites, stan4bart,
bartCause and treatSens making no `setData` call. The reviewer did not rerun the suite, the compares, the check
or the builder's mutants m01 to m21, and ran every probe on one chain, one thread and a Gaussian response.

Size against the plan: planned about 530 lines and stopped at 150 of R or 1000 in all. The range adds 1178
lines and removes 36 over 14 files; without the plan and the ledger, 568 added: R 125, tinytest 404 (386 in the
new file), help 3, NEWS 8, design notes 15, TODO 13. Inside both stops.

Left: `code-matrix-category-count` and `na-level-predictor-doors` in TODO. A new level is refused by name in the one
function `refuseNewDataLabels`; accepting one needs a ruling on the side it takes in each rule and engine work,
since a categorical rule's direction mask has no bit for it and the rule prior depends on the level count
(dec-A213). The plan's lines on a state moved between samplers whose levels differ stay as dec-A210 has them.
