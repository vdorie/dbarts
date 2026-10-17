# factor-levels-across-setdata: setData reads a factor by label against the sampler's levels

Status: PLANNED 2026-10-17, written to the tree of f514c1d8 (origin/bartcore). Not critiqued. Closes TODO
`factor-levels-across-setdata`. Rulings it stands on: dec-A115 (declared levels), dec-A85 and dec-B155
(a level matched by label, an unseen one refused), dec-B279, dec-B286 and dec-B359 (every door codes a
factor column by label and refuses numbers for it), dec-B346 (a level with no rows for a while stays),
dec-B231 (setData changes derived data only unambiguously), dec-B195, dec-B418 and dec-A210 (a state
carries no data), dec-B449 and dec-B450 (a state installs no grid), dec-B230; it settles what dec-A209
scoped out. One call is open for the maintainer ([For the maintainer](#for-the-maintainer)).

agent: tier "Changes draws" ([Process by risk](README.md#process-by-risk)): the diff is R only, but four
call sequences fit another model after it. A blind critique of this plan; one opus implementer for the R
code and tests (a slip here reads one level as another with no message, the defect being fixed); one opus
review that runs the mutants under [Tests](#tests) and checks the fix round. Nothing under src/ changes,
so no new numerics are reachable and no sanitizer run is owed.
rng: NEUTRAL, bit for bit, for every `setData` whose factor columns carry the sampler's own level tables,
and for every fit that never calls `setData`. POSTERIOR-CHANGING for four sequences, none in a benchmark,
a snapshot, a help example or a consumer: (1) `setData` handed a factor whose shared levels sit at other
positions (a level dropped before the last, or the levels reordered): its rows now carry their own
level's code, where today they carry the code of the level at their position; (2) a `copy()` or a reload
after `setData` handed a factor declaring other levels: it now has the sampler's level count and
continues as the original does, where today it is refused or runs on another count; (3) a factor whose
labels the sampler lacks, on a row (renamed levels): refused, where today it is read by position; (4) an
ordered factor whose shared levels come in another order: refused, today read by position.
window: R surface, file-disjoint from the engine slices; after forest-kind-by-class push 2 or before it,
rebase only (both touch [R/dbarts.R](../../R/dbarts.R) and the sampler Rd). Before TODO
`state-frame-prior`'s record slices, which it leaves less to do.
budget: planned ~470 lines (R ~70, tests ~330, help, NEWS and notes ~70); upper figure 900. Stops:
[Budget and stops](#budget-and-stops).

## Goal

`$setData(newData)` reads each factor column of `newData` (training and test) by label against the level
table the sampler was created with, as `setPredictor`, `setTestPredictor` and `predict` already read one,
and stores the data object re-coded onto that table. A factor declaring fewer levels, the same levels in
another order, or more levels that no row has, goes in with each row at its own level; a row at a label
the sampler lacks is refused by name; a level with no row now stays one the trees may split on. The
engine's level table and the record a re-creation reads are then always the same, so a copy and a reload
continue as the original does and every install of the original's stored state goes in.

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
  2026-10-08: "matrices should only work when matrices were given as an input"; the rule: "setPredictor
  with no column takes a data frame, coded by label as the column form is; a numeric matrix is taken only
  where every predictor column is numeric ... and refused by name on a sampler with a factor column". It
  names setTestPredictor and predict as not ruled; neither ruling names `setData`.
- dec-B346, the maintainer on 2026-10-08: "Take both on a swap"; a factor level no row holds is taken on
  a swap, "the forest's part on it drawn under its prior until a row returns".
- dec-B231, the maintainer on 2026-10-05: "You don't want to update derived data in an ambiguous case. It
  needs to be explicit, because a caller might not understand that setData is changing a different part
  of the sampler." Its rule re-derives at a plain `setData` the values derived from the rows (split
  points, standardization, response range); a level table is not derived from the rows (below).
- dec-B195, the maintainer on 2026-10-02: "To the degree possible, state != model", and "state (and
  scratch) describe the current literal state of an MCMC sampler as it walks the posterior, as defined by
  a model and data." dec-A210 reads dec-B231 and dec-B418 as "a state carrying no data" and has a state
  between samplers whose levels differ in order go by position. dec-B449: "a live setState never writes
  the sampler's cut points"; TODO `state-frame-prior` moves the grid off the state.
- dec-B230: a state from other data is "neither supported nor prevented".
- dec-A209 scoped this out: "whether setData gives the live sampler the new table or a re-creation takes
  the table from the state (about 30 engine lines, a draw change for such copies)".

Code as it stands:

- A data object made from a data frame keeps each factor as one column of 0-based codes into its own
  level table, `attr(x, "factor.levels")` ([`makeCategoricalModelMatrix`](../../R/utility.R)).
- A creation reads each factor's level count from that table
  ([`readDeclaredCategoryCounts`](../../src/R_interface_bartcore.cpp)); an ordered factor's grid is its
  count less one. The count is fixed for the engine's life.
- [`bartcore_setData`](../../src/R_interface_bartcore.cpp) refuses a code at or past the engine's count
  and installs the codes as given. The R side ([`bartcoreSamplerSetData`](../../R/bartcore.R)) stores
  the handed data object as the record. So the live engine reads the new rows by position in the new
  table against its old count, and the record declares the new table.
- A copy, a reload and a dead pointer's `setState` re-create from the record
  ([`recreatePointer`](../../R/dbarts.R)), so they take the new table's count, then install the stored
  state, whose trees hold category masks and ordered thresholds by code.
- `predict` re-codes a data frame against the record's table by label
  ([`mapFactorColumnsToTrainingLevels`](../../R/utility.R)); `setPredictor` by column codes labels the
  same way ([`codeCategoricalColumnUpdate`](../../R/bartcore.R)), and a whole-frame `setPredictor` too
  ([`codePredictorFrame`](../../R/bartcore.R)). `setData` is the one door that does not.
- [`heldAfterAdoption`](../../R/dbarts.R) refreshes a re-created sampler's `state` field where its
  ordered column's grid differs from the recorded one; the `setData` sequence below is its only
  documented reach.
- No consumer calls `setData`: stan4bart (branch bartcore) re-creates its samplers from data it built
  and calls `setState`; bartCause (dbarts-1.0) and treatSens call neither; bairrtt calls `setPredictor`
  by column. Nothing in [dbarts.h](../../inst/include/dbarts/dbarts.h) swaps data. Outside the tests,
  dbarts calls `setData` only in the composition matrix, on subsets of a sampler's own data object.

### Measured: the tip beside 0.9-34

Scripts in scratch/fls of this worktree (p2.R, sequences; p3.R to p7.R, corners), run on a build of
f514c1d8 (scratch/libs/fls) and on 0.9-34 (scratch/libs/cran). 200 rows, a numeric column `a` held
fixed across data sets and a factor `f` of levels a to f, the response `a` plus 0, 2, 4, 0, 2, 4 by level,
every data set on the first one's response range; 20 trees, 30 sweeps of burn-in and 10 kept draws, one
chain; then `setData` with the levels in the first column. "Draws as the original" means one sweep from
the same R seed gives the original's training fits to rounding (at most 1.3e-15, the fits re-summed) or
identically. "By position" means the kept draws' prediction at each label after the swap is the one they
gave the label at that position before it (differences 3.9 to 4.1 against the same label, 0 against the
position); where the two readings agree the row says so. `setState(stored)` hands back the state
`setData` stored.

| Levels handed over | 0.9-34 (both kinds, with or without NA) | tip, unordered (with or without NA) | tip, ordered | tip, ordered with NA on it |
|---|---|---|---|---|
| a..f, the same | taken; copy and reload draw as the original; `setState` NULL everywhere | taken; copy and reload draw as the original; `setState` TRUE everywhere | the same | the same |
| a..d, last two dropped | refused: "maximum number of cuts must be of length equal to 5" (a factor is indicator columns) | taken (label = position); copy and reload refused, "state is not consistent with this sampler"; live `setState` TRUE | taken; copy and reload made on 3 points, draws differ by up to 0.135; their `setState(stored)` refused except a reload's first call | as unordered: copy and reload refused |
| a, c..f, b dropped | refused (length 6) | taken, by position; copy and reload refused | taken, by position; copy and reload made, draws differ by up to 0.98; `setState` as above | taken, by position; copy and reload refused |
| a..h, g and h on no row | taken (0.9-34 drops a level no row has); copy and reload draw as the original | taken; copy and reload made on 8 levels, draws differ by up to 0.27 (0.57 with NA) | taken; copy and reload made; their `setState(stored)` refused except a reload's first call | the same |
| a..h, rows at g and h | refused (length 9) | refused: "categorical predictor values must be existing category codes" | refused: "ordered factor predictor values must be existing level codes" | the same |
| f..a, reordered | taken, by position; copy and reload as the original | taken, by position; copy and reload as the original | the same | the same |
| A..F, renamed | taken, by position (label = position) | taken, by position; copy and reload as the original | the same | the same |

On every tip row: `predict` reads a new frame by label against the new table, so it agrees with the live
sampler's positional reading, and refuses a label the new table lacks ("test data factor 'f' has levels
not present in the training data") although the engine still holds that level; a reload predicts as the
live sampler. After the two "dropped" rows a later `setData` back to the six levels, with rows at all six,
is taken by the live sampler and refused by its copy (the copy has fewer levels, or cannot be made).

An example of "by position" (p7.R): with b dropped, the kept draws predict 2.43 at level c after the swap
(the value they learned for b; 4.58 before it). The live chain re-learns: after 100 further sweeps the
mean fit by level is within 0.07 of the truth, so the harm inside a larger sampler is a chain pushed off
its posterior at every swap whose table moves, and kept draws that predict the wrong level.

Corners on the tip:

- A data object made from a matrix of codes declared categorical (`attr(x, "varTypes")`), handed to a
  sampler made from a data frame: taken, the record then has no level table, and `copy()` is refused
  (p3.R). A sampler made from such a matrix, handed codes 0 to 3 of its 0 to 5: `copy()` refused, the
  re-creation inferring the count from the codes present (p4.R).
- The test set of a data object with other levels is coded against that object's own training table, so
  it reads as its training rows do (p3.R).
- A predictor factor with an NA level is created; `predict` at a row of that level is refused as a level
  training never saw (p5.R): the by-label helper matches `as.character`, which loses an NA level.

Emulating the recommended rule on the tip (p2.R `recode`: the new factor re-coded by label onto the
sampler's table before `dbartsData`), on all three kinds: every taken row draws as the original in copy
and reload, `setState(stored)` is TRUE on the live sampler, the copy and a reload as first and later call,
the kept draws predict by label, `predict` at a level the new rows lack is taken (2.58 at e, the value
learned for e), and the rows at g and h, and the renamed levels, are refused.

### What base R does

Measured with lm (scratch/fls/lm.R). `predict.lm` reads a new factor by label against the fit's levels:
fewer levels, another order (ordered factors included) and an extra declared level no row has all
predict as the full factor would; a new label on a row is refused ("factor f has new levels e"), renamed
levels too. `update(fit, data = )` refits from scratch and takes the new data's levels: that is a new
model, the analogue of a new sampler, not of a swap.

## What a user means

`setData` gives the same model new observations: inside a larger Gibbs sampler that imputes a predictor or
moves rows, once per sweep. A factor there is a variable whose possible values were fixed when the model
was stated; a row's value is its label. The level table a data object carries is not a statement the
caller makes about the model: every factor carries one, and it shrinks or reorders whenever a factor is
rebuilt from characters (`factor(x)` on a draw that missed a rare level drops it; `factor(x, levels =
unique(x))` reorders). A user handing such data means "these rows, at these levels", never "level c is now
level b" or "the variable now has four values". So the correct model keeps the sampler's levels, reads
each row by label, keeps a level no row has as a value the trees may still split on (dec-A115, dec-B346),
and refuses a label the model never had, as `predict` does. A different set of levels, or another order of
an ordered factor, is another model, which is a new sampler.

dec-B231's rule re-derives at a plain `setData` the values derived from the rows, and uses values a data
object carries. A level table is neither: it is not computed from the rows, and taking the one every
factor carries as a deliberate record would make every `setData` change the sampler's levels, the
ambiguity the ruling forbids. Where dec-B231 moves a numeric split onto a re-derived grid by value, so
that a threshold keeps its meaning, the factor's analogue is reading by label, which this plan does by
re-coding the rows and not the trees. The critique should test this reading.

## Options

**A. setData gives the live sampler the new data's table, rows by position (the TODO's first).** As a
model: no, while rows are read by position; a level dropped before the last moves every later one.
Draws: the live sampler's change on every swap whose count differs (the tree prior counts levels, and
rules on a vanished level are repaired). Copy, reload and state: agree with the live sampler. Engine
~200 lines: a per-column count through `SamplerBase::setData` (internal virtual, `--preclean`), the
store's counts, codes and an ordered grid resized, a mask's word count and the missing-value code moved
where the count crosses 64, live trees with rules on a vanished code repaired through the forced-update
path, kept draws' masks trimmed. Bridge and R ~30. dbarts.h untouched. Consumers: none calls `setData`.
Inside a larger sampler, a level whose rows vanish for one sweep loses its rules, and comes back with
none when its rows return.

**A2. The same, with the trees re-mapped by label.** As a model: yes, for a variable whose values really
changed. Draws: as A, and rules move with their labels. Engine ~300 to 450 lines beyond R's map: A's, plus
every live and kept category mask permuted, ordered thresholds carried through a map that must be
increasing (another order refused), and a rule for the side a level new to the table takes, which needs
a ruling. Inside a larger sampler it destroys a rare level's rules whenever a sweep's rows miss it.

**B. A re-creation takes the level count from the state (the TODO's second).** As a model: the live
sampler is left as it is, so rows stay read by position. Draws: copies and reloads after such a swap now
continue as the original. Engine ~30 lines (an adopting install sets a factor column's count from a
record), bridge and R ~40 (an unordered factor's count recorded in the state, which records nothing for
it today; an ordered one's is its recorded grid). It adds data to the state against dec-B195, dec-B449
and TODO `state-frame-prior`, and the state's format changes: stan4bart stores states in its fits, so a
fit saved before must still restore (an absent record read as today).

**C. setData reads a factor by label against the sampler's table (recommended).** As a model: yes, the
one above. Draws: the four sequences under `rng`; bit for bit elsewhere, the tables being equal. Copy,
reload and state: the record's table is the engine's, so all agree with the live sampler and the stored
state installs everywhere. Engine 0 lines; bridge 0; R ~70. dbarts.h untouched. Consumers: none; NEWS
states the change from 0.9-34 (renamed levels now refused, other level counts now taken).

**D. setData refuses a factor whose level table differs from the sampler's.** As a model: safe, no
wrong reading reachable. Draws: none move; sequences that run today are refused. Copy, reload and state:
never reach a mismatch. R ~15 lines. The cost is the larger sampler's: a loop that rebuilds a factor each
sweep stops the first time a rare level is not drawn, and the caller must write `factor(x, levels =
<the sampler's levels>)` on every call, which is what C does for them; `setPredictor` would take what
`setData` refuses.

## Recommendation

C. It is what a user means and what base R does with a factor's levels; it gives `setData` the rule every
other door already has (dec-A85, dec-B155, dec-B279, dec-B359), keeps dec-A115's declared levels and
dec-B346's empty level through a swap, and needs no data in the state (dec-B195, dec-A210). It is the
cheapest correct option: no engine, bridge or ABI change, and the TODO's draw change for copies comes with
it. A and B leave a level read as another; A2 is correct for a different variable but costs an engine
slice and a ruling, and harms the larger sampler; D is safe but refuses what the other doors take.

## The rule

At `$setData(newData)`, for each predictor column the sampler records a level table for (`T`, from
`attr(sampler$data@x, "factor.levels")`), in the training and the test predictors of `newData`:

1. `newData`'s column must be a factor column with its own table `N`. A numeric matrix of codes, a column
   with no table, is refused by name in the words `setPredictor` uses for a matrix (dec-B359): "give a data
   frame, whose factor columns are coded by label".
2. Where `N` is identical to `T` nothing is touched: the codes and the object go in as handed.
3. Otherwise each level of `N` is matched to `T` by label, an NA level to an NA level of `T`. A level of
   `N` that `T` lacks and that some row has is refused, naming the column and the labels, in the words of
   `setPredictor` by column ("column 'f' has label 'g' not among its training levels"); one no row has is
   dropped. A missing value stays missing.
4. For an ordered factor, the levels of `N` that rows have must be in `T`'s order; otherwise refused
   ("lists its levels in another order than the sampler's; another order is another model, made with a
   new sampler").
5. The column's codes become each row's level's position in `T`, and its table `T`: the data object the
   sampler stores carries the sampler's tables.

Nothing is installed or stored before every column passes. A column the sampler records no table for (a
sampler made from a matrix of codes) is left as today ([Out of scope, and where it goes](#out-of-scope-and-where-it-goes)).

## Constraints

- `bartcore_setData`, the engine and dbarts.h have no diff; the bridge's code check stays the backstop.
- An identical table is the code that runs today: the record is the handed object, bit for bit.
- A refusal leaves the sampler, its record and its `state` field as they were.
- The handed object is the caller's: R's copy semantics leave it untouched; only the stored record is
  re-coded.
- `heldAfterAdoption` stays; with C it is reached only by a record assigned by hand or an object saved by
  an earlier build of this branch ([Calls made in planning](#calls-made-in-planning)).
- No new argument, no warning. The by-label match and its refusal reuse the words of
  `codeCategoricalColumnUpdate`; the column match follows [`alignContainerFactorLevels`](../../R/utility.R)'s
  dense arm, but matches levels and not `as.character` of rows, so an NA level keeps its name.

## Steps

1. The re-coding (~70 lines, [R/bartcore.R](../../R/bartcore.R), or [R/utility.R](../../R/utility.R)
   beside its siblings). A function (working name `alignDataFactorLevels(newData, x.train)`) applying
   [The rule](#the-rule) to `newData@x` and `newData@x.test` (a dense matrix or the dense-backed mixed
   container; `setData` refuses sparse ones before) and returning the re-coded object. Called in
   [`bartcoreSamplerSetData`](../../R/bartcore.R) after `checkBasisRecord` and before the `.Call`, so
   every refusal precedes the swap; the record stored is its result.
2. Tests ([Tests](#tests); ~330 lines). A new `inst/tinytest/test-setdata-factor-levels.R`. Rewrite, not
   delete: in [test-state-grid-install.R](../../inst/tinytest/test-state-grid-install.R) the block "an
   ordered factor after a whole-data replacement declaring fewer levels" (about 105 lines) pins today's
   rule (a copy on 3 points, the original's state refused). Under C that sequence gives a copy on 5 points
   that takes the original's state; the block's refresh pins move to the same fixture with the record then
   assigned by hand (`sampler$data <- <the data handed to setData>`, `n.cuts` and `sigma` copied), which
   recreates the mismatch C closes. In [test-bartcore.R](../../inst/tinytest/test-bartcore.R) the refusal
   of a row at level 11 expects "existing category codes"; it becomes the by-label text.
3. Help, NEWS and notes (~70 lines). The `setData` docstring in [R/dbarts.R](../../R/dbarts.R): one
   sentence on the rule and its two refusals. [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd),
   Saving: the sentences on `setData` to an ordered factor declaring fewer levels, and on the original's
   state refused by its copy, go; in their place, a copy and a reload hold the sampler's levels.
   [dbartsData.Rd](../../man/dbartsData.Rd): a swap reads the factor by label. NEWS: the entry "New
   predictors for a fit with a factor column come as a data frame, coded by label" gains `setData`, with
   the change from 0.9-34 (other level counts taken, reordered levels read by label, renamed refused).
   [install-surface.md](../design/install-surface.md): the paragraph on a factor column's level table
   across `setData` is rewritten to the rule; [data-store.md](../design/data-store.md) one sentence.
   The ledger entry for the calls made, at landing. TODO: drop `factor-levels-across-setdata`; add the
   two lines under [Out of scope, and where it goes](#out-of-scope-and-where-it-goes).

## Tests

tinytest only; tests/cpp has nothing to pin. Each row: what it pins; the mutant it must kill. Fixtures:
p2.R's data sets, 20 trees, kept draws; "draws as the original" is `expect_equal` of five sweeps' fits
from one seed (the fits are re-summed at a re-creation, 9e-16) and `expect_identical` of the trees.

| Test | Pins | Mutant |
|---|---|---|
| same table neutral | the record is identical to the handed object; next draws identical to a twin's given the object directly | the object re-coded where the tables are identical (non-identical attributes) |
| fewer, tail | unordered and ordered: taken; the record's table is the sampler's six; copy and reload made and draw as the original; `setState(stored)` TRUE on the live sampler, the copy, a reload's first and a later call | the record keeps the handed table (copy refused, unordered) |
| fewer, middle | b dropped: each row's code is its label's position in the sampler's table; kept draws predict at c what they predicted at c before | codes left as handed (by position) |
| reordered, unordered | f..a: by label as above; copy as the original | tables compared by length, not identity |
| reordered, ordered | refused naming the column; sampler, record and field untouched | the order check dropped |
| ordered, a level dropped | b dropped from an ordered factor, the others in order: taken | the order check run on the declared levels rather than the shared ones |
| more, no rows | g and h declared, no row: taken; record table six; copy as the original | an unused unknown label refused |
| unknown on a row | rows at g and h, and renamed A..F: refused naming the labels; sampler untouched | the refusal dropped (unknown labels become missing values) |
| test set | a test set declaring c and e only: test fits as by label | x.test not re-coded |
| back again | after the tail case, `setData` of all six levels with rows at each: taken by the live sampler and by its copy | (pins the closed TODO sequence) |
| NA level | a sampler made with an NA level, handed the levels reordered: the NA level's rows keep it, a missing value stays missing | match through `as.character` (NA level rows become missing) |
| missing values | the ordered NA fixture of p2.R, tail case: copy and reload made and draw as the original | (pins the TODO's missing-value reach) |
| code matrix refused | a sampler made from a frame, handed a categorical code matrix: refused in dec-B359's words, sampler untouched | the matrix check dropped |
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
| 1 R | ~70 |
| 2 tests | ~330 |
| 3 help, NEWS, notes, TODO | ~70 |

Estimates here have run 1.5 to 2 times low. Stop and report, building nothing further: if a bitwise
compare prints a `max |z|` line or a snapshot file moves; if any change under src/ or dbarts.h appears
needed; if a tinytest file outside step 2's three fails (list it, do not edit it); if stan4bart's suite
moves; at 150 R lines or 900 in all.

## Interactions

- forest-kind-by-class push 2: [R/dbarts.R](../../R/dbarts.R), the sampler Rd and NEWS; rebase only. Its
  `checkBasisRecord` call in `bartcoreSamplerSetData` comes first, as now.
- TODO `state-frame-prior`: its record of the grid on the data object is the numeric columns' business;
  after this slice a factor column's grid already follows the record, and the record the sampler's table.
- TODO `warm-start-salvage` and the warm start: a donor's state goes in by position, dec-A210 unchanged.

## Out of scope, and where it goes

- A sampler made from a matrix of codes records no table, and a re-creation infers its count from the
  codes present, so `copy()` is refused after `setData` of fewer codes (p4.R). New TODO line
  `code-matrix-category-count`: record the count such a sampler was made with.
- `predict` and `setPredictor` match a factor through `as.character`, so a predictor's NA level is
  refused at `predict` (p5.R), and at `setPredictor` a missing value would read as that level. New TODO
  line `na-level-predictor-doors`, beside dec-A210's refusal of an NA level for a basis.
- A state moved between samplers whose levels differ: by position (dec-A210), unchanged.

## Calls made in planning

- An ordered factor whose shared levels come in another order is refused. The alternative reads it by
  label in the sampler's order, as `predict` reads a frame; the cost taken is one refusal, the reason that
  an ordered factor's order is how it splits, so another order is another model.
- A declared level no row has and the sampler lacks is dropped, not refused, as `predict.lm` takes one.
  The alternative refuses it, at the cost of a loop whose factor declares a superset of levels.
- A numeric matrix of codes is refused on a sampler with a recorded table, extending dec-B359 to
  `setData`. The alternative carries the sampler's table onto the codes, as dec-A166 did for
  `setPredictor`, which dec-B359 then replaced with the refusal.
- An NA level is matched by name. The alternative refuses a row at an NA level with dec-A210's remedy;
  the cost taken is none at this door, creation accepting such a factor.
- `heldAfterAdoption` is kept and its pins move to a mismatch made by hand. The alternative removes it
  (about 10 lines and five pins), reached then by nothing documented; left for the critique.
- The record stored is the re-coded object, not the handed one, where the tables differ. The alternative,
  the handed object, is what makes a copy differ today.

## Open calls

One, for the maintainer: which option ([For the maintainer](#for-the-maintainer)). The plan is written
for C.

## For the maintainer

A sampler can be handed new data in the middle of a run, which a larger model does every sweep, for
example when it imputes a predictor. When a predictor is a factor, the new data's factor carries its own
list of levels, and that list need not be the one the sampler was made with: a level no row has may have
been dropped, the levels may come in another order, or a new one may appear. Today the running sampler
keeps its original list but reads each new row by its position in the new list, not by its name. So if
level "b" is dropped, every later level moves up one and rows at "c" are read as "b": a stored draw that
predicted 4.58 at "c" before the swap predicts 2.43, what it learned for "b", after it. Separately, a copy
or a saved and reloaded sampler is rebuilt from the new list, so it either cannot be made or runs a
different chain from the original. 0.9-34 refused a factor with another number of levels and also read
reordered levels by position.

The question is which list governs when new data brings a factor with other levels.

1. The new data's list replaces the sampler's, rows still read by position. Copies then match the
   original, but "c" is still read as "b", and the running sampler must cut tree rules that used a level
   now gone, which changes its draws. About 200 lines of engine code.
2. The same, with the tree rules moved by name. Right for a variable that really changed, but each time a
   rare level is missing for one sweep its rules are destroyed, and it needs a new ruling on which way a
   newly added level goes. About 300 to 450 lines of engine code.
3. A copy or a reload takes the level count from the saved state. About 30 lines of engine code and a
   change to what a saved state holds; copies then match the original, but "c" is still read as "b", and
   it puts data into the saved state, against the rulings that a state carries none.
4. Refuse a factor whose list differs. About 15 lines and no draw changes, but a loop that rebuilds its
   factor each sweep stops the first time a rare level is not drawn, and the user must restate the
   sampler's levels on every call, while the other ways of handing in new predictors accept it.
5. Read each new row by name against the sampler's own list, as replacing one predictor column and
   predict already do. Fewer levels, another order and extra levels no row has are all taken; a level the
   sampler never had, on a row, is refused by name; a level with no rows stays one the trees can use, as
   at creation. Copies and reloads then match the original and its saved state installs everywhere. About
   70 lines of R, no engine change, nothing for the consumer packages. The costs: draws change where a
   level used to be read as another; renamed levels, which 0.9-34 read by position, are refused; and an
   ordered factor with its levels in another order is refused, since that order is the model.

Recommended: 5. It is what a user means, "b" is "b", and what base R's predict does with factor levels; it
gives this method the rule every other way in already has, and it costs no engine code.
