# cut-points-undo: a set cut grid can be put back, and leaves nothing behind

Status: PLANNED (dec-B231).

agent: opus implementer, one (the store and the bridge); opus reviewer.
rng: POSTERIOR-CHANGING on two call sequences: a grid longer than `n.cuts` is installed on a column (by
`setCutPoints`, or by `setState` with a state stored on such a grid) and a later `setData` derives that
column's grid. Today it derives at the longer grid's length, afterwards at `n.cuts`; the set of split points
is part of the tree prior. NEUTRAL for every other fit and call: draws are bit for bit unchanged. Proved by
the bitwise gates below (no baseline, exact gate or snapshot file calls `setCutPoints`) and by a seeded
digest of fits that never set a grid, taken on the base and the slice builds.
window: pre-release; any time. Serial with any other work in the store's cut builders or in
[`bartcore_setCutPoints`](../../src/R_interface_bartcore.cpp). The quantile rule's own change
([quantile-grid-spread.md](quantile-grid-spread.md)) has landed, so nothing is waiting on that side.
budget: ~290 lines (C++ store ~25, bridge ~20, R ~10, tests/cpp ~75, tinytest ~120, manual, store reference
and TODO ~40), upper figure 450. Plans have run 1.5-2x low; the design this comes from estimated 180.

## Goal

Code that runs the sampler inside a larger sampler can propose a change of cut points with `setCutPoints`
and, on rejection, put the old grid back with a second `setCutPoints` and restore the stored state, on any
design: one with a constant numeric column, and one with factor columns. After the grid is back, nothing
remembers that it was ever changed: a later `setData` derives the grid `n.cuts` names, as an untouched
sampler, a copy and a reload do.

## Context

All numbers below were run on the tip's build; a column is numeric unless said otherwise.

- The hidden count. The store keeps, per column, a cap on the cut count, a field this slice removed
  (retired: [`maxNumCuts`](../../src/bartcore/data.hpp)). It starts as the count asked for (`n.cuts`).
  [`ColumnStore::setCutPointsForColumn`](../../src/bartcore/data.hpp) raises it to a longer grid's length
  and nothing lowers it. [`ColumnStore::setData`](../../src/bartcore/data.hpp) derives through
  [`ColumnStore::buildCutsForColumn`](../../src/bartcore/data.hpp), which reads the cap, under the uniform
  rule and the quantile rule alike. The cap is not `data@n.cuts`, no state carries it, and a sampler made
  again from its R objects starts it at `n.cuts`.
- Measured, `n.cuts = 20`, two columns, 200 rows, under each rule. Set 50 points on column 1, then set the
  original 20 back: the grid is identical to the original and `data@n.cuts` still reads 20. A plain
  `setData` then derives 50 points on that column where an untouched twin derives 20, and the next five
  draws differ from the twin's by up to 0.058 (uniform) and 0.158 (quantile). The sampler's own `copy()`
  and a reload derive 20. The same residue is left by `setState` with a state stored on a 50-point grid
  followed by `setCutPoints` back to 20: `setState` installs the state's grid through the same store call.
- The two counts that are not defects, pinned as they are. A grid shorter than `n.cuts` (3 points) is kept
  at 3 by `setPredictor(updateCutPoints = TRUE)` and replaced by 20 at `setData`. A 50-point grid is kept at
  50 by `setPredictor(updateCutPoints = TRUE)` (accepted, both rules). After this plan a `setData` on the
  50-point column derives 20 where today it derives 50: that is the change of the `rng:` line.
- Grids the sampler holds and `setCutPoints` refuses. Under the uniform rule a constant column's grid is its
  one value repeated (100 points, 1 distinct), and a column whose spread is below what doubles resolve
  repeats too (100 points, 5 distinct). [`bartcore_setCutPoints`](../../src/R_interface_bartcore.cpp) holds
  a caller's grid to the strict form of [`cutGridIsValid`](../../src/bartcore/data.hpp), so handing a
  sampler its own grid for such a column fails with `$setCutPoints: requires strictly increasing cut points,
  none of them NaN`. The state reader holds a stored grid to the non-strict form, because the store builds
  these grids itself. Under the quantile rule a constant column gets one point and is accepted.
- The whole list. With `column` missing, [`bartcoreSamplerSetCutPoints`](../../R/bartcore.R) names every
  column, and the bridge stops at the first factor: `cannot set cut points for a categorical predictor`, or
  `cannot set cut points for an ordered factor predictor`. The grid a sampler reports has an entry per
  column (empty for an unordered factor, the level midpoints for an ordered one), so the list read from a
  sampler cannot be handed back on any design with a factor. Naming only the numeric columns works.
- What already works. On a design with neither case, store, `setCutPoints(G1)`, `setCutPoints(G0)`,
  `setState(stored)` returns `TRUE`, the grid is back bit for bit, and the continuation differs from an
  untouched twin's by 6.7e-16, the size a plain store and restore already leaves.
- The engine already runs on a grid with equal neighbours: a state whose column grid repeats every point
  (200 points, 100 distinct) installs, 200 sweeps are finite, splits fall on grid values, and a copy is made.
  The door that refuses such a grid is `setCutPoints` alone.
- The grid in force is read in tests as the `cutPoints` attribute of the stored state, as
  ["cutPointsOf"](../../inst/tinytest/test-quantile-grid.R) and
  [test-sampler-degenerate-cuts.R](../../inst/tinytest/test-sampler-degenerate-cuts.R) do.
- Existing pins that move:
  ["strictly increasing"](../../inst/tinytest/test-bartcore.R) expects `setCutPoints(c(0.5, 0.5), 1L)` to
  be refused. Existing pins that stay: a named factor column is refused
  (["categorical predictor"](../../inst/tinytest/test-bartcore.R),
  ["cannot set cut points for an ordered factor predictor"](../../inst/tinytest/test-data-categorical.R));
  an empty grid, a grid past 65533 points and a `NaN` are refused
  (["at least one cut point"](../../inst/tinytest/test-sampler-degenerate-cuts.R)).

## The rule

- The store remembers the count asked for, per column, for its whole life. Every derivation (creation,
  `setData`) uses it. No cap is kept beside it: a refresh re-cuts at the count the column holds, so a set
  grid leaves nothing that a copy and a reload do not also have.
- `setCutPoints` takes a strictly increasing grid of one to 65533 points, none `NaN`, and one grid more:
  the grid the column holds at the call, bit for bit (a -0 for a 0 is another grid), repeated points
  included. Given the whole list, the entries of factor columns are not read, whatever they are.

## Constraints

- A sampler that never has a grid longer than `n.cuts` installed draws exactly what it draws now.
- A refresh keeps the count the column has (`setPredictor(updateCutPoints = TRUE)`), as now, at 3 and at 50.
- An ordered factor's count is still its level count less one
  ([`ColumnStore::fillCutsAtLevelMidpoints`](../../src/bartcore/data.hpp)); `n.cuts` does not reach it.
- A refused [`Sampler::setState`](../../src/bartcore/sampler.hpp) leaves every later refresh of a column
  served as before it: it puts the grid and its count back, and the asked count never moved.
- Naming a factor column stays refused with today's two messages. A whole list of the wrong length stays
  refused: `$setCutPoints: requires one cut point vector per column`.
- No change to the flat C header (it has no entry for cut points), to the state format, or to `data@n.cuts`.
- No argument is added and no default changes.
- Out of scope: recording the grid in force on the data object, `setData`'s arguments, and `setState` no
  longer installing a stored grid. They are the later work under TODO `state-frame-prior`; the manual's
  step-by-step recipe for undoing a grid change is written there, when the data object can report the grid.

## Steps

1. Store. [`ColumnStore`](../../src/bartcore/data.hpp) gains the asked count per column
   ([`requestedNumCuts`](../../src/bartcore/data.hpp)), filled by
   [`ColumnStore::build`](../../src/bartcore/data.hpp) and copied by
   [`ColumnStore::buildFromParent`](../../src/bartcore/data.hpp) (the column-subset arm included), and
   loses the cap. The quantile collectors take their bound as an argument:
   [`ColumnStore::buildCutsForColumn`](../../src/bartcore/data.hpp) gives the asked count, a refresh and
   its feasibility check the count the column holds. tests/cpp, in
   [`testRequestedCutCount`](../../tests/cpp/test_data.cpp): set longer, shorter, then the original back,
   and the store is as one never changed; set longer, keep it, `setData`: the derived count is the asked
   one (fails today: the longer length), under each rule; a refresh onto a column with enough distinct
   values keeps a set count above the asked one and below it, each rule, on a dense column and on a
   CSC-backed one, whose refresh and feasibility check count for themselves; a view of a parent whose columns
   ask for different counts carries each through the column map
   ([`testColumnStoreView`](../../tests/cpp/test_data.cpp)); an ordered factor's count is what it is today.
2. Bridge and R. [`bartcore_setCutPoints`](../../src/R_interface_bartcore.cpp) holds a grid to the strict
   form of [`cutGridIsValid`](../../src/bartcore/data.hpp) unless it is bit for bit the grid the column
   holds and, when the column argument is `NULL`, takes one entry per predictor and skips those of factor
   columns; with nothing left to install it returns having changed nothing.
   [`bartcoreSamplerSetCutPoints`](../../R/bartcore.R) passes `NULL` when `column` is missing instead of
   naming every column and takes a data frame as the list of its columns. A list of another length goes to
   the bridge, which refuses it, with no entry read. Of a list of the right length it drops the entries of
   factor columns unread and refuses any other entry that is not numeric, as
   `$setCutPoints: 'cuts' must be numeric`: a function, `NULL`, and what `as.double` would have taken at the
   tip, a character vector, a logical, a factor (as its codes) and a Date (as its day count), none of which
   a test pinned. The rest it coerces to double. The bridge keeps the same refusal for another caller of
   the entry; no call through R reaches it. The refusal for an unusable grid becomes
   `$setCutPoints: 'cuts' must be strictly increasing and not contain NaN, unless it is the grid the column
   holds`.
3. tinytest, a new file `test-cut-points-undo.R`; "fails today" names what the tip does.
   - The undo leaves no residue: `n.cuts = 20`; sampler and twin run and store; the sampler sets 50 points,
     sets the stored grid back and restores (`TRUE`); the twin restores its own state; both take the same
     plain `setData`. Grids identical and the next five draws identical. Fails today: 50 points against 20.
     Under each rule.
   - A sampler, its `copy()` and its reload derive alike after the same sequence. Fails today: 50, 20, 20.
   - A 50-point grid left in place, then `setData`: `n.cuts` points, `data@n.cuts` unchanged. Fails today.
   - A state stored on a 50-point grid installed, the grid set back, `setData`: 20 points. Fails today.
   - Pins of today's counts: a 3-point grid through a refresh (3) and through `setData` (20); a 50-point
     grid through a refresh (50, accepted, each rule).
   - A constant column and a narrow column, uniform rule: the column's own grid is accepted and the grid
     after is identical; the whole own list too. Fails today: refused.
   - An unordered factor and an ordered factor beside numeric columns: the whole own list is accepted and
     changes nothing; a whole list with the factor's levels and a function in the factor entries is accepted
     in silence and sets the numeric columns, and so does a data frame with a column per predictor; a bad
     numeric entry of a whole list is refused; naming the factor column is refused with each of the two
     messages; a list shorter or longer than the design is refused with no warning. The first two fail
     today. On a design of factor columns alone the whole list returns and the stored state is unchanged.
   - Any grid with equal neighbours that is not the one the column holds is refused, on a column whose own
     grid repeats a point too, as a decreasing grid, a `NaN` and an `NA` are, with the new message, by
     column and as an entry of the whole list, and the grid is left as it was; over a column of zeros the
     held zeros are accepted and as many negative zeros refused; a strictly increasing grid is accepted; a
     function, a character vector, a logical, a factor, a Date and `NULL` are refused by name, and whole
     numbers taken; a column named three times has each entry read. After a constant column's grid is
     changed its old grid is refused and `setState` brings it back.
   - A state on a shorter grid installed over a 50-point grid, the 50 points set again, then a state
     refused: the column still refreshes at 50, each rule.
   - A sparse column with 50 points set, then 3: each refresh returns `TRUE` and keeps the count, each rule.
   - The whole undo on a design with a constant column and a factor: store, set another grid on a numeric
     column, hand the stored state's whole list back, restore (`TRUE`), and the continuation is identical to
     a twin's that restored its own state. Fails today: refused at the second call.
   ["strictly increasing"](../../inst/tinytest/test-bartcore.R) stays a refusal of `c(0.5, 0.5)`, with the
   new text.
4. Mutations (Verification): apply, install with `--preclean`, run, report the failing counts, revert, `touch`.
5. Records. [`dbartsSampler$setCutPoints`](../../man/dbartsSampler-class.Rd): the `cuts` item (strictly
   increasing, the held grid excepted, the whole list does not read factor entries, a named factor column
   is refused) and one sentence that a later `setData` derives at most `n.cuts` points whatever grid was
   set; the method's docstring in R/dbarts.R to match. [data-store.md](../design/data-store.md): the
   `numCuts` and `cutPoints` items and the new count in place of the cap, with the measured defect in two
   lines; that file is this change's design record. TODO: rewrite the `state-frame-prior` entry (below),
   and add `repeated-cut-restore`.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library (the store gains a field); `tests/cpp` builds and
  passes, clean under ASan and UBSan; the full tinytest suite; the new file under ASan on the R-loaded path.
- On a reference build: the four `test-reproducibility-*.R` files pass unchanged, and the three compares are
  bitwise, every scenario reporting identical draws, counted per scenario with no `max |z|` line: 55 against
  `equivalence-1b7d730c.rds`, 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. Nothing is re-recorded: the two scenarios that call `setData` never
  set a grid, so their cap is the asked count before and after.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick` mode, unchanged. None calls
  `setCutPoints` or `setData`; no new exact gate, the changed sequence having no posterior of its own to
  state (it is the posterior of the grid `n.cuts` names, which those gates cover).
- One script on the base and slice builds digesting seeded draws of fits that set no grid, a `setData`
  under each rule, a refresh, a copy and a reload included: equal.
- Mutations, each expected to fail the named test. The plan's first, a cap that keeps its old raise, and
  the reviewer's, a refused `setState` that does not put the cap back, have nothing to change once no cap
  is stored; while one was, the first failed tests/cpp alone, a cap that is too high being invisible from R.
  - a derivation counts from the count the column holds, not the asked one, on both arms and on each alone:
    tests/cpp "set longer, keep it, `setData`" and tinytest "a 50-point grid left in place";
  - a refresh counts from the asked count: tests/cpp "a refresh keeps the set count" and tinytest's
    50-point refresh pin and its refresh after a refused state, under the quantile rule (refused); the
    same on the two CSC arms: tests/cpp "a CSC refresh keeps the set count" and tinytest's sparse refresh;
  - a view does not copy the asked count, on each arm, and the column-subset arm reads it by its own index:
    tests/cpp's view checks;
  - the bridge refuses the held grid too: tinytest "a constant column and a narrow column"; the bridge
    takes any non-decreasing grid: tinytest's refusals and the pin in test-bartcore.R;
  - the whole list does not skip factor columns, in the bridge and in R: tinytest "an unordered factor and
    an ordered factor";
  - a whole list's numeric entries are not checked: tinytest's whole-list refusals;
  - a whole list over factor columns alone raises: tinytest's design of factor columns alone;
  - a list longer than the design is taken: tinytest's wrong-length refusals;
  - the held grid is matched by value, not bit for bit: tinytest's negative zeros over a column of zeros;
  - R drops entries by position when columns are named: tinytest's column named three times.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status; `R CMD check
  --as-cran` on a tarball from a clean copy (R/ and man/ change).
- Not a hot-path change: one more integer per column, read at a grid install and at a derivation.

## Out of scope, and where it goes

TODO `state-frame-prior`, restated (the entry now says "Not designed"; the coordinator has the text). What
this plan leaves to it: the grid in force written to the data object by every call that changes it,
`setCutPoints` included; `setData(updateCutPoints = )`; a state install that compares the stored grid and
never installs it; the manual's table of what puts each derived value back before `setState`.

## Calls made in planning

- A caller's grid stays strictly increasing, and the one grid taken with equal neighbours is the grid the
  column holds, bit for bit (decided after the first review). The alternative, this plan's first call,
  took any non-decreasing grid. A stored split names its cut by value and a restore puts it on the first
  index holding that value, so on a grid with every point tripled a store and restore moved the next 30
  draws by 0.569, and beside missing values it installed trees the prior gives probability zero; the undo
  needs none of that. Cost: once the grid of a column whose own grid repeats a point has been changed,
  `setCutPoints` does not take the old one back, and `setState` brings it. The restore is left alone: the
  same move happens on every grid that repeats a value, and the store's two rules, a refresh that keeps a
  constant column's grid and a state install each produce one, which is TODO `repeated-cut-restore`.
- The bridge skips the factor entries of a whole list by the store's own column kinds, and R drops them
  unread before it coerces the rest, by `data@varTypes`, the record the bridge builds those kinds from
  (after the first review). Before it R coerced every entry, so a factor's levels in its own place warned
  and a function failed with a message naming neither the method nor `cuts`.
- No cap is stored (after the first review); this plan first kept it as a field beside the asked count.
  A cap that is too high cannot be seen from R; one that is too low refuses every later quantile refresh
  of the column, and a stored one is left low by any refused install that does not put it back. Cost: the
  tests/cpp checks that read the cap by name read the asked count instead.
- After a 50-point grid is set and kept, `setData` derives `n.cuts` points, not 50: `n.cuts` is the one
  documented count, and the other reading keeps a number no R object shows. Released 0.9-34 had the same
  cap (read from its source, not run), so that one sequence changes against it. No NEWS entry: the manual
  never said which count a later derivation used, and the sequence is rare. The alternative is one line
  under the user-visible changes; the coordinator may prefer it.
- The manual does not yet print the undo recipe for a grid change: the documented way to read the grid in
  force arrives with the data object's record. The tests read the stored state's attribute, as three test
  files already do.
- The tip against the design, whose measurements were re-run and stand. The design ordered this work in
  series with the quantile rule's change, which has since landed. It names `setCutPoints` as what raises the
  cap; `setState` raises it too (measured), and the rule covers both. It did not list the pin in
  test-bartcore.R that step 3 rewrites. Nothing in this slice was found done already.
