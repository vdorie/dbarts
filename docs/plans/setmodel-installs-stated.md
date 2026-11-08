# setmodel-installs-stated: $setModel installs each forest's stated number

agent: opus
rng: neutral for every call taken today with a model that states nothing
different; draw-changing for the calls ignored or refused today
budget: ~150 lines of R, ~250 of tests, help and NEWS sentences

## Goal

`$setModel(m)` installs every number a forest's record in `m` states that
differs from the sampler's: the variance forest's power and base, and a mean
forest's own sd, power and base, each through the code `$setVarianceForest`
or `$setForest` runs, so the result is bitwise that door's. A number the
record does not state, a copy from the fitting function's level among them,
is left as it is. Structure is still refused before anything is stored, and
the later write wins between `$setModel` and the named door (dec-B652).

## Context

- dec-B652 (the rule); dec-B524, revised in part; dec-B254, dec-B518,
  dec-B523, dec-B551, dec-B554; dec-A258's calls 11, 19 and 31, numbered in
  setmodel-parameters-not-structure.md, Calls made.
- `setModel` (R/dbarts.R); `setModelOnForests`, `levelDepths`,
  `levelLeafPrior`, `refuseShareDifference` (R/fittingLevel.R);
  `keepModelRecords`, `movedStructureRefusal` (R/records.R); the variance
  record in `resolveSamplerSpec` (R/spec.R); `writeForestDoor`,
  `writeVarianceForestDoor`, `modelWithDepth`, `modelWithForestSd`
  (R/forestDoor.R).
- A variance forest exists only on a model of one forest (creation refuses
  `variance` beside several forests), so its install is on the one-forest
  route alone.

## Where each statement is recorded

- Variance forest: `dbartsVarianceRecord@stated`, today `"n.trees"` only.
  Creation adds `"power"` and `"base"` where `varianceForest()` states them;
  an unstated one holds the mean depth copied at creation and is not stated.
  `$setVarianceForest` adds the depth names it writes.
- Mean forest on a model of several: `dbartsForestRecord@stated` names a
  stated power or base (creation and `$setForest` already write it); an sd
  is stated where the record's `sd` is not NA. Forest 1's depth is read from
  the tree prior (call 32), its record holding the same.
- A model of one forest and a multinomial model: the forest's own depth and
  size are the fitting function's (the tree prior and leaf prior), which
  `$setModel` installs today; nothing new is recorded.

## What $setModel installs, by sampler kind

- One forest: as today, plus the variance forest's stated power and base
  where they differ from the sampler's, by `$setVarianceForest`'s
  validation, record write and engine call (`setForestTreePrior` on the
  variance forest).
- Several forests, BCF included: per depth name, the level first, as today:
  it reaches the candidates that state none of their own on the sampler, and
  the one candidate whatever it states (dec-B523); its value is read from the
  reached forests the given model does not state. Then each forest whose
  given record states the name at a value other than the one it then holds
  takes it as `$setForest` writes it: the record, its `stated`, forest 1's
  tree prior, and the engine's per-forest depth. Then the leaf prior as
  today, and each forest whose given record states an sd other than its
  record's then takes it as `$setForest(f, sd = )` writes it
  (`modelWithForestSd`, `setForestSd`). A number the given record does not
  state is left, so a basis forest's copied depth and a forest stating its own
  on the sampler no longer refuse the call.
- Multinomial: the forests take one statement, read from the tree prior,
  which `$setModel` installs today; an sd is refused at creation. No change.
- Heteroscedastic: one forest; see the first item.

## Order

Every check, then every engine write. The checks: structure, the leaf scale,
split probabilities, coefficients (still refused, naming `$setForest`), the
parts of a shared total compared with those the sampler holds after the own
installs, the variance numbers validated as the door validates them, and an
sd install's unit read from the engine's row, refused in the bridge's words
where the response has no spread. The writes: the route's engine call (the
level, depths, residual prior, move settings), then the mean forests' sd
installs, then the variance forest's depth; each of the later two can no
longer refuse once checked. The model is stored last.

## Copy, reload, older models

A copy and a reload re-create from the model, whose records then hold the
installed numbers and statements, so they run what the sampler runs. A model
saved by this development build carries no variance depth statement: it
reads as not stated and is left, as today. 0.9-34 had no variance forest and
no per-forest settings, so nothing shipped is owed a migration.

## Constraints

- No draw moves for a call taken today with a model stating nothing
  different: the trio bitwise against benchmarks/baselines/MANIFEST, and
  the inert pins on every sampler kind.
- No engine or bridge change; no change to `$setForest`'s or
  `$setVarianceForest`'s draws.
- Out of scope: the coefficient's prior through `$setModel`, the
  multinomial per-category settings.

## Steps

1. The variance record's statement at creation and at its door.
2. The variance install in `setModel`.
3. The mean forests' own depth and sd installs in `setModelOnForests`.
4. Tests in test-setmodel-forests.R: the pins below; the hand-edit pins
   rewritten on models built through doors (dec-B554).
5. Help on setModel, setForest, setVarianceForest; the NEWS item on the
   forest door.

## Verification

- B built with `varianceForest(power = 4)` (and `base = 0.8`):
  `A$setModel(B$model)` draws bitwise as `A$setVarianceForest(power = 4)`.
- A mean forest's stated sd, power and base: `$setModel` of a model built
  with them equals `$setForest(f, ...)` bitwise on several forests.
- `dbartsSpec(data, power = 2.5, variance = varianceForest())$model` moves
  the mean depth and leaves the variance forest.
- A sampler's own model is inert on every kind.
- Another variance tree count is refused and the sampler is unchanged.
- A copy and a reload run what the sampler runs.
- A later `$setVarianceForest` after `$setModel` wins.
- The gates in the prompt: both libraries, tests/cpp, the suite, the trio
  `--bitwise --strict-coverage`, the four snapshots, the exact gates in
  quick, lint, air, codoc, win drift and doc freshness.

## Calls made

1. A variance depth statement lives in the existing `stated` slot, as a mean
   forest's does. Alternative: a slot of its own.
2. A stated number is installed where it differs from the value the forest
   holds once the level is written; an equal one writes and states nothing,
   as the door does. Alternative: record the statement whatever the value, as
   a new sampler from the model would.
3. A number the given record does not state is left, so the own-difference
   refusals on unstated numbers and the refusal of a withdrawn statement
   (call 35, for `$setModel`) go: dec-B652's "a number the record does not
   state is left as it is". Alternative: keep refusing a withdrawal, naming a
   new sampler.
4. The level first, the own installs after, as `$setTreePrior` then
   `$setForest` would; the level's value read from the reached forests the
   given model does not state. Alternative: own installs first, when a
   forest the given model states at the sampler's value would then follow
   the level.
5. The one candidate keeps dec-B523's rule, the level reaching it whatever
   it states; a number the given record states there is an own install,
   which records the statement. Alternative: the level's write alone, the
   flags kept (call 28).
6. A coefficient difference stays refused naming `$setForest`: dec-B652
   names sd, power and base. Alternative: install it through
   `setForestCoefficient`.
7. The parts of a shared total are compared after the own sd installs, so a
   model from a sampler whose forest stated an sd live is taken, and one
   created with that forest stating it, whose sharers took other parts, is
   refused as today. Alternative: none that keeps the parts set once.
8. A hand edit that changes a value and its statement together is read as a
   statement (dec-B554: a record slot is not a door, and the record cannot
   tell the two apart); one that changes a value alone is left.
