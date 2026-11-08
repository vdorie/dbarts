# state-basis-levels: a stored state records each factor forest's levels and installs by name

Status: ACTIVE, 2026-11-07. Ruling: dec-B647, which revises dec-A210's call 6. No call is open.

agent: opus engine-builder, one; one opus review that runs the mutants
([Process by risk](README.md#process-by-risk): the bridge's state install changes, and the draws after
an install in another level order change).
rng: stated per call sequence.
- NEUTRAL, bit for bit: every sampler without a forest multiplied by a factor, whose states carry no new
  record and install through the same code; and an install of a state in the receiving sampler's own
  level order (a sampler's own state, a copy, a reload, a sampler made from the same factor), whose
  coefficients go in unpermuted.
- POSTERIOR-NEUTRAL, as a correction: a state installed in a sampler whose factor has the same levels in
  another order. The coefficients go onto their own levels where today they are exchanged, so the
  chain's first draws after the install change; the chain's stationary distribution does not.
- A state installed where the levels' names differ is refused where today it installs by position.
window: pre-release, before 1.0-0.
budget: about 250 lines: R about 70, bridge about 40, tests about 110, help and plan the rest.

## Goal

A state stored from a sampler with a forest multiplied by a factor records that forest's level names,
and `$setState`, `copy()` and a reload put each level's coefficient on the receiving sampler's level of
the same name. A state whose recorded names differ from the receiver's is refused, the sampler left as
it was. Help and the refusal say so.

## Context

- A factor forest. The forest's coefficients are each chain's amplitude block for that forest, one entry
  per level, in the order of the forest's indicator columns, which is the order of
  `data@basis.levels[[f]]` ([`dataBasisLevels`](../../R/data.R)). The block rides each chain's `glue`
  entry of the state, `[K, q_1..q_K, amplitudes, K prior variances]`, written by the bridge's
  [`storeState`](../../src/R_interface_bartcore.cpp) and read by
  [`readAmplitudeGlue`](../../src/R_interface_bartcore.cpp); the engine installs a forest's block only
  where it draws it ([`AmplitudeForestCombiner::restoreGlue`](../../src/bartcore/combiner.hpp)). The
  leaves are scalars and the kept draws in a state are trees only: a run's per-draw coefficients are the
  fit's `glue` channel, never the state's, so there is no kept coefficient in a state to permute. The
  prior variance is one per forest and does not depend on the order.
- The levels are fixed for a sampler's life: `$setForestBasis` matches a replacement by name and keeps
  the creation's levels, and `$setData` and `$setModel` refuse a sampler that carries amplitudes. So a
  sampler's own state, its copy and its reload always meet the same levels; only a state handed to
  another sampler can meet another order.
- Records already on a state: the bridge writes `cutPoints`, `cutMass`, `missing.columns` and others as
  top-level attributes, read by name, an absent optional one defaulted (the registry rule at
  [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)). A stored tree is moved by value onto the
  receiving sampler's grid (dec-B459); this does the same for a coefficient, by name.
- Today: [`installStateOnto`](../../R/dbarts.R) hands the state to the bridge, which installs by
  position. The help ([`dbartsSampler$setState`](../../man/dbartsSampler-class.Rd), its Saving section
  and `setForestBasis`) says "A state carries no levels, so ... goes in by position".
- Probe on the base build (shipped), 150 rows, a true difference of 0.92 between the arms, the state
  taken after 300 draws, 50 trees: the donor's last difference 1.03; the first draw after an install in
  the same order 0.92, in the reversed order 0.03.

## Calls made in planning

1. Where the record lives: an R attribute `basisLevels` on the state, written R-side from
   `data@basis.levels` by one helper every store goes through, beside the bridge's `cutPoints`. The
   alternatives: the engine's state record, which would put R's level names into an engine that keeps no
   names (the engine stays caller-neutral); or the bridge writing it, which needs the data object at
   every `storeState` call, which the bridge does not hold.
2. The record is written only where some forest multiplies a factor, so the state of every other sampler
   is the object it is today. The alternative, always writing it, changes every stored state's
   attributes for nothing.
3. Who permutes: R matches the names and hands the bridge one integer order per forest (NULL where
   nothing moves); the bridge, which owns the glue's layout, reorders the block of each chain before the
   engine sees it. The alternative, R rewriting `state[[c]]$glue` itself, ties R to the wire layout.
4. A state without the record installs by position, as today: a state of a sampler with no factor
   forest, or one stored by this development build before the record. Compatibility is owed only to
   shipped versions and 0.9-34 had no factor bases, so refusing such a state was the alternative; it
   would refuse the reload of a development-build fit for no gain, its levels being its own.
5. The match: for each forest, both sides multiply no factor, or both multiply a factor with the same
   set of names, which is then matched by name. Anything else is refused by name: different names, a
   different count, a factor forest against a numeric one or one with none, and a forest the record
   names that the receiver does not have. A held forest (`amplitude = fixed()`) is matched the same way,
   though its block is not installed: one rule for every forest. The alternative skipped held forests.
6. forceUpdate: a level mismatch is refused under either value. forceUpdate repairs trees to the grid;
   there is no repair of a coefficient that has no level to go to. The alternative, a forced install by
   position, brings back the silent exchange the ruling removes.
7. The refusal is an error, as a state of another number of forests or another leaf model is, and not
   `FALSE`, which says a tree would have to change. `$setState` checks the levels before it rebuilds a
   sampler whose engine is gone, so a refusal leaves even such a sampler as it was. The words:
   "the state's forest 2 multiplies a factor with the levels 'control', 'treated', and this sampler's
   forest 2 one with the levels 'a', 'b': a state's coefficients go only onto levels of the same names,
   in any order", with "no factor" for a side that has none.
8. `copy()` of a live sampler stores the engine's chain through the same helper, so the copy's `state`
   carries the record; a copy of a dead sampler and a reload install the stored state, which carries
   the record it was stored with. The `state` field after `$setState(st)` is `st` with its own record,
   so a later copy or reload matches it by name again.

## NEWS

None: factor bases are new in 1.0-0, and 0.9-34 is not affected.

## Constraints

- Draws do not move for any sampler without a factor forest or any install in the same order: the
  equivalence trio bitwise with strict coverage on the reference build, every scenario identical; the
  four reproducibility snapshot files unchanged.
- No engine change; the bridge's `bartcore_setState` takes one more argument.
- Out of scope: `installTrees`, which refuses a sampler of several forests.
- Owed later: a pin of the permutation's direction on a factor of three or more levels, when a forest
  may multiply one. On two levels a permutation is its own inverse, so an inverted order passes every
  test today.

## Steps

1. Plan (this file).
2. R: the store helper on every `storeState` call; the level match in `installStateOnto` and before
   `$setState`'s rebuild; the bridge argument and reorder; the help sentences in `setState`,
   `setForestBasis` and the Saving section, the refusal added to the list of states refused.
3. Tests, in `inst/tinytest/test-forest-kind.R` where the positional pin is: the record's shape; the
   probe (each level's coefficients after install, the first draw near the donor's in either order); a
   names mismatch refused with the sampler unchanged under both forceUpdate values; a copy's and a
   reload's state carrying the record and installing by name; a state without the record by position.

## Verification

- Both libraries installed with `--preclean` (shipped and `--enable-reference-build`).
- tests/cpp plain and under `OPT="-O2 -g -fsanitize=address,undefined"`.
- The full tinytest suite in one process, 0 failed.
- The trio `--bitwise --strict-coverage` on the reference build against the baselines MANIFEST names.
- The four reproducibility snapshot files on the reference build; the exact gates in `quick`.
- lintr, air, the codoc, win-drift and doc-freshness checks, and `mutation-battery.R verify-anchors`.
- Mutants: the permutation skipped in the bridge; the record dropped from a live copy's state. Each
  fails a pin.
