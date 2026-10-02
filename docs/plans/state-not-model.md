# state-not-model: a saved state holds the chain, not the model

Status: PLANNED 2026-10-02 under dec-B195, dec-B196, dec-B197 and dec-B200 in [decisions.md](../decisions.md).
Starts after [fit-stores-k.md](fit-stores-k.md) lands.

agent: opus implementer, one; opus reviewer.
rng: NEUTRAL. A state installed in a sampler under the model it was saved with gives the draws it gives today.
Only an install across a model change behaves differently, and no recorded baseline contains one.
window: pre-release, before the 1.0-0 merge.
budget: ~1400 lines (engine ~260, bridge ~110, R ~120, manual ~70, tests ~780, records ~60). Plans have run
1.5-2x low.

## Goal

A state holds what the chain is: trees, leaf values, the quantities the sampler draws, the generator, and the
units those numbers are stored in. It holds no prior parameter and no value the sampler holds fixed. A sampler
has one set of units, fixed by its model and data, and a state stored in other units is converted into them as
it is installed. Installing a state through `setState`, `copy`, a reload or a warm start never changes the
sampler's model, `getLeafPrior` reads the same before and after and never reports NA, and every chain of a
sampler runs under one prior.

## Context

- The writer and the two readers: [`storeState`](../../src/R_interface_bartcore.cpp),
  [`setState`](../../src/R_interface_bartcore.cpp), [`readWarmStartState`](../../src/R_interface_bartcore.cpp).
  The engine side: [`Chain::setState`](../../src/bartcore/chain.hpp),
  [`installForest`](../../src/bartcore/chain.hpp), [`Sampler::setState`](../../src/bartcore/sampler.hpp),
  [`installForests`](../../src/bartcore/sampler.hpp) and its undo,
  [`restoreTouched`](../../src/bartcore/sampler.hpp), which snapshots and re-installs a chain's state. The
  exchange structs: [`ForestStateData`](../../src/bartcore/combiner.hpp),
  [`ChainStateData`](../../src/bartcore/combiner.hpp). The R callers:
  [`setState`](../../R/dbarts.R), [`getPointer`](../../R/dbarts.R), [`copy`](../../R/dbarts.R),
  [`installTrees`](../../R/dbarts.R).
- What each block is, against the package's original division of data, model, state and derived scratch:

  | block | class | today | after |
  |---|---|---|---|
  | trees, leaf values, tree params, saved draws, latents, thresholds, generator, DART split weights and counter, drawn amplitudes, variance-forest trees, the saved-draw cursor | state | installed | unchanged |
  | `k`, `sigma`, `resid.df`, `shape`, `dart.alpha` | state where the sampler draws it, model where it holds it fixed or its family pins it | installed either way | written and installed only where drawn |
  | `leaf.scale` | model | installed; marks a calibration-map forest foreign | not written, not installed |
  | glue: amplitude prior variance | state on a scale-mixture forest, model on a fixed-variance one | installed on both | installed on a scale-mixture forest only |
  | glue: amplitudes of a forest created with `update.amplitude = FALSE` | model | installed | not installed |
  | `fit.scale` | the units the chain's numbers are stored in | installed, and a k-named prior's centre and width move with it | compared with the sampler's own; stored numbers converted when they differ |
  | `cutPoints`, `leaf.covariate.center` and `.scale`, a heuristic gp lengthscale | scratch, but frozen while the data moves, so only the state has them | installed | unchanged |
  | a supplied gp lengthscale | model | installed | the sampler's (dec-A146) |
  | weights and survival digests | data, by digest | compared | unchanged |

- No model value is a unit of a stored number. Leaf values, slopes and amplitudes are stored raw; the scale and k
  enter only the next draw's prior. Measured on 41 sampler configurations: an install with the `leaf.scale`
  block deleted gives the reader output and draws today's install gives, provided a named sd is re-stated
  afterwards where the transform is frozen away from the one the data now gives - the step
  [`reissueNamedLeafSd`](../../R/dbarts.R) already takes after every other channel that moves the transform.
- A drawn k is relative to the anchor, which is model. Installed under another anchor it is the stored k, not
  the stored spread: the chain starts one sweep at anchor over k. Same-model installs are unaffected.
- Precedents for keeping the sampler's own value, re-expressed in the installed transform:
  [`GaussianResponse::restoreScale`](../../src/bartcore/model.hpp) for the sigma prior and
  [`calibrateVarianceLeaf`](../../src/bartcore/chain.hpp) for the variance leaf.
- The registry rule at [`stateFormatVersion`](../../src/R_interface_bartcore.cpp): no format has shipped, so
  neither constant moves.
- The response transform and the prior. A k-named leaf prior's width is a constant times the range of the
  chain's transform and every leaf prior's centre is its midpoint; the calibration map's forest scales and the
  negative-binomial centre follow it too. Measured: installing a state saved on 3y + 10 moves a k-named
  sampler's anchor from 3.92 to 11.77 and its centre from 3.36 to 20.08, and a sampler given chains in two
  transforms runs them as two posteriors. No R field holds the transform: after a response swapped without the
  scale update, a sampler rebuilt from its own control, model and data derives another one, and only the state
  brings the saver's back.
- Each chain holds its own rescaled response, so the engine has no single transform today;
  [`GaussianResponse::restoreScale`](../../src/bartcore/model.hpp) moves a chain to an installed one and carries
  the sigma prior, inexactly: at an unchanged range the carry moves the prior by an ulp in about a third of
  cases, so one call more or fewer than today changes draws.
- Conversion, emulated in R on a constant leaf by rewriting a state before today's install: every leaf value,
  saved draws included, times the ratio of the ranges, plus the shift difference over the range spread across
  the trees. Replayed predictions agree with today's to 5e-16 relative and the converted chain then runs as a
  sampler wholly in those units does. By reading, a monotone leaf converts the same way, a linear leaf's slopes
  scale and its intercept shifts, the negative-binomial has a shift only, and variance-forest factors take the
  squared ratio spread over the trees as [`reanchorVarianceForest`](../../src/bartcore/chain.hpp) does to a live
  forest. A gp leaf's saved draw has no mean to shift, and with amplitudes a shift has no forest to own it.
- Consumers. stan4bart restores, for replay of kept trees, a state holding one chain from each of its per-chain
  samplers, into a sampler created from its specification's model; the chains' transforms differ from each
  other and from that sampler's, its sigma is held fixed and written during sampling, and the replay reads
  neither sigma nor the leaf prior. Its suite compares replays to 1e-10. bartCause and treatSens install no
  state. rbart's per-chain samplers come back from workers with a transform frozen away from the data's.

## Decision

Ruled: the leaf prior is the sampler's (dec-B195); a value held fixed is model and a drawn one is state, sigma
included (dec-B196); `setSigma` on a sampler that does not draw sigma rewrites the model's fixed value, recorded
on the R object as `setLeafPrior` records a spread, so `copy` and a reload keep it (dec-B197). A write through
the C header cannot reach the R object; a reloaded sampler then holds creation's value until the next write.
Applied under those rulings by the agents, recorded as dec-A146 for the maintainer's mark:

1. A warm start follows the same rule: the recipient keeps its model, and the donor's trees and the values the
   recipient draws seed it. Today `bart(leaf.prior = normal(k = 2), warm.start = fit)` with a fit made at k = 4
   runs at 4 while its model says 2.
2. A supplied gp lengthscale is the sampler's. A state holding saved draws made under another is refused, since
   a saved gp draw cannot be replayed under another kernel; without saved draws it installs.

Ruled after the first two were applied:

3. The response transform (dec-B200). The maintainer, asked what an installed transform does to the prior:
   "Well, prior != state, as we agreed."; and, offered conversion against each chain keeping its units with the
   prior re-expressed per chain: "Convert on install." The range a k-named prior is anchored to is a model
   value, recorded on the R object, set when the sampler is created or deliberately re-anchored and by nothing
   else. A re-anchor is therefore a model change that restoring a state does not undo: the restored values are
   converted into the re-anchored units.

## Constraints

- The invariant: for every state the suite and the consumers install today into a sampler under the saving
  model, the draws after the install are bitwise what they are now. Bitwise means restore against restore; a
  restore has never matched the uninterrupted chain to the last bit.
- Whether a value is installed is decided by the destination's own switch - does it draw k, sigma, the df, the
  shape, the concentration; is this forest a scale mixture; does it update its amplitudes - and then by
  presence: a destination that draws a value and finds no block keeps the value it has. No install is refused
  for a missing k, df, shape or concentration block.
- The predicates "draws its df" and "draws its shape" are new and additional.
  [`carriesResidualDf`](../../src/bartcore/model.hpp) and [`carriesShape`](../../src/bartcore/model.hpp) also
  mark the family for the latents and shift checks and keep that job.
- A block the reader no longer wants is ignored when present, so a state written before this change installs.
- A named sd is re-stated with creation's own arithmetic after the transform is installed, through the pointer
  the install used: [`reissueNamedLeafSd`](../../R/dbarts.R) fetches the pointer itself and would recurse inside
  [`getPointer`](../../R/dbarts.R).
- The undo of a failed warm start keeps working: what it snapshots and puts back must still return the
  recipient to where it was. With no model value installed, there is none to put back.
- The record is an attribute on the model object, the engine's exact pair; a slot would break samplers saved
  before it existed. `initialize` writes it from the engine on a first creation, whatever the model handed in
  carries, so a sampler made from another's model on other data anchors to its own data as today. The
  re-anchoring calls - `setResponse` and `setOffset` with the scale update, `setData` - refresh it; `setModel`
  carries it over; `setState`, `copy`, a reload and `installTrees` never write it. A re-anchor through the C
  header cannot reach it, the gap dec-B197 states for sigma.
- A sampler re-created from its own R object is in the record's units before anything reads it. Where an
  install follows the creation - a reload with a state, `copy` - the install's own move does this, exactly as
  today; where none follows, creation moves the chains itself. About ten entries read the transform before any
  run, the reader, `storeState`, predict and the prior draws among them, so the move is never deferred.
- Conversion is a pass over the incoming state ahead of the install, comparing each chain's `fit.scale` with
  the sampler's pair by exact equality. Equal: nothing is touched and today's install runs, its one
  `restoreScale` included, so the invariant above holds by construction. Different: the stored numbers are
  rewritten and the chain's `fit.scale` set to the sampler's, and the same install runs. The caller's state
  object is not modified.
- What cannot be converted is refused by name, with the sampler unchanged: a gp leaf or forests with amplitudes
  when the shift differs. Their scale converts.
- No install is refused because its chains carry different transforms from each other: stan4bart's restored
  samplers are such chains, and dec-B191's refusal is withdrawn. They are converted, each from its own units.
- A warm start converts the donor's values into the recipient's units and no longer adopts the donor's
  transform, so a donor on another range seeds the function it held.
- The reader takes the anchor, the prior mean and the response scale and shift from the sampler's one
  transform. Its branch for chains that disagree goes.
- The cut grid, the leaf standardization and a heuristic lengthscale install as today. They shape a prior too -
  the slope prior, the gp kernel, the split rule - and are left for a later item, named in the Landing note.
- No new entry or field in the shipped header; no change to either state-format constant.
- Out of scope: a family check in the state validity test (a gaussian state installs into a probit sampler
  today and still will, less its sigma); `copy()` reading the stored state rather than the live chain; recording
  `setCutPoints` or a header-set store capacity on the R object; moving the forest specifications off the
  control.

## Steps

1. Engine. Predicates for a drawn df and a drawn shape. The state writer emits k, sigma, the df, the shape and
   the concentration only where drawn and no leaf scale. [`Chain::setState`](../../src/bartcore/chain.hpp) and
   [`installForest`](../../src/bartcore/chain.hpp) install each only where the destination draws it and the
   block is present, an amplitude variance only on a scale-mixture forest (the hand-written BCF arm of
   [`restoreGlue`](../../src/bartcore/combiner.hpp) included), amplitudes only where the forest updates them,
   and never the leaf scale. The validity check requires a df or shape only of a sampler that draws one. Remove
   [`noteInstalledLeafScale`](../../src/bartcore/chain.hpp),
   [`adoptInstalledAmplitudePriors`](../../src/bartcore/chain.hpp),
   [`nodeScaleIsMapDerived_`](../../src/bartcore/chain.hpp) and their part of
   [`InstallMarks`](../../src/bartcore/chain.hpp). tests/cpp follows.
2. Bridge. `leaf.scale` leaves both parsers and the writer; `k`, `sigma`, `resid.df`, `shape` and `dart.alpha`
   become optional in both, `dart.alpha` no longer required beside `dart.probabilities`.
3. R. Re-state a named sd after the install in the three paths that install a state. Drop every NA branch of
   [`reportLeafPrior`](../../R/dbarts.R) and the NA-sd refusal in
   [`resolveForestSpreads`](../../R/dbarts.R).
4. `setSigma` on a sampler that does not draw sigma records the value on the model field.
5. The warm start through the same install rule; the lengthscale refusal, with its own flag beside those the
   install already reports.
6. The anchor. The engine holds the sampler's pair apart from each chain's transform; the bridge passes the
   record at creation and reads the pair back; the record attribute and its four writers in R; the creation
   route that moves the chains when no install follows, with `copy` routed as an install.
7. Conversion: the pass over the incoming state for each leaf model and for variance factors, saved draws
   included, ahead of [`Sampler::setState`](../../src/bartcore/sampler.hpp) and
   [`installForests`](../../src/bartcore/sampler.hpp); the two refusals.
8. Manual: `setState`, `storeState`, `copy`, `getLeafPrior`, `setLeafPrior`, `setSigma`, `installTrees` and
   `warm.start`, and the docstrings they mirror: what a state holds, that the model is the sampler's, that a
   state in other units is converted and to what accuracy, that a re-anchor is a model change a restore does
   not undo, with the sequence that rolls one back, and that the grid and standardization come with the state.
9. Tests. Rewrite the cases that pin a model value riding the state:
   ["FOREIGN CALIBRATION"](../../inst/tinytest/test-forest-basis-r5.R),
   ["a pre-write state"](../../inst/tinytest/test-multiforest-leaf-prior-writer.R),
   ["the per-forest leaf scale rides the state"](../../inst/tinytest/test-bcf.R),
   ["divergedScale"](../../inst/tinytest/test-calibration-midchain.R) and its neighbour "divergedK", the
   `leaf.scale` format cases and ["donor.sf"](../../inst/tinytest/test-sampler-state-format.R), the fixed-shape
   state oracle in inst/tinytest/test-shape-channel.R, the fixed-df and fixed-shape state checks and the
   leaf-scale case in tests/cpp, and the two parser anchors in benchmarks/R/mutation-battery.R. Add: for each
   of k, a named sd, a forest spread, sigma, the df, the shape, the concentration, a fixed amplitude variance
   and fixed amplitudes, store, write, restore, with the reader and the next draws matching a twin that only
   wrote; a state from one model installed under another leaves the recipient's reader unchanged; drawn values
   still install; a state with the old blocks present installs; a gaussian state leaves a probit sampler's
   sigma at 1. For the anchor and conversion: the reader identical across an install from a sampler on a
   rescaled response, for each naming and leaf model; prior-only draws at the sampler's anchor after it; the
   converted state's replayed predictions equal to the donor's to rounding and the live fit preserved, for
   constant, monotone and linear leaves, a variance forest, the negative-binomial and a gp leaf on an equal
   shift; the two refusals, sampler unchanged; chains from two samplers on different ranges combined into one
   state, installed, and run as one posterior; a re-creation after a swap without the scale update bitwise
   today's with its state, and under the saver's prior without one; a copy and a reload whose stored state
   predates a re-anchor; the rollback sequence; a warm start from a donor on another range.
10. Records: a design note carrying the block inventory and where the four-way division holds and does not, the
   state paragraph of docs/architecture.md, the ledger entries for the calls made here, the index row, the TODO
   item, and the Landing note naming every install path covered.

## Verification

Against a private library, installed with `--preclean` (step 1 changes virtuals):

- `cd tests/cpp && make && ./test_bartcore` passes.
- `tinytest::test_package("dbarts")` passes with no new warning. Mutation, run once and reported: putting the
  leaf-scale install back fails the new store, write, restore tests.
- The equivalence compare on a reference build reports identical draws for every scenario against the current
  baseline in `benchmarks/baselines/MANIFEST`, and the four seeded-drift snapshot files pass there.
- Every gate `.github/workflows/exact-gates.yaml` lists passes in `quick` mode: what a fit's stored state
  carries changes.
- stan4bart's suite, its store-trees test with several chains included, bartCause's and treatSens's pass
  against a private-library chain built on this tip.
- Mutation, run once and reported: with the conversion pass skipped, the combined-chains test fails.
- Sanitizers on tests/cpp and on the R-loaded path for the new state tests, as Gate hygiene describes: the
  conversion pass writes through every stored tree.
- `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`, `tools/check-win-drift.R` and
  `tools/check-doc-freshness.R` pass, each on its own exit status.
