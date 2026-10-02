# state-not-model: a saved state holds the chain, not the model

Status: PLANNED 2026-10-02 under dec-B195, dec-B196 and dec-B197 in [decisions.md](../decisions.md).
Step 3 waits on the transform design under Decision.

agent: opus implementer, one; opus reviewer.
rng: NEUTRAL. A state installed in a sampler under the model it was saved with gives the draws it gives today.
Only an install across a model change behaves differently, and no recorded baseline contains one.
window: pre-release, before the 1.0-0 merge.
budget: ~900 lines (engine ~150, bridge ~90, R ~50, manual ~50, tests ~500, records ~40). Plans have run 1.5-2x
low.

## Goal

A state holds what the chain is: trees, leaf values, the quantities the sampler draws, the generator, and the
frame those numbers are stored in. It holds no prior parameter and no value the sampler holds fixed. Installing
a state through `setState`, `copy` or a reload never changes the sampler's model, and `getLeafPrior` reads the
same before and after.

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
  | `fit.scale`, `cutPoints`, `leaf.covariate.center` and `.scale`, a heuristic gp lengthscale | scratch, but frozen while the data moves, so only the state has them | installed | unchanged (Decision 3) |
  | a supplied gp lengthscale | model | installed | Decision 2 |
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
- Consumers. stan4bart restores, for replay of kept trees, a state holding one chain from each of its per-chain
  samplers; the chains' response transforms differ, its sigma is held fixed and written during sampling, and
  the replay reads neither sigma nor the leaf prior. bartCause and treatSens install no state.

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

Open:

3. The frame: the response transform, cut points and leaf standardization stay in the state as the units the
   chain is stored in. The maintainer, asked about the transform: "Well, prior != state, as we agreed." A
   k-named leaf prior takes its centre and width from the transform, so an installed transform must not move
   it: the sampler's prior stays its own in response units and is re-expressed in the installed units, as a
   named sd and the sigma prior already are. The mechanism, and where the anchor is recorded so a re-creation
   reproduces it, are being designed; step 3 grows by it.

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
- No install is refused because its chains carry different response transforms: stan4bart's restored samplers
  are such chains, and dec-B191's refusal is withdrawn. What the reader reports for them follows Decision 3.
- The transform, the cut grid, the leaf standardization and a heuristic lengthscale install as today.
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
3. R. Re-state a named sd after the install in the three paths that install a state. Drop the reader's
   calibration NA branches in [`reportLeafPrior`](../../R/dbarts.R) and narrow the NA-sd refusal in
   [`resolveForestSpreads`](../../R/dbarts.R) to the disagreeing-chains case.
4. `setSigma` on a sampler that does not draw sigma records the value on the model field.
5. The warm start through the same install rule; the lengthscale refusal, with its own flag beside those the
   install already reports.
6. Manual: `setState`, `storeState`, `copy`, `getLeafPrior`, `setLeafPrior`, `setSigma`, `installTrees` and
   `warm.start`, and the docstrings they mirror: what a state holds, that the model is the sampler's, and that
   the transform, grid and standardization come with the state.
7. Tests. Rewrite the cases that pin a model value riding the state:
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
   sigma at 1.
8. Records: a design note carrying the block inventory and where the four-way division holds and does not, the
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
- `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`, `tools/check-win-drift.R` and
  `tools/check-doc-freshness.R` pass, each on its own exit status.
