# snapshot-tests

agent: sonnet
rng: neutral engine-wise (test values regenerate, draws unchanged)
budget: ~400 lines moved/rewritten across 5 test files + 1 script

## Goal

The RNG-locked literals stop depending on whole-file execution history
and stop pretending to be correctness tests: per-family reproducibility
files, reseeded per block, compared against script-regenerated
reference objects, labeled as seeded-drift tripwires.

## Context

- ~65 hardcoded values across 5 files:
  test-continuousResponse-regression-singleThreaded.R (~37),
  -multithreaded.R (~11), test-binaryResponse-regression.R (~10),
  test-rbart-reproducibility.R (~5), test-xbart-reproducibility.R (~2).
- Regeneration today replays whole files ("draws depend on the
  preceding fits through R's stream position and last-ulp arithmetic",
  test-rbart-reproducibility.R:118-119).
- Their original oracle (BayesTree cross-package values) is dead; their
  remaining value is a CI-run tripwire for unintended seeded-result
  drift. Correctness lives in the exact-posterior gates and component
  tests, not here.

## Constraints

- Coverage must not shrink: every model family / feature combination
  currently pinned keeps a pinned block.
- Values move from hand-pasted literals to .rds references only if
  tinytest ergonomics stay good (diff readability on failure);
  otherwise keep literals but reseed per block - the per-block
  reseeding is the non-negotiable part.
- Out of scope: deleting statistical/behavioral tests in the same
  files; those stay where they are.

## Steps

1. Extract the pinned blocks into test-reproducibility-<family>.R
   files; each block does set.seed + fit + compare, independent of
   file position.
2. tools/regenerate-snapshots.R: re-runs each block and rewrites the
   reference values (literals via deparse or a common/*.rds), so
   regeneration is one script, not a file replay.
3. Header comment in each file: drift tripwire, not correctness; when
   an intentional RNG-shifting change lands, run the script and eyeball
   that shifts are plausible (magnitudes, not signs of bugs).
4. Delete the superseded pinned blocks from the original files.

## Verification

- Full tinytest passes before and after (values identical - this item
  must not itself shift anything: same seeds per block means NEW
  reference values; regenerate once, then assert stability by running
  the suite twice).
- Running tools/regenerate-snapshots.R twice in a row produces a zero
  diff.
