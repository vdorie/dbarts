# test-data-parity

agent: opus (planner drafts the full plan against the code when the
  item is picked; this file records the directive and its seams).
rng: expected neutral (storage and ingestion parity; prediction values
  unchanged) - the picked plan pins this per step.
budget: TBD at pick time.

## Goal

Test data objects are handled the same way as training data: stored as
a frame / columnar source rather than a forced dense matrix, with
per-column kinds spanning {dense, sparse} x {ordinal, categorical} end
to end - ingestion, the engine's test codes, prediction, and the
setTestPredictor mutation surface. After it lands, the training and
test sides share one data model (VD directive, 2026-07-13).

## Context (seams the full plan must ground in code)

- Training-side precedent: the owned container and frame ingestion
  (docs/plans/data-ownership-2-ingestion.md), reference-install
  mutation (-3), column views (-4), CSC-categorical kernel (-5).
- Test side today: x.test is a dense matrix; the engine holds dense
  testCodes (data.hpp); sparse-categorical test columns densify at
  ingestion (data-ownership-5-sparse.md decision 2 and Q4 - this item
  SUPERSEDES that densification interim).
- Subsumes the "sparse x.test" bullet of docs/plans/
  sparse-extensions.md (the other three extensions stay there).
- Surfaces: dbartsData test slots, setTestPredictor / test-offset
  paths, xbart's fold views (buildFromParent test rows), predict.

## Constraints (standing)

- dbarts.h frozen; internal .Call surface only.
- Existing dense-matrix test input remains accepted and byte-identical
  (equivalence + tinytest gates per the usual classes).
- Sequencing open: queued by directive, need established; the pick
  decides where it lands relative to the rest of the backlog.
