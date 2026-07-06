# data-ownership

agent: opus
rng: neutral
window: pre-release (user-visible aliasing semantics freeze at 1.0-0)
budget: ~300 lines

## Goal

ColumnStore owns its raw predictor storage (copied at ingestion)
instead of borrowing REAL(x). The bridge's protection-slot pinning for
predictors and the const_cast write-through into R vectors disappear.

## Decision

Question: keep classic's borrow-and-alias semantics (per-observation
updates visible through the originating R matrix,
R/bartcore.R:487-488) or own copies? Recommendation: own. The design
doc already allows it ("allowed to own",
docs/design/core-generalization.md cold layer), the mutation paths
already copy for rollback, the data-handle views already own, and the
pinning apparatus (src/R_interface_bartcore.cpp:35-50, PROT_* slots)
plus in-place writers (src/bartcore/data.hpp:778,794;
src/bartcore/sampler.hpp:869) exist only to service borrowing. Evidence
that would change it: a measured ingestion-copy cost that matters, or a
consumer relying on the aliasing (bairrtt should be checked). VD signs
off before implementation.

## Constraints

- R-side contract: setPredictor and friends already swap data@x on
  success R-side; that stays. What changes: writes are no longer
  visible through a matrix the caller kept a reference to.
- The weights/offset/response borrows can move in the same pass or stay;
  plan covers predictors and test predictors only - smallest coherent
  change. Note which PROT_ slots remain.
- Out of scope: sparse column storage (already owned); dbarts.h changes
  (the header documents borrowing - update the comment only).

## Steps

1. ColumnStore: own x/x_test (std::vector<double>), fill at build;
   remove the const_cast writers in favor of writes to owned storage.
2. Bridge: drop PROT_PREDICTORS/PROT_TEST_PREDICTORS retains and the
   post-success retain dance around setPredictor; simplify rollback.
3. R/bartcore.R: delete the aliasing note; confirm data@x swap-on-
   success semantics unchanged.
4. Update dbarts.h's creation-contract comment ("borrows for the
   sampler's lifetime" no longer true for predictors).
5. Component test: sampler results identical when the source R matrix
   is modified or GC'd after creation.

## Verification

- Component tests + full tinytest; equivalence exact (neutral).
- gctorture smoke run of the mutation surface
  (test-sampler-setPredictor*.R under gctorture(TRUE), abbreviated).
- bench-sampler compare: creation cost may rise by one copy; sampling
  metrics unchanged.
