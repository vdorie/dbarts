# data-ownership

Status: COMPLETE, 2026-07-14. All five plans landed; plan-5's landing note
records "the data-ownership program (plans 1-5) is COMPLETE"
(docs/plans/data-ownership-5-sparse.md). The design note
(docs/design/data-ownership.md) froze 2026-07-11 - that FROZEN refers to
the design converging, not to the implementation being unfinished; the
program itself is done.

agent: fable (design); implementation plans follow from the design note
rng: neutral (design phase)
budget: design note + fresh implementation plans

## Goal

Predictor storage is owned by the engine as a typed, quantized,
XGBoost-DMatrix-style container ingested directly from a data.frame;
the borrow-and-alias apparatus (PROT_* pinning, const_cast
write-through) disappears with the borrow. The full running record -
decision history, panel synthesis, open considerations - lives in
docs/design/data-ownership.md (draft).

## Status

- Original copy-raw plan rejected (VD, 2026-07-06); owned-quantized
  direction approved, design proceeding.
- Three-lens opus panel rendered 2026-07-06; synthesis in the design
  note: never-retain default, explicit creation-time column flags
  (re-cuttable; mutable), read-only borrow fast path kept 2/3.
- APPROVED FOR PRE-RELEASE IMPLEMENTATION (VD, 2026-07-11): "we'll do
  it before release" - 1.0-0 is NOT imminent, so the program runs at
  full scope. VD's 2026-07-11 design review sharpened the never-retain
  reading: the R layer keeps the ingested frame/snapshot and supplies
  raw values AT CALL TIME for re-cut/setData, so engine-side raw
  retention narrows to mutable-flagged columns only (the CoW-free
  hot-path updates that cannot round-trip through R per iteration);
  the re-cuttable flag may dissolve into call-time supply. To be
  converged with the remaining Still-open items before the note
  freezes.

## Next

1. DONE 2026-07-11: design converged and FROZEN (resolutions in the
   design note: call-time raw supply, no epoch, construction-only
   borrow, @x drop attempt in plan 2, sparseFactor exported,
   extract(sampler, "predictors")).
2. Write and implement the five plans sequentially (one implementer
   at a time): (1) container + dense ingestion + construction-only
   borrow (engine core); (2) data.frame-direct ingestion +
   sparseFactor + the @x drop attempt (R surface); (3) mutable flag +
   mutation-surface rewire + call-time supply plumbing; (4)
   views/sharing + standalone handle (blocks forest-split-bcf); (5)
   sparse categorical kernel + test-side sparse.
