# data-ownership

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

1. Converge the design note's "Still open" list with VD (in progress
   2026-07-11).
2. Freeze the design note; split into the five implementation plans
   (container/dense ingestion/borrow fast path; data.frame-direct
   ingestion + wrapper; flags + mutation rewire + snapshot semantics;
   views/sharing + standalone handle; sparse categorical).
3. Implement sequentially after cran-readiness lands (one implementer
   at a time).
