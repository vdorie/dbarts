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
  direction approved, design proceeding, nothing time-critical.
- Three-lens opus panel rendered 2026-07-06; synthesis in the design
  note: never-retain default, explicit creation-time column flags
  (re-cuttable; mutable), read-only borrow fast path kept 2/3.
- VD open considerations under discussion: updateX must stay
  copy-on-write-free (proposed: engine-owned raw for declared-mutable
  columns; data@x becomes a creation-time snapshot); multi-model data
  sharing (proposed: shared standalone handle + column-subset views).

## Next

1. Converge the design note's "Still open" list with VD.
2. Freeze the design note; split into implementation plans (container,
   ingestion/R surface, mutation surface, views/sharing).
