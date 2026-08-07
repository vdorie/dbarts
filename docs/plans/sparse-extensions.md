# sparse-extensions

agent: opus
rng: neutral per extension (sparse landed with dense-identical gates)
budget: per extension; consumer-gated

## Goal

The deferred sparse-column extensions land when a consumer needs them,
each preserving the landed invariant that sparse builds match dense
builds of the same values exactly.

## Context

- docs/design/sparse-columns.md landing notes defer: in-place mutation
  on sparse columns (raw-x surface fixed at creation today), sparse
  x.test (DELIVERED by docs/plans/test-data-parity.md - resident sparse
  test storage; removed from this scope), a streaming range kernel, and
  mutation on dense-backed mixed columns. The per-column u8 width entry
  moved to hot-layer-u8.
- rbart_vi and linear leaves refuse sparse inputs; lifting those
  refusals rides whichever extension motivates it.

## Constraints

- Consumer-gated: none of these start without a named workload.
- Each extension keeps the dense-equality component-test device from
  the landing (dgCMatrix vs dense bitwise identity).
- Out of scope: order-preserving partitions (measured and rejected in
  the design doc; re-open only with new numbers).

## Steps

1. On demand, pick the extension, write its short plan as a section
   appended to sparse-columns.md (the design doc already carries the
   analysis), implement under the standing process.

## Verification

- The landing gates from sparse-columns.md, per extension.

## Status

- Extension (i), sparse-column in-place mutation, LANDED 2026-07-22
  (sparse-columns.md "Extension (i)" section): setPredictor/updatePredictor
  mutate CSC-backed columns, in-place when the nonzero pattern holds and
  rebuilding the layout otherwise, with dense-identical codes/draws and
  transactional roll-back. Of the extensions deferred here, two have since
  been delivered elsewhere (2026-08-07 correction): sparse x.test by
  docs/plans/test-data-parity.md (typed sparse test store, CLOSED there),
  and dense-backed mixed-column mutation, R mirror included, by
  typed-ingestion slice 2a (95db5d7). Still deferred and consumer-gated:
  a streaming range kernel, rbart_vi/linear-leaf on sparse, per-column
  u8 widths.
