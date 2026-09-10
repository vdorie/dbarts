# setpredictor-partition

agent: opus (S1 the engine change, its oracle, the re-record); sonnet
  (S2 the kernel bench width defect, file-disjoint, runs beside S1).
rng: shifting for S1. The root partition on the setPredictor path keeps
  the inherited member order instead of rewriting the identity, so leaf
  sums reassociate at the last bit on every fit that calls
  `setPredictor`; tree structure, varcount and the RNG stream are
  unmoved. Every equivalence scenario that never calls `setPredictor`
  must reproduce the current baselines BITWISE; those that do move at
  max |z| = 0.00 and are re-recorded under the MANIFEST's P17 rule with
  the oracle below. S2 is neutral.
window: engine, one function; pre-release (VD 2026-09-10, dec-B117).
budget: S1 ~40 engine + ~80 tests/cpp + the record; S2 ~10 C.

Decision: dec-B117 in [docs/decisions.md](../decisions.md).

## Goal

`dbartsSampler$setPredictor` updates run 19 to 37 percent faster with
no change to sampling speed or to which trees grow, because the
re-partition that follows a predictor swap no longer rewrites index
spans that are already in partitioned order.

## Context

- [`revalidateTrees`](../../src/bartcore/chain.hpp) re-partitions every
  tree from the root through
  [`repartitionSubtree`](../../src/bartcore/tree.hpp) after a
  predictor changes. The root partition it reaches uses the dense
  identity rewrite (`misc_partitionRange`) for a dense column, which
  writes the whole span; the span it receives is already partitioned
  under the live rules (measured: 7400 of 7400 root partitions on the
  accept path had no misplaced element), so the in-place kernel
  (`misc_partitionIndices`) does one comparison scan and no writes.
  The measurement, replicated on x86 and arm64, is recorded in the
  TODO under the engine constants audit item.
- During sampling the root's span is almost never already partitioned,
  and the identity rewrite wins there on sequential column reads, so
  the change is confined to the revalidation path; the sampling-time
  dispatch is untouched.
- Scenarios calling `setPredictor` exist in all three equivalence
  harnesses ([equivalence.R](../../benchmarks/R/equivalence.R),
  [bcf-equivalence.R](../../benchmarks/R/bcf-equivalence.R),
  [multinomial-equivalence.R](../../benchmarks/R/multinomial-equivalence.R))
  and in several tinytest files; snapshot-valued tests among them are
  replayed as whole files.

## Constraints

- Gates (shifting): tests/cpp; full tinytest with replayed snapshots;
  the three equivalence compares against the current baselines showing
  every non-setPredictor scenario "identical draws (same RNG stream)"
  and every setPredictor scenario under the statistical mode at |z| < 4;
  re-record all three with the P17 oracle named in the MANIFEST rows;
  bench-sampler compare on a quiet machine (the setPredictor rows are
  the claim; every run row must sit at 1.00 within the floor); local
  sanitizer legs on tests/cpp and the R-loaded path; `--preclean`.
- Oracle (P17), a cross-implementation check: a tests/cpp test that,
  after a predictor update, compares the member SET of every leaf of
  every tree under the in-place root against the identity-rewrite
  root over the same state, equal as sets, and the leaf sufficient
  statistics equal to within a stated summation tolerance; poisoned by
  a mutation that misplaces one member.
- Out of scope: the sampling-time root dispatch; the sparse layout
  (stays, by measurement); any change to which trees grow.

## Steps

S1:
1. In [`repartitionSubtree`](../../src/bartcore/tree.hpp)'s path only,
   route a dense root through the in-place kernel; the sampling-time
   dispatch keeps the identity rewrite. State the constraint in the
   comment: the span arriving here is already partitioned under the
   live rules, so the in-place scan writes nothing, and the identity
   rewrite would write the whole span.
2. The oracle test above, plus its poison run recorded in the landing
   note.
3. Replay snapshot tests that move; re-record the three baselines from
   the shipped build; MANIFEST rows naming the oracle and the partition
   (which scenarios moved, all at max |z| = 0.00, which reproduced
   bitwise); demote the predecessors; update the workflow pins and the
   feature matrix per the MANIFEST's re-record rule.
4. bench-sampler compare on a quiet machine; record the setPredictor
   rows and the run rows in the landing note.

S2:
5. benchmarks/kernels/bench.c: index arrays typed `misc_index_t` so the
   correctness check passes on both ISAs and the range-vs-indices rows
   read whole index words; run it on the Mac and record the corrected
   rows in the bench's own header or README where it keeps them.

## Verification

    cd tests/cpp && make && ./test_bartcore
    R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'
    R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-f0236082.rds
    (bcf and multinomial likewise; then record, then compare against the new files: 52/12/11 identical)
    R_LIBS=<lib> Rscript benchmarks/R/bench-sampler.R compare benchmarks/baselines/bench-sampler-127f04ee.csv
