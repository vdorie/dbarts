# constant-leaf-fits

agent: Opus
rng: neutral (structural: identical values consumed in identical order)
budget: ~700 changed lines (stop and report past ~1050)

## Goal

Replace the constant-leaf per-tree fits slab (numTrees x n doubles) with
per-tree node-indexed mu tables plus a per-tree uint32 leafOf map. The
draw-time scatter (misc_setIndexedVectorToConstant) disappears, the hot
sweep reads 4 bytes per obs-tree instead of 8, and slab-wide ops become
mu-table ops. Linear (vector-param) and GP (function-param) forests keep
the dense slab unchanged.

## Context

- The scalar draw already lands in forest.paramByNode (node-indexed,
  chain.hpp:2050-2057) and then scatters it into the slab
  (chain.hpp:2066-2072); the test path already reads paramByNode
  directly (chain.hpp:2075-2080). The compact form persists that table
  per tree and deletes the scatter.
- Hot consumers, all value-preserving under mu[leafOf[i]] gathers with
  loop order unchanged: the fused residual rolls (chain.hpp:618-634 and
  the no-recording twin near chain.hpp:861-898 - keep the two in
  lockstep), the totalFits rebuild (chain.hpp:669-677), prior refresh
  (chain.hpp:956-990), the mutation add/sub pair (chain.hpp:1068-1071),
  the restore rebuild (chain.hpp:1117-1153).
- Slab-wide consumers that become mu-table ops, bitwise-identically:
  BCF glue scale (combiner.hpp:488 - scaling mu then gathering equals
  scaling each entry), level-centering tree-0 shift (combiner.hpp:756),
  drawGlue leaf-value read (combiner.hpp:459-467 - the entry IS mu).
- Dense materialization stays available for exports: the memcpy
  accessor (chain.hpp:495) and any other raw-slab reader gather into
  the caller's buffer (identical bytes). Audit the standalone call at
  chain.hpp:1872 (local fits buffer) and storeSavedTreeRecord
  (chain.hpp:1268), which runs while params are live and can read them
  instead of fits.
- leafOf invariant: leafOf_t[i] holds the bottom-node index owning
  observation i, correct whenever a consumer reads it. Maintenance is
  incremental in the accepted-move paths (members of the nodes a
  birth/death/change repartitions - work proportional to move size, not
  n) and by rebuild wherever partitions change wholesale (data
  mutation/refreshSubtree, state restore, construction: all-root is all
  zeros). Seam choice (MoveContext member vs Tree method parameter) is
  the implementer's; rejected moves that restore partitions exactly
  must leave leafOf untouched.

## Constraints

- Constant leaf only: dispatch on the existing hasVectorParams /
  hasFunctionParams traits (if constexpr; an empty dense vector on the
  scalar path is acceptable - no Forest template specialization).
- leafOf is uint32 holding node indices. Do NOT narrow to uint16/uint8:
  a tree's node count is bounded only by n, and the width buys little
  once the scatter is gone; narrowing is a C6-profile-gated follow-up.
- Empty leaves keep their forced-zero params in the mu table; consumer
  guards (occupancy skips) stay as they are.
- No dbarts.h or inst/include change; no state-format change (fits and
  leafOf are both derived state, rebuilt on restore).
- No R-file change expected; if one becomes necessary, stop and report.
- bench-sampler is orchestrator-run; do not run it.

## Steps

1. Forest storage: persistent per-tree mu tables (the paramByNode flow,
   sized to each tree's node array) + the leafOf slab; scalar path only.
2. leafOf maintenance per the invariant above, plus a tests/cpp check:
   after a randomized move sequence, leafOf matches a map freshly
   derived from fillBottom + index segments.
3. Rewrite the sweep rolls, totalFits rebuild, prior refresh, mutation
   add/sub, and restore rebuild as order-preserving gathers.
4. Scalar sampleParametersAndSetFits / setTreeFitsFromParameters write
   the mu table and drop the scatter; audit the listed call sites.
5. Combiner: glue scale, tree-0 shift, drawGlue read via mu tables.
6. Dense-export accessors materialize by gather.

## Verification

R CMD INSTALL --preclean .; tests/cpp make clean && make &&
./test_bartcore; tinytest full suite (3050, zero failures);
equivalence.R compare vs equivalence-ac6ec2c.rds (22/22 identical
draws); bcf-equivalence.R vs bcf-equivalence-99205ee.rds (5x6 bitwise);
multinomial-equivalence.R vs multinomial-equivalence-2bd34db.rds (3x5
bitwise). Any divergence: stop, report scenario/channel, change nothing.

## Landings

C2a be5091d (2026-07-17). The compaction: muByTree + uint32 leafOf land
for the constant leaf, every consumer converted to order-preserving
gathers (sweep rolls, totalFits rebuilds, prior refresh, mutation
add/sub, restore, BCF glue scale, level shift, drawGlue read, dense
exports); linear/GP keep the dense slab via the leafIsConstant trait.
All gates bitwise (suite 3050/0, 22/22 + bcf 5x6 + multinomial 3x5);
reviewer re-ran component binary, equivalence, and bcf independently.
Diff 388 lines. DEVIATION HELD OVER: leafOf is rebuilt by an O(n)
scattered-write pass per tree per sweep inside
sampleParametersAndSetFits / setTreeFitsFromParameters - a correct
gate-green base, but scattered writes are cache-line-granular, so this
retains most of the scatter share the plan targets.

C2b 252c625 (2026-07-17). The seam move: metropolisJumpForTree grows a
default-null changedNode out-parameter written at the five accept
sites; on stepTaken the chain patches leafOf below that node only
(updateLeafOfBelow, one fillBottom walk - O(move size), not O(n));
rejections write nothing; grow-from-root, data mutation, and the
restore paths keep the full rebuild; the per-sweep scattered-write
rebuild is deleted from the draw and set paths. One subtlety beyond
the sketch: sampleTreesFromPrior replaces partitions wholesale but the
next sweep's residual roll still reads the old (mu, leafOf) pair as
the cached-fits evaluator, so it cannot rebuild in place - it marks
trees in a per-tree leafOfStale flag and the sweep rebuilds lazily
after the move settles, before the draw; every full-rebuild site
clears the flag. Test asserts leafOf against a freshly derived map
after each of 25 sweeps plus a prior-reset round. Diff 145 lines.
Gates all bitwise (suite 3050/0, 22/22 + bcf 5x6 + multinomial 3x5),
run twice: by the implementer and end-to-end by the reviewer (install
--preclean, component tests from clean, suite, all three compares).
bench-sampler compare PENDING the batched quiet window (with C1/C3).

Compare DISCHARGED 2026-08-04 on dbarts-bench (x86 Ryzen 3700X, full
mode, same-machine A/B vs 9047e05, the parent of the first commit;
attribution run against 252c625 isolates the refactor; arm64 remains
unquantified). Reproducible across all three runs: the intended win,
run-n10000-p10-t75 14-18% FASTER; and one real REGRESSION,
setPredictor-accept-n1000-t75 +22-28%, present at 252c625 - consistent
with the note above that data mutation keeps the full O(n*m) leafOf
rebuild. Follow-up item: docs/plans/archive/setpredictor-leafof-rebuild.md.
The other flagged metrics are noise: sub-millisecond metrics on that
box drift +-10% between identical invocations (embedded-offset spanned
0.99-1.10 across three runs of the same comparison), above the 5%
gate.
