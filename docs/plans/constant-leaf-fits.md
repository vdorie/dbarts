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

C2b (pending): move the seam - metropolisJumpForTree exposes the
changed node index; on stepTaken the chain rewrites leafOf for the
changed node's member segment only (one fillBottom walk from that
node covers birth/death/change); rejected proposals write nothing
(accept-time only, after the move settles); grow-from-root, data
mutation, and restore keep the full rebuild; the per-sweep rebuild is
deleted; the leafOf consistency test asserts after every sweep of a
multi-sweep run. Budget ~250 lines; full gate battery, all bitwise.
If a session break leaves an uncommitted diff in the worktree, it is
unverified C2b work: re-run the battery before trusting it.
