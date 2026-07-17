# review-perf-followups

agent: Opus (every commit; engine/numerics)
rng: C1-C4 neutral (C1 empirically - see Constraints), C5 shifting on
  the multinomial channels only, C6 measurement-only
budget: C1 ~350 lines (widened, below); C3 ~150; C4 ~150; C2/C5 carry
  their own sub-plans

## Goal

Land the package review's Tier 5 performance findings and Tier 4 engine
notes (source record: package-review-remediation.md Tiers 4-5), neutral
commits before the draw-moving multinomial restructure, then close the
parked NT-store question with a real large-n profile.

## Context

- C1 (P2): the leaf marginals carry sum wz^2 (z'Wz), a model-free data
  constant that is additive over any partition of a fixed member set and
  so cancels in every MH comparison; scan.hpp:23-27,68-70 already omits
  it. Drop it from ALL THREE marginals - constant
  (model.hpp:109-126,146-154, concept at :42), linear (model.hpp:299,
  318, computeProjection output :371-412), and GP (:336-338) - because
  the GP leaf delegates oversized nodes to the constant marginal
  (model.hpp:693-698): a constant-only drop would mix z'Wz-free and
  z'Wz-carrying marginals across the maxLeafSize_ boundary and change
  the GP posterior. Also: tree.hpp:177,497-523,800-801; moves.hpp:
  215-216; the four kernel signatures moments.c:311-410 /
  include/misc/stats.h:42-45 (internal, scalar, undispatched). The
  cross-model marginal tests (test_model.cpp) must pass UNEDITED once
  all marginals agree; an edit there means the drop is inconsistent.
  No serialization field carries it.
- C2 (P1): per-leaf mu plus a per-observation leaf id replace the
  per-tree n-vector scatter (treeFits slab ~8x smaller; kills the
  recorded ~18.5% setIndexedVectorToConstant share). Sub-plan first.
- C3 (P4): birth-move child-by-subtraction and partition/suffstat gather
  fusion, tree.hpp:757-783. Order-careful; rng-neutral.
- C4 (Tier 4): E3 route combinedFits (combiner.hpp:668-682) through
  softmaxLocationMajor (combiner.hpp:561), bitwise-neutral only if the
  arithmetic is identical - verify; E5 guard the n==0 OOB read at
  model.hpp:2074-2079 (or prove the upstream reject); E1 confirm the
  exact gate exercises unequal per-forest shrinkage in the
  level-centering conditional (combiner.hpp:729-748) - investigation,
  change only on evidence.
- C5 (P3): multinomial margins O(nK^2) -> O(nK) via prefix/suffix
  log-sum-exp (NOT subtract-exp: cancellation when the excluded category
  is argmax), plus K-contiguous fit locality and a shared exp cache.
  Sub-plan first.
- C6: n=1e6 profile on the bench box after C2/C3, then the recorded
  x86-simd rule: restore NT-store variants behind a size threshold if
  the store-bound elementwise share is material, else delete the
  commented bodies, numbers in the landing note either way.

## Constraints

- Sequence C1 -> C2 -> C3 -> C4 -> C5 -> C6; one implementer at a time.
- C1's neutrality is a review claim, not structural: the sum wz^2
  rounding sits inside each node's marginal and cancels only
  algebraically in the branch-sum ratio (moves.hpp:47-63,507-517). If
  any anchor moves: STOP and report. Reclassification to shifting
  (posterior unchanged by the algebra) and any re-record are the
  orchestrator's call, never the implementer's.
- No dbarts.h diff anywhere in this arc.
- C5 re-records only multinomial-equivalence; the ac6ec2c and bcf
  anchors stay bitwise through the entire arc.
- Hot-path commits (C1-C3) need bench-sampler compare on a quiet
  machine before their records close; the orchestrator requests the
  window - implementers never run it.

## Steps

1. C1: drop sumWeightedResponseSq end to end - node field, the four
   kernel signatures, marginal signature, concept, move copy, birth
   collapse; update scan.hpp's omission comments to match.
2. C2: write and land the sub-plan, then implement.
3. C3: fuse birth partition and child stats; child-by-subtraction.
4. C4: E3 and E5 changes; E1 investigation recorded in the landing note.
5. C5: write and land the sub-plan, then implement with the re-record.
6. C6: profile, decide by the recorded rule, restore or delete.

## Verification

Per commit: R CMD INSTALL --preclean .; cd tests/cpp && make clean &&
make && ./test_bartcore; full tinytest::test_package("dbarts");
Rscript benchmarks/R/equivalence.R compare
benchmarks/baselines/equivalence-ac6ec2c.rds (22/22 identical draws);
bcf-equivalence.R compare vs bcf-equivalence-99205ee.rds (5x6 bitwise);
multinomial-equivalence.R compare vs multinomial-equivalence-2bd34db.rds
(3x5 bitwise; C5 instead re-records with a neutrality trail on the
untouched channels). air format --check . on any R touch. bench-sampler
compare vs bench-sampler-60a13b6.csv for C1-C3 (quiet window, no
concurrent load).

## Landings

C1 033c3cc (2026-07-17). Drop landed for ALL THREE marginals, not the
constant leaf alone: the first, constant-only attempt failed the
component tests because the GP leaf's maxLeafSize_ fallback delegates
oversized nodes to the constant marginal - mixing a z'Wz-free fallback
with z'Wz-carrying GP marginals would have changed the GP posterior.
The widened shape: constant and linear marginals drop the term
structurally (linear's marginal becomes + 0.5 b'M^-1 b / sigma^2); GP
computes z'Wz and adds it back explicitly because z'V^-1 z entangles it
(Woodbury), which keeps the fallback coherent. Node loses the field,
the four suffstat kernels lose their third output, moves.hpp loses the
restore copy. The cross-model marginal tests passed UNEDITED after the
widening (the built-in consistency check the amendment demanded);
absolute-value oracle tests updated by adding the analytically exact
dropped term in-test. Gates: install --preclean, tests/cpp from clean,
tinytest 3050/0, equivalence 22/22 identical draws (the empirical
neutrality question resolved: no knife-edge accept flipped), bcf 5x6
and multinomial 3x5 bitwise; reviewer re-ran the component binary and
the 22-scenario compare independently. Diff 316 lines. bench-sampler
compare PENDING the batched quiet window (with C2/C3).

C3 (2026-07-17): NO-GO, no code landed; the record is this note.
Child-by-subtraction diverged at the component gate (right-child
marginal off by the FP reorder; parent - left is not the direct
reduction) and was reverted at first divergence per protocol - it
also couples birth() to a live parent-stats invariant several tests
do not satisfy. Pass fusion audited as structurally non-neutral: the
suffstat kernels (moments.c:331-405) reduce remainder-first in 5-way
groups over FINAL ascending buffer positions, while every partition
variant (SIMD/scalar two-pointer, mask, wide-mask, sparse, MIA) is
swap-based and does not finalize elements in that order, so in-loop
accumulation cannot reproduce the kernel's bracketing bitwise; this
holds uniformly across weighted/categorical/missing paths.
refreshSubtree additionally has stale parent stats during a full
refresh. A shifting-class rescue (re-record everything for an
estimated 1-3%) judged not worth it with C2's scatter kill landed.
Both P4 fusion ideas are hereby closed; do not reopen without a
fixed-lane redesign of the suffstat kernels (simd-survey.md #3
territory).

C4 8a6a0c2 (2026-07-17). E3 ROUTED: combinedFits now gathers into
combined_ and applies softmaxLocationMajor in place, the documented
combinedTestFits pattern; arithmetic verified identical and the
multinomial train channel bitwise confirms. E5 GUARDED: rescale()
early-returns the identity scale (min 0, max 0, range 1) at n == 0 -
the OOB seed read was reachable from construction and
setResponse/setData/setOffset, none of which reject n == 0. E1
CLOSED, NO ACTION, by evidence: multinomial-exact.R arms 1-5 all use
symmetric priors (equal numTrees/nodeScale/k per category), so v_k is
equal in every arm and invV cancels from the level-centering
conditional; categorical-exact.R is single-forest categorical-split,
never reaching afterCombine; and the gate is blind in principle -
every arm compares identified softmax probabilities, invariant to the
grand shift by construction, while the test surface cannot even
build unequal per-category configs (one n.trees for all K). The
conditional therefore affects only the unidentified redundancy
direction (mixing), the design comment already documents the
observation-count approximation as unbiased for identified
quantities, and a leaf-count-exact refinement would be a
shifting-class change with no measured mixing deficit - not taken.
Gates: suite 3050/0, 22/22 + bcf 5x6 + multinomial 3x5 bitwise
(implementer full battery; reviewer re-ran multinomial + component).
Diff 25 lines.

C6 (2026-07-17, delete landed with this note's commit). The bench-box
empirical profile was BLOCKED by machine permissions
(perf_event_paranoid=4, ptrace_scope=1, no passwordless sudo, no
valgrind; Rprof cannot see below the single .Call frame) - the box
also lacks the libtirpc dev symlink tests/cpp needs, though the
22-scenario equivalence compare passed statistically there (max |z| =
0.00, the cross-host contract). The decision fell to enumeration,
which is the stronger argument: after the C2 fits compaction, NO
store-bound elementwise misc_ kernel runs per-sweep on the
constant-leaf path (rolls and totals are fused scalar gathers); the
surviving large-n callers are RMW (test-fit accumulation, BCF glue
scale - NT stores lose on RMW by construction) or cold/once-per-cycle
(restore fills, rescale passes). The consumers NT stores would
accelerate no longer exist, so per the recorded rule the commented
intrinsics bodies are DELETED (~820 lines across
linearAlgebra_{sse2,avx}.c; the live dispatch arms and transpose
kernels stay; x86 cross-target syntax check clean; these TUs do not
compile on arm64, CI x86 runners validate the real build).

Speed record for the arc (arm64 host, quiet window 2026-07-17).
Compare vs 60a13b6 pre-fix flagged setPredictor-accept-n1000-t75
0.252 -> 0.278 ms/update (+10%): the mutation transaction's
SIMD-contiguous passes became scalar gathers at the compaction. The
fused-walk commit fb4d76a (map rebuild + total add in one node walk)
RECOVERED it: the cold re-record reads 0.256 (+1.6%, in band). The
intermediate warm compare that also read 0.278 post-fix was a
measurement artifact - it drifted every arm including controls upward,
and the warm drift masked the fix; the cold record is authoritative
(the control-difference rule cut both ways here). All other arms
within +/- 2 percent cold. A/B at n=1e5 (60a13b6 vs tip, median of
5x40): ntree=75 19.375 -> 18.400 (0.950), ntree=200 48.325 -> 47.225
(0.977) - the compaction wins where VD's workloads live. Baseline
re-recorded cold at the tip as bench-sampler-b9d53c7.csv. ARC CLOSED
2026-07-17: P1 = C2ab, P2 = C1, P3 = C5, P4 = C3 no-go, Tier 4 = C4,
NT-store = C6 delete; every item dispositioned with its evidence.
