# block-fusion-stage-b (b>1 fusion, the win, the one re-record)

DRAFT. Stage B turns the landed b=1 atom machinery into the block-fused
sub-sweep at b>1: build the block-static field once per block, carry S per atom,
and regroup the constant-leaf suffstat/draw/roll to O(atoms). It is the WIN
(design section 7, Stage B) and carries the single VD-approved re-record.

agent: opus
rng: shifting for the shipped default. Commits (i)-(iii) are NEUTRAL at the
  shipped default (b=1 stays exactly Stage A, bitwise) and add only EXACT or
  dormant b>1 machinery; the one reduction-order change (the affine O(atoms)
  regrouping) and the default flip to b>1 land together in commit (iv), the SOLE
  re-record carrier. Same shifting-class discipline as the suffstat/threading
  re-records: draws move, the posterior does not.
window: Stage B of docs/design/block-fusion.md (frontier 3.4). Branch from the
  current bartcore tip (Stage A landed: d9a5d39 landing note; b=1 atom path
  DEFAULT ON, bench-neutral). Stage A is the de-risk anchor this builds on: the
  flag becomes the b-dispatch, atomOf/A-cache/S-carry come back on, `members`
  becomes owned, and the suffstat/roll/draw regroup to O(atoms). Constant-leaf
  Gaussian ONLY (design 5.3); linear/GP stay on the legacy writer. No dbarts.h
  change: b is an internal runtime knob, not new public surface.
budget: ~1200-1800 lines across ~5 files. atoms.hpp grows the joint b>1 map +
  the affine suffstat + the owned buffer (~500-700); chain.hpp the block driver +
  the block-boundary g/O bookkeeping (~250-400); moves.hpp the b>1 move
  maintenance (~150-250); test_atoms.cpp the b>1 fuzz + cross-ISA suites
  (~400-600); bench-sampler.R the large-n/b grid (~80-150). Clean abort at every
  commit: default blockSize = 1 reverts the shipped engine to Stage A.
bench: VD grants the quiet machine (MEMORY: bench-sampler on dbarts-bench, x86
  AVX2). The ~6x DRAM drop is a memory-bandwidth effect that shows on the x86
  box, not the arm64 laptop. Headline gate at commit (iv) over the extended grid;
  see section 9 for the baseline decision.

---

## 0. What Stage B is, in one screen

Stage A funnels every constant-leaf per-leaf suffstat through the b=1 atom path
(AtomMap::buildAggregateWrite, chain.hpp:798-806): an atom is one leaf, its
members alias tree.indices, and its (A, G, Q) are the node cache computed
BITWISE by the same misc kernel. Stage B lifts this to a BLOCK of b consecutive
trees. An atom becomes a cell of the joint partition of the b trees (design 2.1);
its state is (A, G, S) plus the b-tuple of leaves, and the constant-leaf math
collapses per-observation field maintenance to per-atom algebra:

- ONCE per block (O(n)): build g_i = w_i(y_i - O_i) from the frozen outside-block
  fit O, scatter into per-atom G(c), A(c); seed S(c) from the block trees' fits.
- INTERIOR (b trees x O(atoms), no O(n) pass): every leaf suffstat, leaf-mean
  draw, residual roll, and rejected-move restore is atom algebra over (A, G, S).
- ONCE per block (O(n)): scatter the b drawn leaf means into treeFits and update
  the running total so the sigma/latent stages and the next block see a correct
  residual.

The suffstat regroup (design 3.1) is the ONLY floating-point change: a leaf's
sumWeightedResponse becomes D(L) + mu*W(L) with D(L) = sum_{c in atoms(L)} r(c),
r(c) = G(c) - A(c)S(c). This equals today's index-order gather in EXACT
arithmetic (G, A, S are exact partial sums of the same terms) but groups the
floating-point sum differently, so draws SHIFT (same posterior, different bytes).
Everything else Stage B adds is EXACT (g, S, scatter) or dormant at the shipped
default until commit (iv) flips it on.

## 1. The b-dispatch and the block driver

b is a runtime knob. It replaces the role the compile-time `BARTCORE_BLOCK_FUSION`
flag played in Stage A (which only chose atom-path-vs-legacy). Keep the compile
flag as the master enable of `useAtomSuffstatSource` (chain.hpp:379-382); add a
runtime `std::size_t blockSize` per Forest (next to `numTrees`, chain.hpp:261),
defaulted at Chain construction by the n-adaptive rule (section 8, decision 1).
`blockSize == 1` is EXACTLY Stage A's path (bitwise); `blockSize > 1` engages
fusion. The atom path only exists when `useAtomSuffstatSource` is true
(constant-leaf + compile flag), so a linear/GP forest ignores blockSize entirely.

run()'s tree loop `for (size_t t = 0; t < forest.numTrees; ++t)`
(chain.hpp:771) becomes a block loop over `ceil(numTrees / blockSize)` blocks,
block `[t0, min(t0+blockSize, numTrees))`. At blockSize == 1 each block is one
tree and the body is byte-for-byte today's per-tree sequence (roll ->
buildAggregateWrite -> metropolisJumpForTree -> sampleParametersAndSetFits ->
per-tree roll bookkeeping). At blockSize > 1 the body is: block-entry g/O/S build
(section 3) -> for each tree in the block (interior suffstat -> move ->
draw -> S-carry roll, section 4) -> block-exit fit scatter + total update
(section 3). The last block may be short (numTrees not divisible by b); the atom
math is agnostic to b, so a short final block just fuses fewer trees.

The Stage-A run() wiring that parks the Stage-B bookkeeping at b=1
(`aCacheBypass = true`, `trackAtomOf = false`, chain.hpp:754-755, and the inert
S carry note chain.hpp:821-825) becomes conditional on blockSize: OFF at
blockSize == 1 (Stage A), ON at blockSize > 1.

## 2. The owned per-block atom map (design 2.3, 2.4)

Stage A resolved DESIGN A: at b=1 `members` aliases tree.indices and no owned
buffer exists (atoms.hpp:53-56, 74-75). b>1 needs the owned atom-ordered buffer:
one `atomMembers` per block (a permutation of 0..n-1 grouped by atom, design
2.3), NOT b per-tree buffers. Add it to AtomMap as an owned `std::vector<size_t>`
alongside the aliased `members` pointer; at blockSize == 1 keep aliasing
tree.indices (the shipped path is unchanged), at blockSize > 1 point `members`
at the owned buffer.

Re-enable the three Stage-A-parked mechanisms for b>1 (they are Stage B's, kept
default-on for the component tests, disabled only on the shipped b=1 path):

- `atomOf` (obs -> atom, atoms.hpp:61, `trackAtomOf`): the move bookkeeping that
  lets a split/merge find an observation's atom in O(1). Required at b>1.
- the cross-sweep A cache (atoms.hpp:98-142, `aCacheBypass`): at b>1 the joint
  map is PATCHED not rebuilt each sweep (design 4.5, persistence), so A(c) served
  from the memcmp-validated cache saves the O(n) re-scan the b=1 monolithic
  kernel could not. Re-enable it (aCacheBypass = false at b>1).
- the S carry (atoms.hpp:72, `setInBlockFits`): S(c) = sum over block trees of
  their fit on c, the Gauss-Seidel coupling carried per atom.

Joint-map build (buildForBlock, new): route each observation through the b block
trees, form its leaf b-tuple, and bucket it into an atom. This is the O(bn)
block-init / block-composition-change / periodic-rebuild path (design 4.5); the
steady-state map PERSISTS and is patched by moves. The b>1 move maintenance
generalizes the landed b=1 kernels: splitAtom (atoms.hpp:479-565), mergeAtoms
(631-657), refreshSubtree (766-792), undoSplit (573-588), snapshotSubtree/
restoreSubtree (680-725) each currently touch a single leaf slot; at b>1 a move
in block tree t_j slices only the t_j coordinate of the leaf-tuple, so an atom
splits/merges against the OTHER b-1 coordinates held fixed (design 4.1-4.3). The
leafTuple SoA (atoms.hpp:66) widens from one int per atom to b ints per atom.
The atom member partition still reuses the tree's OWN partitionChildren over the
owned buffer (design 5.2, atoms.hpp:497, 772); no second partitioner is forked.

## 3. Block boundaries: the two O(n) passes (design 2.1, 3.5) -- EXACT

Block ENTRY builds the block-static field, all EXACT (no reduction-order change):

- O_i = F_i - sum_{t in block} treeFits[t*n + i], where F is the running
  full-forest fit (maintain F across blocks; it is today's `totalFits` kept
  current per block instead of rebuilt once per sweep at chain.hpp:843-851).
- g_i = w_i(y_i - O_i) (unweighted: g_i = y_i - O_i). One O(n) scatter into the
  per-atom G(c) = sum_{i in c} g_i and A(c) = sum_{i in c} w_i.
- S(c) = sum over block trees of the tree's OLD fit on c (from treeFits), so the
  interior starts from the correct in-block fit.

Block EXIT is the one O(n) scatter per block (design 3.5): for each block tree
and leaf, write the drawn mu into treeFits[t*n + i] over the atom's members (one
walk of atomMembers writes all b trees' fits, cache-friendly; or call the
existing setTreeFitsFromParameters chain.hpp:2105 b times). Maintain a running
full fit F incrementally (F += new block fits - old block fits) so the NEXT
block's O is available without an O(bn) rescan. Do NOT feed incremental F to the
sweep-end stages: keep today's once-per-sweep fresh totalFits rebuild
(chain.hpp:843-851) so the sigma/SSR draw (chain.hpp:865-866) and refreshLatents
(chain.hpp:856) read the residual materialized EXACTLY from treeFits (design 3.7),
never an incrementally-accumulated sum that can drift over sweeps. Incremental F
is a block-local O accelerator only; if its drift ever matters, recompute a
block's O exactly as (fresh totalFits - the block trees' fits).

At blockSize == 1 this section is BYPASSED: keep the existing per-tree fused roll
(chain.hpp:774-787) and the post-loop total rebuild (chain.hpp:843-851). The b=1
residual path is untouched, so Stage A stays bitwise. The block-boundary g/O/S
build engages only at blockSize > 1 and is component-tested for exact equality
against a reference per-tree computation.

## 4. The interior affine O(atoms) reformulation (design 3) -- the re-record

This is the sole floating-point change. For block tree t_j and its leaf L:

- SUFFSTAT (3.1, replaces the per-tree buildAggregateWrite gather): iterate the
  atoms of L in the FIXED order (section 5) and accumulate
  D(L) = sum_{c in atoms(L)} r(c), r(c) = G(c) - A(c)*S(c); then W(L) = sum A(c)
  and sumWeightedResponse(L) = D(L) + mu_{t_j}(L)*W(L). Write (W(L),
  sumWeightedResponse(L)) into the node suffstat cache -- the SAME seam every
  constant-leaf consumer reads (atoms.hpp:365-382, writeNodeCaches). Feed the
  UNCHANGED model.hpp math: logIntegratedLikelihood (model.hpp:109-122) and
  drawFromPosterior (model.hpp:128-143) are byte-for-byte the same functions of
  the suffstat.
- DROP Q (fact 1.2): sumWeightedResponseSq cancels in every move ratio and is
  unused by drawFromPosterior, so the b>1 atom carries only (A, G, S). The node
  cache's sumWeightedResponseSq is fed 0; logIntegratedLikelihoodForNode still
  reads it (model.hpp:146-153) but the centeredSumOfSquares term is additive over
  any partition, so it cancels identically in logLikelihoodForBranch's
  new-minus-old difference whether Q is the true raw sum (Stage A) or 0 (b>1).
  Dropping Q is itself draw-NEUTRAL; it rides in this commit only because it is a
  b>1-only change (Stage A keeps Q for the b=1 bitwise anchor).
- DRAW (3.2): drawFromPosteriorForNode reads the node cache, so
  sampleParametersAndSetFits (chain.hpp:2146-2206) is UNCHANGED -- one standard
  normal per non-empty leaf, in fillBottom leaf order, same RNG consumption.
- ROLL (3.3, replaces the O(n) treeY roll): after tree t_j's draw, for each leaf
  L and each atom c in atoms(L), S(c) += (mu_new(L) - mu_old(L)). O(atoms) total;
  no per-observation field touched. Tree t_{j+1}'s suffstat then reads S(c) and
  subtracts mu_{t_{j+1}}(c), automatically seeing t_j's new fit (Gauss-Seidel).
- RESTORE (3.4): a rejected structure move only dirties the atom partition and
  the derived leaf suffstats, never S or G of an unsplit atom (params are drawn
  after the move resolves). Restore = undo the atom split/merge (undoSplit /
  restoreSubtree, already landed for b=1; generalized to b>1 in section 2). No
  O(n) pass.

Activation: this commit sets the default blockSize to the n-adaptive value
(decision 1), which is the FIRST time the shipped default reaches the affine
path. That default flip is the reduction-order change for the shipped engine and
the sole re-record trigger. Clean abort: set the default back to 1.

## 5. Machine independence at b>1 (hard gate)

The invariant (MEMORY: dbarts simd reproducibility): same seed => same draws on
ANY CPU. It holds because the draw-path reductions are SCALAR fixed-order:
misc_computeIndexedSufficientStatisticsFast (moments.c:334) is a plain mod-5 +
stride-5 loop, NOT SIMD-dispatched (only the mean/variance helpers dispatch,
moments.c:1330+). The per-atom (A, G) already come from that scalar kernel, so
they are byte-identical across ISAs.

Stage B adds ONE new reduction that feeds the draw: the tree-level sum
D(L) = sum_{c in atoms(L)} r(c). It MUST be a scalar, fixed-order accumulation:

- FIX the atom-iteration order for D(L). Enumerate atoms(L) in a deterministic
  order that is identical on every ISA. Atom ids are assigned by the same
  RNG-driven move sequence on every machine, so any id-derived order is
  cross-ISA stable; pin it explicitly (e.g. ascending atom id, or an atom-ordered
  per-leaf list built in fillBottom order) and never off a hash/pointer order.
- Enumerate atoms(L) efficiently: maintain a leaf -> atoms grouping (per-leaf
  atom list, or atoms sorted so a leaf's atoms are contiguous) so the sum is
  O(atoms in L) not O(all atoms) per leaf. The grouping order IS the fixed
  reduction order.
- Never SIMD-dispatch or thread the D(L) sum. It is O(K) and cheap (K ~ 10^3);
  keep it scalar so FMA/vector-tree reassociation cannot reorder it.

Gate: a tests/cpp component test forces each instruction set via
misc_stat_setSIMDInstructionSet (moments.c:1330) -- scalar / SSE2 / AVX2 / NEON
as available -- runs a b>1 sweep from a fixed seed, and asserts the drawn
parameters are byte-identical across all forced ISAs. This is the invariant the
whole no-backwards-compat reproducibility contract protects.

LANDING NOTE: the equivalence anchor (section 7) is a within-host bitwise
check. Cross-host comparison is only statistical, because the scenario data
itself (friedman's sin(), platform libm) differs in the last ULP across OSes --
true already at b=1, not a fusion regression.

## 6. Commit-by-commit

Each compiles, gates, and aborts cleanly (default blockSize = 1 => Stage A).

### (i) Block driver + runtime blockSize knob. NEUTRAL. ~200 lines.
- Add `std::size_t blockSize = 1` to Forest; the n-adaptive default rule (returns
  1 for now, so nothing changes yet). Rewrite run()'s tree loop (chain.hpp:771)
  as a block loop; at blockSize == 1 the body is byte-for-byte today's per-tree
  sequence. Make the Stage-A b=1 parking (aCacheBypass, trackAtomOf,
  chain.hpp:754-755) conditional on blockSize.
- GATE: NEUTRAL -- equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds,
  tinytest full pass NO snapshot regen, tests/cpp clean. Confirmatory bench
  (harness/machine pin). The block loop is a pure refactor of the tree loop.
- LANDED: 9d2d3dd. GATE held as planned -- equivalence 22/22 IDENTICAL,
  tinytest 2728 pass with NO snapshot regen.

### (ii) Joint b>1 atom map + owned member buffer + re-enabled bookkeeping. NEUTRAL at default. ~450 lines.
- Owned atomMembers buffer; widen leafTuple to b per atom; buildForBlock (route
  n obs through b trees into atoms); generalize splitAtom / mergeAtoms /
  refreshSubtree / undoSplit / snapshot/restore to slice one leaf-tuple
  coordinate (design 4.1-4.3). Re-enable atomOf/A-cache/S-carry at b>1.
- GATE: NEUTRAL at default (b=1 path untouched). tests/cpp: a b>1 mutation fuzzer
  (extend testAtomFullMoveFuzz) asserts the patched map matches a from-scratch
  buildForBlock rebuild after EVERY move (bitwise A/G, matching leaf-tuples,
  atomOf round-trip), and a rejected move restores the map bitwise. No draw path
  yet, so equivalence/tinytest stay identical.
- LANDED: f494c72. Fuzzer gate held -- the patched map matches a from-scratch
  buildForBlock rebuild after every move and on rejection, bitwise. DEVIATION:
  death/change/swap maintain the atom map by an O(n) regroup rather than an
  in-place splice; flagged as a perf suspect for the large-n bench (section 9).

### (iii) Block-boundary g/O/S build + fit scatter-back + running total. EXACT; NEUTRAL at default. ~300 lines.
- Block-entry O/g/G/A/S seed and block-exit scatter + F update (section 3),
  engaged only at blockSize > 1; blockSize == 1 keeps the per-tree roll
  (chain.hpp:774-787) and post-loop rebuild (chain.hpp:843-851).
- GATE: NEUTRAL at default. tests/cpp: at b>1, assert g/G/A/S seeds and the
  block-exit treeFits + running total equal a reference per-tree computation in
  EXACT arithmetic (these steps are exact, not a reduction-order change).
- LANDED: 9d16711. Component gate held -- entry seeds and exit scatter matched
  the per-tree reference exactly. DEVIATION found later: blockStaticField
  computed g = w*(y - O), double-counting weights against the weighted leaf
  kernel (which applies w once); correct is g = y - O. Never shipped -- the
  method was dormant at blockSize 1 -- fixed in (iv-a) below.

### (iv) Interior affine regroup + Q-drop + order-fix + default flip. THE RE-RECORD. ~350 lines.
- Interior suffstat (3.1), S-carry roll (3.3), rejected-move restore (3.4) in
  atom terms; drop Q (1.2); fix the atom-iteration + reduction order for cross-ISA
  bitwise (section 5); set the default blockSize to 4 / n-adaptive 8 (decision 1).
- This is the sole draw-shifting commit; gates in section 7. Clean abort: revert
  the default to 1 (leaves all machinery in place, ships Stage A).
- SPLIT: landed as (iv-a), the neutral machinery, and (iv-b), the default flip;
  (iv-b) has NOT landed (see section 7 BENCH status).
- LANDED (iv-a): c4fc3f3. Neutral only -- affine O(atoms) interior, run()
  wiring, drop Q at b>1, cross-ISA fixed-order reduction -- default blockSize
  stays 1. GATE: affine-identity max |diff| 3.55e-13; cross-ISA bitwise
  (scalar/NEON/SSE2/AVX2); equivalence 22/22; tinytest green. Also fixed the
  (iii) g-weighting deviation. Confirmed section 3's plan: sweep-end totalFits
  is rebuilt fresh from treeFits, never the incremental running fit.
- FOLLOW-ON bee4354: DBARTS_BLOCKSIZE env override -- the scenario knob that
  section 7's FORCE-b>1-ON-THE-ANCHORS precondition needs. Neutral when unset.
- FOLLOW-ON a694ec8: forcing b>1 at small n (n.trees=1) via bee4354 surfaced a
  real bug -- the six atom-map move hooks dispatched on blockWidth > 1, so a
  fused block of WIDTH 1 (tree count < blockSize, or a short final block) fell
  into the b=1 kernels, which assume atom members alias tree.indices and that
  the node cache carries raw Q; buildForBlock owns the members buffer and
  writeAffineNodeCaches seeds Q = 0. Wrong leaf offsets plus an uncancelled Q
  penalty on births under-fit the forest (categorical-exact gap 0.093 vs 0.004
  at b=1). Fix: dispatch the hooks on an explicit blockMode flag set by
  buildForBlock, cleared by the per-tree/aggregate builders. New gate
  testBlockAffineIdentityUnderMoves (widths 1, 2, 3; max |diff| 4.4e-13, fails
  within one birth if the old dispatch is restored); categorical-exact gap at
  b=4 recovers to 0.0042. Width >= 2 was always correct -- the bug was the
  width-1-fused edge the eventual (iv-b) flip would hit. b=1 path unaffected
  (equivalence 22/22).

### (v) Harden for the Stage-B gates. ~250 lines.
- Small-n fallback to blockSize = 1 below a cutoff (like the SIMD toggles);
  finalize the n-adaptive rule; cross-ISA CI test (section 5); extend
  bench-sampler.R with the large-n/b grid and record the fresh baseline
  (section 9). Confirm b-invariance is only distributional (equivalence anchor
  fixes one b; decision 4).
- GATE: bench across the grid on the quiet machine; equivalence/tinytest stable
  (b is a perf knob once the fixed order is set). Deeper b-sweep, rebuild cadence,
  and fragmentation tuning are Stage C (design section 7), not this plan.
- LANDED (partial): b74fec9 adds the bench-sampler.R biggrid mode (n in
  {1e4,1e5,1e6} x numTrees in {75,200} x blockSize in {1,4,8} via
  DBARTS_BLOCKSIZE); default invocation unchanged. Still pending: the small-n
  fallback, the finalized n-adaptive rule, cross-ISA CI wiring, and the fresh
  baseline recording -- all wait on the headline bench verdict (section 7).

## 7. Gates for the re-record commit (iv) and beyond

Shifting-class gates (docs/plans/README.md; design 7 Stage B):

- EQUIVALENCE, statistical: re-record the anchor at a FIXED b (decision 4) and
  gate matched-RNG DISTRIBUTION equality (z-mode, |z| < 4) of the new build vs
  equivalence-ac6ec2c.rds. NOT bit-identical (expected). Name the new baseline
  equivalence-<hash>.rds and update benchmarks/baselines/MANIFEST (mark ac6ec2c
  historical, the new file current).
- EXACT-POSTERIOR: the known-posterior / SBC scenarios already in the suite must
  hold -- benchmarks/R/sbc.R (gaussian, probit), categorical-exact.R,
  logistic-reference.R, and bcf-exact.R (constant-leaf BCF). These prove the
  regroup did not move the stationary distribution.
- FORCE b>1 ON THE ANCHORS (gate-validity precondition): the equivalence and
  exact-posterior scenarios run at small n, where the n-adaptive default
  (decision 1) and the small-n fallback (commit v) would silently revert them to
  blockSize = 1 -- making the re-record a no-op that hides the shift and proves
  nothing about the regroup. Both gates MUST pin blockSize to the fixed anchor b
  explicitly (a scenario/control knob that overrides the size heuristic), so they
  exercise the b>1 affine path at small n. Without this the entire re-record is
  untested.
  LANDED: bee4354 provides exactly this override (DBARTS_BLOCKSIZE).
- MACHINE INDEPENDENCE: the cross-ISA bitwise component test at b>1 (section 5)
  -- scalar / SSE2 / AVX2 / NEON byte-identical. Hard gate.
- TINYTEST: regenerate the RNG-locked snapshots ONCE by replaying whole test
  files (values depend on each file's full execution history, not just the
  preceding seed). Then full pass from a --preclean install.
- tests/cpp: the b>1 fuzzers ((ii)) + the exact block-boundary checks ((iii)) +
  the cross-ISA test, clean; delete stale binaries first (header edits, no dep
  tracking).
- BENCH (headline): bench-sampler on dbarts-bench (x86, quiet) over
  n in {1e4, 1e5, 1e6} x m in {75, 200} x b in {4, 8}, targeting the E1/E2
  ~6x DRAM drop -> multi-x wall-clock at n >= 1e5. Kill/scope-back if the
  realized speedup is far below the model (move-scan re-partition cost dominating
  at shallow trees, design 8.1). See section 9 for the baseline.
  STATUS: not yet run to a verdict. Early biggrid numbers show b=4 ~25x SLOWER
  than b=1 at n=1e4 -- expected small-n overhead (per-block fixed cost not yet
  amortized); the win must appear at n >= 1e5 or this is a kill signal.

## 8. Resolved decisions (VD-approved; recorded, not re-opened)

1. DEFAULT b = 4, or n-adaptive 8 at n >= 1e5. blockSize is defaulted at Chain
   construction from numObservations; a small-n cutoff falls back to 1.
2. v1 SCOPE = constant-leaf Gaussian only (covers stock BART, binary via latents,
   weights, BCF/grouped-RE forests). Linear/GP leaves stay on the legacy writer
   (they read the same node cache; useAtomSuffstatSource is false for them).
3. ONE one-time re-record, carried by commit (iv). Bundle it with
   within-chain-threading's re-record IF both land in the same window (one
   snapshot regeneration; both are shifting-class), else keep them independent.
4. The EQUIVALENCE ANCHOR fixes ONE b (record the new baseline at that b) and
   gates b-invariance (b=1 vs 4 vs 8) STATISTICALLY (z-mode), not bitwise.
5. GO STRAIGHT TO ATOMS: no standalone contiguous-layout item; the owned
   atom-ordered buffer + block-static g/w fields ARE the contiguous-per-block
   layout (design 5.1). Stage A already delivered the plumbing bitwise.

## 9. Decision for VD: bench harness + headline baseline

NEW (not covered by the five above). The headline gate needs a grid the current
harness does not have: bench-sampler.R caps at n = 10000, m in {75, 200}, no b
axis (benchmarks/R/bench-sampler.R:67-69), and the current speed baseline
bench-sampler-32fc7c8.csv has no large-n or b>1 rows to compare against.

Question: how to gate the ~6x headline claim.
Recommendation: (a) extend bench-sampler.R with an OPT-IN large-n/b grid
(n in {1e4, 1e5, 1e6} x m in {75, 200} x b in {1, 4, 8}) behind a CLI/env flag so
the default CI grid stays fast; (b) record a FRESH baseline on dbarts-bench (x86
AVX2 -- the DRAM-bandwidth-bound regime where the drop manifests; the arm64
laptop will not show it) and gate Stage B against THAT, reporting b>1 speedups vs
the b=1 column of the same recording (an in-recording A/B, immune to
cross-machine noise); (c) mark the new CSV current in MANIFEST.
Evidence that would change it: if the b=1 rows in the fresh recording differ
materially from 32fc7c8 on the shared small-n points, investigate before trusting
the large-n rows; if the x86 box is unavailable, the gate degrades to relative
b>1-vs-b=1 on whatever quiet machine is available (the DRAM claim is then
under-measured and must be flagged, not asserted).

STATUS: (a) has landed (b74fec9, the opt-in biggrid mode). (b) the fresh
dbarts-bench recording and (c) the MANIFEST update are still pending on the
headline verdict above.

## 10. Risks and out-of-scope

- MOVE-SCAN re-partition cost (design 8.1, the main way the win underperforms):
  shallow-tree change/swap re-partition near the root, ~2x the modelled 7% scan.
  The bench gate at (iv)/(v) measures it directly; the rebuild cadence is a
  Stage C lever.
- COMPLEXITY: the joint map is a second partition structure maintained in lockstep
  with b trees. Mitigation: Stage A de-risked the b=1 plumbing bitwise; the b>1
  generalization reuses the landed split/merge/snapshot/restore and the
  memcmp-validated A cache; the fuzzer ((ii)) fingerprints the map against a
  from-scratch rebuild after every move.
- LinkingTo CONSUMERS: stan4bart's mu forest is constant-leaf Gaussian, so its
  draws shift under the default flip (expected, shifting-class, within the
  re-record). Lockstep release (MEMORY: stan4bart walnuts plan); note, not a gate.
- OUT OF SCOPE: linear/GP atoms (design 5.3), within-chain-threading composition
  (design 6, but see decision 3 on shared re-record), the GPU seam, and the deep
  b-sweep / rebuild-cadence tuning (Stage C). No dbarts.h / ABI change in v1.
