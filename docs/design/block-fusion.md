# Block-fused sub-sweeps over persistent atom maps (frontier 3.4): engineering design

Status: CLOSED, WONT-DO. Stage A (b=1 exact refactor) landed then the whole
project was closed: Stage C measured the fused path 4-6x slower than b=1 (no
win, section 10, tip 10c59bc), and the machinery was excised at 5caf990 (kept
on archive/block-fusion). Original proposal (2026-07-11) follows. Turned the
falsifier-cleared CONCEPT in docs/design/parallel-bart-frontier.md section 3.4
into an implementable engine design. VD approved pursuing the flagship engine
perf project and accepting its one-time re-record (2026-07-11).

This note is grounded in the current engine (bartcore HEAD 2e2b1c9) and stays
consistent with the MEASURED falsifier numbers (docs/plans/parallel-falsifiers.md,
DONE 2026-07-08): E1 real-forest atom census and E2 field-fraction profile. It is
the companion build-out of the substrate note docs/design/data-layout.md, which
found the single-tree contiguous reorder is only a ~10% standalone lever and
explicitly named THIS construction ("the LARGE win on this primitive") as where
the real perf lives.

## 0. TL;DR of the mechanism (one screen)

Sweep the m trees in consecutive BLOCKS of b (default target b in {4,8}). For a
block B of b trees, an ATOM is a maximal set of observations that all b trees
route to the same b-tuple of leaves (the joint partition of the b trees). E1
measured atoms far below n at b in {4,8}: ~74-113 atoms at b=4 and ~1120-1596 at
b=8 for n=1e5, and atoms/n FALLS as n grows.

The engine's constant-leaf math (fact 1.2) is affine in the residual: a leaf's
only residual-carrying suffstat sumWResid_t(L) = D(L) + mu_t(L) sumW(L), where
D(L) is the leaf-sum of the full-forest residual field. Because every observation
in an atom shares the same in-block leaves, that whole per-observation field
collapses to a per-ATOM residual mass r(c). Then, inside a block:

- entry (once per block, O(n)): scatter the block-static field g_i = w_i(y_i -
  O_i) (y minus the frozen OUTSIDE-block fit) into per-atom masses G(c), A(c);
  seed per-atom in-block fit S(c);
- interior (b trees x O(atoms), NO O(n) pass): every leaf suffstat, leaf-mean
  draw, residual roll, and rejected-move restore is atom algebra;
- moves (~7% traffic, O(n_leaf)): a birth/change/swap slices atom interiors, so
  it re-partitions the affected members' raw g/w and splits their atoms;
- exit (once per block, O(n)): scatter the b trees' drawn leaf means back into
  the per-observation treeFits so the sigma/latent stages and the NEXT block see
  a correct residual.

So the ~85% "field maintenance" that E2 measured (residual roll + whole-tree
suffstat recompute + fit scatter) is paid ONCE PER BLOCK instead of once per
tree: an ~m/b reduction in O(n) DRAM passes, landing E1/E2's ~6x (7-10x low end)
modelled drop at b in {4,8}, n >= 1e5. The b=1 special case is atoms == leaves
and reproduces today's draws BITWISE (the de-risking anchor); b>1 regroups the
summation and re-records once (the "shifting" RNG class), machine independence
preserved by a fixed atom-iteration order.

## 1. What the engine does today (verified in code, the baseline to preserve)

The per-forest sweep is Chain::run (chain.hpp:728-787). Per tree t, in order:

1. Residual roll (chain.hpp:731-744): treeY carries a running residual across the
   sweep. Entering tree t, resid_i = y_i - sum_{k!=t} f_k(i) (new fits for trees
   already drawn this sweep, old for the rest). One fused O(n) pass retires the
   previous tree's new fits and admits this tree's old ones. treeFits[t*n..] is
   tree t's per-observation contribution (obs order, persistent); totalFits is
   the running sum, rebuilt after the loop (chain.hpp:779-787).
2. setNodeAverages -> computeLeafStats (tree.hpp:492-521): for each leaf, an
   O(n_leaf) GATHER-reduce of the residual over the leaf's members,
   `misc_computeIndexedSufficientStatisticsFast(treeY, indices+begin, len, ...)`,
   producing (sumWeights, sumWeightedResponse, sumWeightedResponseSq) cached on
   the Node (tree.hpp:172-177). This is E2's dominant "whole-tree suffstat
   recompute" and x86's #1 hotspot (32%).
3. metropolisJumpForTree (moves.hpp): a birth partitions one leaf's members
   (partitionChildren -> misc_partitionIndices, tree.hpp:621-676) and computes 2
   child suffstats; death merges incrementally (orphanChildren, tree.hpp:730-740);
   change/swap snapshot the affected index segment (snapshotSubtree,
   tree.hpp:765-773), refresh the subtree (repartition + recompute stats,
   refreshSubtree tree.hpp:695-704), and restore on rejection (restoreSubtree,
   tree.hpp:775-780). Scoring uses logLikelihoodForBranch -> the leaf's
   logIntegratedLikelihoodForNode over the cached node suffstat.
4. sampleParametersAndSetFits (chain.hpp:2076-2179, constant-leaf branch
   :2105-2136): per leaf draw mu ~ N(posteriorMean, posteriorSd^2) from
   (sumWeights, sumWeightedResponse) and SCATTER-write it into treeFits with
   `misc_setIndexedVectorToConstant` (x86's #3 hotspot, 15%).

The index buffer is `indexBuffer[t*n .. ]`, a persistent per-tree permutation
P_t; node [begin,end) is a contiguous range of P_t whose entries are scattered
observation ids (tree.hpp:167-206). Every hot pass above is memory-bound on this
shuffled buffer; both SIMD investigations concluded width buys ~0-1.15x and the
only escape is structural (make the field per-atom, not per-observation).

The constant-leaf math that makes atoms work (model.hpp:97-163):

- drawFromPosterior needs ONLY (sumWeights, sumWeightedResponse) and sigma^2/k
  (model.hpp:128-137). sumWeightedResponseSq is NOT used in the draw.
- logIntegratedLikelihood uses sumWeightedResponseSq only inside
  centeredSumOfSquares = sumWZSq - sumWZ*mean (model.hpp:117-118). This term is
  additive across any partition of a fixed observation set, so it CANCELS in
  every move ratio exp(newLogL - oldLogL) (fact 1.2; birth splits a leaf, the
  raw sumWZSq of parent = sum over children). It is dead weight for the decision
  (the 3.7 side finding); b>1 drops it.

## 2. Atom representation and the map

### 2.1 Definitions

Fix a block B = {t_0, ..., t_{b-1}} of b consecutive trees. Let

- O_i = sum_{k not in B} f_k(i)  -- the frozen OUTSIDE-block forest fit; constant
  for the duration of B's sub-sweep (the other m-b trees are not touched).
- g_i = w_i (y_i - O_i)  -- the block-static per-observation field: y and weights
  are fixed within a sweep, O is frozen within the block, so g is CONSTANT for
  the whole block sub-sweep. (Unweighted default: w_i = 1, g_i = y_i - O_i.)

An ATOM c is an equivalence class of observations under "same leaf in every tree
of B". Equivalently, c is a cell of the joint partition
leaf_{t_0} x leaf_{t_1} x ... x leaf_{t_{b-1}}. The atom map is the surjection
obs -> atom-id, plus each atom's b-tuple of leaves.

### 2.2 Per-atom state

For each atom c the map holds:

STATIC over the block (residual-independent, the cached "structure"):
- A(c)  = sum_{i in c} w_i               -- weight mass (count when unweighted)
- G(c)  = sum_{i in c} g_i = sum_{i in c} w_i (y_i - O_i)  -- static residual mass
- the member list: a contiguous slice of an atom-ordered index buffer
  `atomMembers[c.begin .. c.end)` (obs ids), and the b-tuple leaf_{t_j}(c).

DYNAMIC over the block (residual-dependent, the "carried" part):
- S(c)  = sum_{j in B} m_j(c)            -- the atom's total IN-BLOCK fit, where
  m_j(c) = mu_{t_j}(leaf_{t_j}(c)) is tree t_j's constant fit on atom c. Updated
  in O(atoms) whenever a leaf mean is redrawn.

The current full-forest residual mass on the atom is then DERIVED, never stored
redundantly:
   r(c) = G(c) - A(c) * S(c) = sum_{i in c} w_i (y_i - F_i),  F_i = O_i + in-block.

Note fact 1.2: NO per-atom sum-w-y^2 is needed for the constant leaf (sumWZSq
cancels). The frontier note's "sum w*y^2" is only for leaf models whose marginal
is quadratic in the response (linear/GP; section 6). The constant-leaf atom is
just (A, G, S, members, leaf-tuple).

### 2.3 The map data structure

Two coupled arrays, sized to atoms not n where possible:

- `atomMembers` : size_t[n]. A permutation of 0..n-1 grouped by atom: atom c owns
  the contiguous slice [atomBegin[c], atomEnd[c]). This is the atom-map analogue
  of the per-tree index buffer, but ONE buffer shared by the whole block (not b
  of them). Within an atom the members are stored in a fixed order (section 4
  fixes it for the b=1 bitwise anchor).
- `atomOf` : uint32_t[n] (or atom-id in obs order) -- the obs->atom map, used
  when a move needs to re-derive which atom an observation is in. Optional; can
  be recomputed from the slices. Keep it for O(1) move bookkeeping.
- per-atom SoA: A[c], G[c], S[c] (double x3), and leafTuple[c] (b x int32, the
  atom's leaf id in each block tree), atomBegin/atomEnd (size_t x2).

Atom count K is bounded by E1: K <= ~1.6e3 at b=8, n=1e5 (occupancy ~89), and
~10^2 at b=4. So the per-atom SoA is ~(3 doubles + b int32 + 2 size_t) * K ~
tens of KB at b=8 -- three orders below the n*m treeFits slab it displaces. The
n-sized arrays (atomMembers, atomOf, and the block-static g/w) are O(n), NOT
O(n*m): ONE atomMembers for the block, versus b per-tree index buffers today.

### 2.4 Footprint accounting (vs. what it replaces)

Per block, added: atomMembers (n size_t) + atomOf (n uint32) + g field (n double)
+ per-atom SoA (~K*(3*8 + 4b + 16) bytes). At n=1e5, b=8: ~0.8 MB + 0.4 MB +
0.8 MB + ~0.1 MB ~ 2.1 MB. This is REUSED across blocks (one workspace, sized
once), so total added memory is ~2 MB, independent of m. It does NOT add an
n*m array (contrast data-layout.md strategy 1a's rejected leaf-order residual
copy at n*m doubles = 160 MB). The persistent treeFits (n*m) and per-tree index
buffers (n*m size_t) stay as they are (needed at block boundaries and for
state/predict). Net: O(n) + O(K) added, O(1) in m.

## 3. The O(atoms) reformulation for the constant-leaf Gaussian (correctness core)

Claim: every in-block operation is exact atom algebra over (A, G, S) and the
O(b) in-block leaf means, with NO O(n) pass, and it equals today's per-observation
computation (exactly at b=1; same distribution, regrouped FP at b>1).

Notation: tree t_j in the block, its leaf L, atoms mapped into L are
`atoms(L) = { c : leaf_{t_j}(c) = L }`. Define S_{-j}(c) = S(c) - m_j(c) =
sum_{s in B, s != j} m_s(c), the in-block fit from every block tree EXCEPT j.

### 3.1 Leaf sufficient statistic (replaces setNodeAverages / computeLeafStats)

The residual entering tree t_j is resid_i = y_i - O_i - sum_{s != j} f_s(i), i.e.
the full residual with tree t_j's own fit added back (exactly what treeY holds
when setNodeAverages runs today). For a leaf L of tree t_j:

  sumWeights(L)          = sum_{c in atoms(L)} A(c)                        (= W(L))
  sumWeightedResponse(L) = sum_{c in atoms(L)} [ G(c) - A(c) * S_{-j}(c) ]

Because all c in atoms(L) share tree t_j's leaf L, m_j(c) = mu_{t_j}(L) is one
constant, so S_{-j}(c) = S(c) - mu_{t_j}(L) and the sum factors:

  sumWeightedResponse(L) = sum_{c in atoms(L)} [ G(c) - A(c)(S(c) - mu_{t_j}(L)) ]
                         = sum_{c in atoms(L)} [ G(c) - A(c) S(c) ]
                           + mu_{t_j}(L) * sum_{c in atoms(L)} A(c),

equivalently

  sumWeightedResponse(L) = D(L) + mu_{t_j}(L) * W(L),
     D(L) = sum_{c in atoms(L)} r(c),  r(c) = G(c) - A(c) S(c).

This is fact 1.2 verbatim, per atom. The cost is one pass over the atoms grouped
by their tree-t_j leaf = O(K) for the whole tree, replacing O(n) gather. Feed
(W(L), sumWeightedResponse(L)) into the UNCHANGED model.hpp math: same
logIntegratedLikelihood and drawFromPosterior, so the leaf draw and the move
likelihood are byte-for-byte the same functions of the suffstat. sumWZSq is
dropped (b>1) or supplied per-atom (b=1 anchor, section 5).

### 3.2 Leaf-mean draw (replaces the constant-leaf branch of sampleParametersAndSetFits)

drawFromPosterior(rng, k, W(L), sumWeightedResponse(L), sigma^2) unchanged
(model.hpp:128-137). Draw one mu per leaf of tree t_j. RNG consumption identical
to today (one standard normal per non-empty leaf, in the SAME leaf order -- see
section 4.2 on fixing the leaf-iteration order).

### 3.3 The residual roll (replaces the fused O(n) treeY update)

Today: after tree t_j's draw, treeY is rolled to tree t_{j+1} by an O(n) pass.
In atom form, the draw changed tree t_j's fit from mu_old(L) to mu_new(L) on each
leaf; the roll is:

  for each leaf L of tree t_j, for each atom c in atoms(L):
     delta = mu_new(L) - mu_old(L)
     S(c)  += delta          // m_j(c) changed by delta; S = sum_s m_s

That is O(K) total (each atom is in exactly one leaf of tree t_j, so the double
loop visits each atom once). No treeY touched. When the next tree t_{j+1} is
scored, section 3.1 reads S(c) and subtracts m_{j+1}(c) -- it automatically sees
tree t_j's NEW fit. This is the Gauss-Seidel coupling, carried as S per atom
instead of a rolled per-observation field.

### 3.4 Rejected-move restore (replaces the index-segment snapshot/restore)

A rejected birth/change/swap must leave (A, G, S, members, leaf-tuple) exactly as
before the proposal. Because structure proposals do not change any drawn mu
(fact 1.1: structure is residual-independent, the leaf means are redrawn AFTER
the move resolves), the ONLY atom state a rejected move can dirty is the atom
partition (splits/merges) and the derived leaf suffstats -- never S or G of an
UNsplit atom. So restore = undo the atom split/merge (section 4.3), an O(touched
atoms) operation, mirroring restoreSubtree's O(index-segment) cost. No O(n) pass.

### 3.5 Fit write-back at block exit (the one O(n) scatter per block)

At block exit, the b trees' drawn leaf means must be materialized into the
per-observation treeFits[t_j * n + i] (obs order) so that: (a) totalFits and the
residual are correct for the sigma/latent draws (chain.hpp:790-802), and (b) the
next block reads a correct O for its own g field. For each block tree t_j and
each leaf L, scatter mu_{t_j}(L) to the leaf's member observations. Two equivalent
routes:
  - per-atom: for each atom c, for each block tree j, write mu_{t_j}(leaf_j(c))
    to treeFits[t_j*n + i] for i in c's members -- one walk of atomMembers writes
    all b trees' fits (b sequential writes per member, cache-friendly).
  - or reuse the existing per-tree setTreeFitsFromParameters (chain.hpp:2035) b
    times.
This is O(bn) writes = O(n) per tree amortized, but done ONCE at block exit, not
interleaved with b residual rolls + b suffstat gathers. That is the DRAM
collapse: E2's field passes (roll + suffstat + scatter, ~85%) go from b per block
to ~1-2 per block.

### 3.6 Why this equals today's computation

- b=1: atoms == leaves of the single tree, S_{-j} = 0 (no other block tree), so
  sumWeightedResponse(L) = G(L) = sum_{i in L} g_i, and g_i = treeY_i. If the atom
  reduces g over its members in index-buffer order with the same kernel, this is
  BITWISE misc_computeIndexedSufficientStatisticsFast(treeY,...). Section 4.4.
- b>1: the value sum_{c in atoms(L)}[G(c) - A(c)S(c)] + mu W(L) equals the true
  sum_{i in L} w_i resid_i in EXACT arithmetic (G, A, S are exact partial sums of
  the same terms). In floating point the grouping differs (atom partials vs a
  single index-order reduction), so draws SHIFT -- same posterior, different bytes.
  The affine identity is not an approximation; it is the exact regrouping the
  frontier note calls "the identical Gauss-Seidel sweep, arithmetic regrouped".

### 3.7 Sigma draw and latents stay O(n), once per sweep

The sigma draw's SSR (chain.hpp:908, misc_computeSumOfSquaredResiduals over the
full residual) and the latent refresh (refreshLatents) are OUTSIDE the block
machinery: they run once per sweep after all blocks, over the correct
per-observation residual materialized at block exits. They are unchanged and
remain O(n) once/sweep (already ~1/m the suffstat frequency; not a target).

## 4. Atom-map maintenance under MH moves and across sweeps

This is the hard part. Moves mutate a block tree mid-sub-sweep and slice atoms.

### 4.1 The move taxonomy in atom terms

- BIRTH in tree t_j: leaf L splits by a new cut on variable v into L_left, L_right.
  Each atom c in atoms(L) is partitioned by the cut: members with code < cut go
  to a child atom c_left (leaf-tuple identical to c except tree t_j's slot ->
  L_left), the rest to c_right. Atoms whose members all fall one side are just
  relabeled (leaf slot updated, no split). This is the ~7% raw-data touch: it
  reads column v's codes for L's members and re-partitions them. Cost O(n_L).
- DEATH in tree t_j: L_left, L_right merge into parent L. Their atoms merge
  pairwise: a c_left and c_right that agree on all OTHER block trees' leaves
  recombine into one atom (A, G add; S recomputes with the merged leaf's future
  mu). Because death is the exact inverse of a birth, and E1 atoms are few, this
  is cheap; in practice death reuses the pre-split atom identity if it is still
  recorded (section 4.3).
- CHANGE / SWAP in tree t_j: the node's rule (and a subtree's routing) changes,
  re-partitioning a subtree near the root. Every atom mapped under that subtree
  can be re-sliced. E2's correction: shallow-tree change/swap re-partition near
  the root, so this touches ~2x the modelled 7% (hence the realized ~6x not
  7-10x). In atom terms: snapshot the affected atoms, rebuild the joint partition
  restricted to the changed subtree's observations by re-routing them through the
  b trees, O(n_subtree). This is the surviving scan the design does NOT remove
  (and should not pretend to).

### 4.2 Splitting an atom exactly (the birth kernel)

To split atom c by tree t_j's new cut on v:
  1. Partition c's member slice `atomMembers[c.begin..c.end)` in place by v's
     code (two-pointer, EXACTLY misc_partitionIndices / the tree.hpp two-pointer
     -- an integer-exact permutation, bitwise across ISAs; section 4.4 keeps the
     order matching the tree's own partition so member order stays canonical).
  2. Accumulate the two children's A and G by reducing w_i and g_i over each
     child's sub-slice. THIS is where the block-static g/w fields are read (the
     only per-observation read a move needs). Order-fixed reduction (section 4.4).
  3. Children inherit S_{-j}(c) (all OTHER block trees' leaves unchanged), so
     S(c_left) = S(c_right) = S(c) initially -- but tree t_j's leaves L_left,
     L_right have no drawn mu yet during scoring (structure first, params after).
     The scoring suffstat adds back tree t_j's fit via mu_{t_j}(L) which is being
     split away, so the two child leaf suffstats are computed from
     (A(child), G(child) - A(child) S_{-j}(c)) directly -- consistent with 3.1.
  4. Register the child atoms (new ids or recycled), update atomOf for the moved
     members, update leafTuple.

The leaf-iteration order and the member order within an atom are FIXED (canonical
= the order the tree's own partition would produce), so at b=1 the atom reduction
visits observations in the same order as today's index buffer, and RNG draws
(one normal per leaf) occur in the same leaf order. This is what makes the b=1
anchor bitwise (section 4.4) and keeps b>1 machine-independent (a fixed atom
order reproduces the same regrouped sum on every ISA).

### 4.3 Rollback (rejected moves)

Mirror the tree's existing snapshot/restore (tree.hpp:758-780):
- BIRTH rejected: discard the two child atoms, restore the parent atom's slice
  order (the two-pointer partition is unstable, so snapshot the affected
  atomMembers segment before partitioning, exactly as snapshotSubtree snapshots
  the index segment), and restore leafTuple/atomOf. O(touched atoms + n_L).
- CHANGE/SWAP rejected: snapshot the affected atoms' (A, G, S, members,
  leaf-tuple) and the atomMembers segment under the changed subtree before the
  re-partition; restore on rejection. This is the atom analogue of
  SubtreeSnapshot; its cost is the same O(n_subtree) the engine already pays.

Because a rejected structure move never consumed a leaf-mean draw (params are
drawn only after the move resolves, once per accepted structure), rollback does
NOT touch RNG and does not perturb S/G of unsplit atoms -- consistent with 3.4.

### 4.4 The b=1 bitwise guarantee, concretely

At b=1 the block is one tree; the atom map is that tree's leaf partition, and
atomMembers for the block IS the tree's own index buffer slice per leaf. Set
g_i = treeY_i (reuse the exact residual the roll already maintains -- at b=1 the
block-entry g build is the current per-tree roll, unchanged and still O(n)). Then:
  - suffstat: reduce g over each leaf's members with the SAME kernel and the SAME
    index order = misc_computeIndexedSufficientStatisticsFast on treeY. Bitwise.
  - sumWZSq: the b=1 path still computes it per leaf over g (same kernel), fed to
    the unchanged logIntegratedLikelihood. Bitwise.
  - draw: same (W, sumWResid, sumWZSq) -> same posterior -> same RNG -> same mu.
  - roll: at b=1, block exit writes tree t's fit and the next block (next tree)
    rebuilds g from the rolled residual -- identical to today's per-tree roll.
  - moves: the atom split/partition IS the tree's partitionChildren (same
    two-pointer, integer-exact); rollback IS snapshotSubtree/restoreSubtree.
So the b=1 path is a pure structural refactor of today's per-tree loop expressed
in atom vocabulary. It exercises the map build, per-leaf atom aggregation, the
draw, the roll, splitting, and rollback -- everything except the b>1 affine
regrouping -- and MUST match the equivalence baseline byte-for-byte. This is the
primary correctness anchor and the milestone gate (section 7, stage A).

### 4.5 Persistence across sweeps (O(changed leaves))

The joint partition of a block changes between sweeps only by the moves accepted
in that block's b trees last sweep -- O(1) accepted moves per tree per sweep
(shallow trees, low accept rates; E3 measured ~7-12% accept). So the atom map is
NOT rebuilt each sweep; it is PATCHED as moves are accepted during the sub-sweep
(4.2/4.3 already maintain it live) and persists into the next sweep intact. A
full rebuild (re-route all n observations through the b trees, scatter into atom
buckets, O(bn)) is needed only at block initialization, at a block-composition
change, or as an optional periodic re-canonicalization to bound fragmentation
(atoms can accumulate near-empty slivers after many splits/merges; a rebuild
every K sweeps re-tightens occupancy -- tunable, measure if needed).

This generalizes the landed U'WU cache (model.hpp:458-539): there, each leaf's
residual-INDEPENDENT crossproduct is cached and re-validated by comparing the
ordered member list against tree.indices[begin..end); a structural move that
alters membership fails the compare and rebuilds. Here the residual-INDEPENDENT
structure (the atom partition + A, G) is cached and patched on membership change,
while the residual-DEPENDENT part (S, hence r) is carried incrementally. Same
template, one axis up (joint over b trees instead of per leaf).

## 5. Relation to the contiguous-layout substrate, and scope

### 5.1 Does 3.4 subsume, need, or ignore data-layout.md?

3.4 SUBSUMES the single-tree contiguous reorder and does NOT need it landed
first. data-layout.md's own verdict (section 8.1, open Q1): treat the contiguous
layout as the SUBSTRATE for 3.4, justify it by 3.4's payoff, and "do NOT ship the
single-tree reorder as a standalone perf item unless 3.4 is deferred." The reason
they are one program:

- The single-tree reorder gathers the residual into leaf-contiguous scratch,
  runs the suffstat/draw sequentially, and scatters fits back -- paying one O(n)
  gather + one O(n) scatter PER TREE for a ~10% win, because the residual's
  cross-tree coupling forbids carrying leaf order across trees (data-layout.md
  section 3). 3.4 is exactly that primitive with the gather/scatter amortized
  across b trees: build g once per BLOCK, carry the coupling as S per atom, pay
  the scatter once per block. Same substrate (contiguous per-cell fields, fixed
  reduction order), the block just makes the amortization win multi-x instead of
  ~10%.
- The atom-ordered atomMembers buffer + block-static g/w fields ARE the
  contiguous-per-node layout, generalized from per-tree to per-block. So building
  3.4 builds the substrate; there is no separate layout milestone to land first.

Recommendation: do NOT land the standalone single-tree reorder. Go straight to
atoms; the b=1 anchor (section 4.4) delivers the layout plumbing as a bitwise
refactor, and b>1 delivers the win the standalone reorder could not.

### 5.2 Interaction with the data-ownership partition rework

data-layout.md flags one shared touchpoint: partitionChildren (tree.hpp:621),
which the data-ownership program is already reworking (u8/u16 width templating).
3.4's atom split (section 4.2) partitions atomMembers by a column's codes -- it
REUSES the same two-pointer/misc_partitionIndices primitive, so it inherits the
width-templated partition rather than forking it. 3.4 should FOLLOW the
data-ownership partition settling (as the substrate note recommends), then call
the settled partition on the atom buffer. No conflict: atoms add a second buffer
the same partition primitive runs over.

### 5.3 Scope: constant-leaf Gaussian first

v1 is the CONSTANT-LEAF GAUSSIAN only -- the ~85% field case E2 measured, the
default sampler, and the case where fact 1.2 collapses the atom aggregate to just
(A, G, S). This is the whole win for stock `bart()`/`dbarts()` and stan4bart's mu
forest. Generalization order, template = the U'WU cache:
- LINEAR / GP leaves: their marginal is quadratic in the response (U'Wz, z'Wz),
  so an atom needs per-atom U'Wg and the U'WU block -- the atom is the joint of
  the leaf's own crossproduct cache across b trees. Heavier per-atom state (p^2),
  fewer atoms help less (these leaves are deeper). DEFER; the U'WU cache already
  gives them most of the residual-independent reuse per leaf. Design later.
- WEIGHTS: fully supported in v1 -- g_i = w_i(y_i - O_i) and A(c)=sum w_i carry
  weights exactly; the atom math is the weighted suffstat. Latent families whose
  weights vary per sweep (workingWeightsVaryPerSweep, chain.hpp:797) force a
  g/A rebuild each sweep, which the block already does at entry.
- PROBIT / LOGISTIC latents: the latent refresh (refreshLatents) rewrites y (and
  weights) once per sweep, OUTSIDE the blocks. The block just rebuilds g at entry
  from the current working response -- no special handling; v1 covers binary BART.
- GROUPED RANDOM EFFECTS / BCF: each forest sub-sweeps its own residual
  (formForestResponse, chain.hpp:700-704); blocks are per-forest, g uses the
  forest's response net of the other forest's scaled contribution. The glue/ridge
  interweave (chain.hpp:804-807) runs once per sweep outside blocks. Compatible;
  block per forest.
- MISSINGNESS (MIA): the partition already routes missing via the rule's
  missing-direction bit (partitionIndicesMIA in the worktree); atom split reuses
  it. No new handling.
- GROW-FROM-ROOT warm start (chain.hpp:938-1048): a separate init path; leave it
  on the per-tree code path (it runs a fixed few sweeps at init, not hot). Blocks
  engage only in the steady-state run() loop.

## 6. Interaction with within-chain threading and a future GPU seam (noted, not designed)

- Within-chain threading (docs/plans/within-chain-threading.md) parallelizes the
  O(n) passes with a FIXED-BLOCK reduction for thread-count invariance. 3.4
  COMPOSES: it makes those O(n) passes rarer (once per block, not per tree), so
  the barrier count drops from ~3m to ~3(m/b) per sweep, AND the surviving passes
  (g build, fit scatter-back) are contiguous atom-ordered buffers that partition
  cleanly into fixed blocks with no false sharing (data-layout.md section 7). The
  O(atoms) interior is serial and cheap (K ~ 10^3), so it need not thread. Both
  changes are shifting-class and should share ONE re-record if landed together.
  Do NOT double-count: the g build and scatter-back are themselves O(n) passes in
  the threading budget.
- GPU seam (frontier section 4): 3.4 reduces host<->device command traffic --
  the device holds x, the sorted orders, and the atom workspace resident; the host
  drives b trees of O(atoms) logic through a command stream, and only the g build
  / fit scatter cross the bus per block. The atom-ordered contiguous buffers are
  exactly the resident shape a seam wants. Design the seam only after the CPU
  experiments pick winners (frontier's explicit sequencing); note only.

## 7. Staged implementation plan

Big change; stage it so each stage has an independent gate and a clean abort.
The gate procedure is the standard one (CLAUDE.local.md): R CMD INSTALL
(--preclean after headers), tinytest, tests/cpp component tests, equivalence
(bit-identical where claimed / statistical + exact-posterior where re-recorded),
and bench-sampler on the QUIET machine as the headline gate (the whole point is
DRAM/throughput). Budgets are rough engineering estimates.

### Stage A -- b=1 exact refactor (the atom machinery, BIT-IDENTICAL). ~2-3 wk.

Build the atom map, the per-leaf atom aggregation, the atom split/rollback, and
the block driver, wired at b=1 so atoms == leaves and g == treeY (section 4.4).
Deliverable: run() routes the constant-leaf sweep through the atom path at b=1.
GATES (all MUST hold, NO re-record):
- equivalence bit-identical vs the current baseline (this is the whole point of
  the stage -- it proves the machinery is a faithful refactor);
- tinytest full pass with NO snapshot regeneration (bitwise draws unchanged);
- tests/cpp: new component tests for atom build / split / rollback / persistence
  (oracle = the live per-tree partition; the sampler-internals testing lessons in
  MEMORY apply -- snapshot-vs-live oracle, populate state for copies);
- bench-sampler: expect NEUTRAL (b=1 does the same O(n) work); a regression here
  means the atom bookkeeping overhead is too high -- a kill signal for the whole
  approach, caught cheaply before any FP change. This is the key de-risk.

### Stage B -- b>1 fusion for the constant-leaf Gaussian (the WIN, re-record). ~3-4 wk.

Turn on b>1: build g once per block from y - O, carry S per atom, regroup the
suffstat/roll/draw to O(atoms) (section 3), keep the ~7% move scans. Drop sumWZSq
from the constant-leaf atom (fact 1.2). Choose the fixed atom-iteration order and
the fixed reduction order so draws are machine-independent across ISAs.
GATES:
- equivalence: NOT bit-identical (expected). Gate STATISTICALLY (matched-RNG
  distribution equality vs baseline) PLUS exact-posterior checks (SBC /
  known-posterior scenarios already in the gate suite) -- the shifting class,
  same discipline as the suffstat/threading re-records;
- machine independence: the fixed atom order must reproduce the SAME bytes on
  scalar/NEON/SSE2/AVX2 -- add a cross-ISA bitwise component test at b>1 (the
  invariant simd-survey.md protects);
- tinytest: regenerate the RNG-locked snapshots ONCE (replay whole files);
- bench-sampler on the quiet machine at n in {1e4, 1e5, 1e6}, m=75/200, b in
  {4,8}: HEADLINE gate. Target the E1/E2-modelled ~6x DRAM drop -> a large
  wall-clock win at n >= 1e5. Kill/scope-back if the realized speedup is far
  below the model (e.g. move-scan re-partition cost dominates at shallow trees).
- This stage carries the single, VD-approved one-time re-record.

### Stage C -- tune b and harden. ~1-2 wk.

Sweep b (auto-pick by n, m, tree depth per E1's occupancy curve: larger n favors
larger b since atoms/n falls; deep trees push atoms up so cap b lower), decide
the periodic-rebuild cadence, tune the cache/rebuild thresholds, and settle the
default (likely b=4 or 8, per E1). GATES: bench-sampler across the grid; confirm
no regression at small n (blocks should fall back to b=1 or skip below a cutoff,
like the SIMD toggles); tinytest/equivalence stable (b is a perf knob, not a
draw change once the fixed order is set -- but VERIFY b-invariance is only
distributional, or fix ONE default b in the equivalence anchor).

Later (separate programs, not this plan): linear/GP atoms (section 5.3),
within-chain threading composition (section 6), GPU seam.

## 8. Risks and open questions for VD

### 8.1 Risks

- Move-scan re-partition cost. E2's correction: shallow-tree change/swap
  re-partition near the root, ~2x the modelled 7% scan, which is why the realized
  drop is ~6x not 7-10x. If atoms fragment or change/swap dominate, the surviving
  O(n_subtree) scans could eat more of the win. Mitigation: the bench gate at
  stage B measures it directly; the atom rebuild cadence (stage C) bounds
  fragmentation. This is the main way the win underperforms.
- Complexity / bug surface. The atom map is a second partition structure with its
  own split/merge/snapshot/rollback, maintained in lockstep with b trees. The
  b=1 anchor (stage A) de-risks the plumbing bitwise before any FP change, and
  the rollback mirrors the existing SubtreeSnapshot exactly. Fuzz the map against
  the live per-tree partition (mutation-fuzzing.md precedent).
- The re-record. One-time, shifting class, machine independence preserved by the
  fixed atom order -- the ask VD already accepted. Snapshots regenerate, the
  equivalence anchor re-records at a fixed default b, z-mode passes.
- Memory. O(n) + O(K) workspace, reused across blocks, O(1) in m (section 2.4) --
  small. No n*m addition (unlike data-layout.md 1a).
- Bench noise. bench-sampler must run on the quiet machine, never concurrent with
  other load (MEMORY: orchestrator pacing).

### 8.2 Open questions for VD

1. DEFAULT b. E1 says atoms stay 2-4 orders below n at b in {4,8} and atoms/n
   falls with n. b=8 gives a bigger amortization but ~15x more atoms than b=4
   (still tiny). Recommendation: default b=4 conservative, or auto-pick b=8 for
   n >= 1e5. VD's call on whether b is a fixed default or n-adaptive (n-adaptive
   complicates the equivalence anchor -- see Q4).
2. SCOPE of response families in v1. Recommendation: constant-leaf Gaussian only
   (covers stock BART, binary via latents, weights, BCF/grouped-RE forests);
   defer linear/GP (their U'WU cache already gives per-leaf reuse). Confirm v1 is
   constant-leaf only.
3. RE-RECORD timing. Stage A is bit-identical (no re-record); stage B carries the
   one re-record. Confirm bundling it with within-chain threading's re-record IF
   both land in the same window (one snapshot regeneration instead of two), or
   keep them independent.
4. EQUIVALENCE ANCHOR under a b knob. If b is n-adaptive, the equivalence gate
   must fix ONE b (or assert only distributional equality across b). Recommend
   fixing b in the anchor scenario and gating b-invariance statistically. VD to
   confirm the gate policy.
5. GO STRAIGHT TO ATOMS vs build the contiguous layout first. Recommendation (and
   data-layout.md's): skip the standalone single-tree reorder; the b=1 anchor
   delivers the layout as a bitwise refactor and b>1 delivers the real win.
   Confirm.

## 9. Stage A landing note (2026-07-11, LANDED, bit-identical)

Stage A is complete and DEFAULT ON. The constant-leaf steady-state run() sweep
sources its per-leaf sufficient statistics through the b=1 atom path
(src/bartcore/atoms.hpp, wired in chain.hpp), and birth/death/change/swap are
atom-routed (moves.hpp). It is a pure structural refactor: the atom aggregation
is bitwise the legacy misc kernel over the same member order, landing the same
node-cache values, so every downstream reader and every draw is unchanged.

WHAT SHIPPED

- The b=1 atom path is the DEFAULT (`#define BARTCORE_BLOCK_FUSION 1`). The knob
  stays: a build defines `-DBARTCORE_BLOCK_FUSION=0` for a one-line abort back to
  the legacy setNodeAverages/computeLeafStats writer, and it is Stage B's
  b-dispatch seam. Both writers stay compiled so tests/cpp drives them side by
  side.
- SCOPE is constant-leaf + steady-state only, provably: `useAtomSuffstatSource`
  gates on `Chain::leafIsConstant` (false for vector-param linear and
  function-param GP leaves - they keep the legacy writer), and the atom path is
  wired only in run()'s steady-state loop, never the grow-from-root warm start
  (which calls setNodeAverages directly). Categorical / sparse / MIA / weighted /
  probit-logistic latents / BCF-grouped are covered (the aggregation just calls
  the same kernel over the same members; g rebuilds from the working response
  each sweep; the A cache drops on any working-weight change).
- DESIGN A (section 2.4, 4.4) is the resolved implementer decision: `members`
  ALIASES the block tree's `tree.indices` at b=1; no owned second n-buffer, no
  member-copy pass. The owned per-block buffer arrives with Stage B.

b=1 OVERHEAD REDUCTION (the milestone's precondition; commit "Reduce the b=1
atom path to a single per-tree pass"). At b=1 nothing amortizes, so the atom
path had to match legacy's O(n) work exactly. Three b=1-only overheads were
removed, all bitwise:
- the cross-sweep A cache (section 4.5) SAVES NO WORK at b=1 - the monolithic
  suffstat kernel rescans G/Q every sweep and drags A along - and only ADDS a
  per-leaf member-list memcmp. It is gated OFF at b=1 (`aCacheBypass`) and stays
  a STAGE-B mechanism (kept default-on for its component tests).
- the obs->atom map `atomOf` is Stage-B move bookkeeping no b=1 consumer reads;
  its O(n) per-tree scatter (and the moves' O(n_leaf) writes) are gated OFF at
  b=1 (`trackAtomOf`), so buildForTree is O(#leaves).
- build + aggregate + write walked the tree THREE times where legacy walks once;
  fused into `AtomMap::buildAggregateWrite` - one fillBottom traversal that
  records topology, aggregates through the kernel, and writes the node cache
  inline. The inert S(c)=mu carry is dropped from the hot loop (kept in
  setInBlockFits for the harness + Stage B). Net per-tree cost == setNodeAverages'
  single traversal + O(#leaves) SoA bookkeeping.

GATES (all held, NO re-record - Stage A is bit-identical by construction):
- equivalence vs benchmarks/baselines/equivalence-ac6ec2c.rds: 22/22 IDENTICAL
  (same RNG stream), exact.
- tinytest: full pass (2728 ok, 0 fail) from a --preclean default-ON install, NO
  snapshot regeneration.
- tests/cpp: clean rebuild all pass + the full-vocabulary mutation fuzzer + the
  atom build/aggregate/split/merge/refresh/rollback/A-cache suites.
- bench-sampler on the x86 box (dbarts-bench), interleaved A/B pinned to one core
  (medians of 5 paired passes), atom default-ON vs legacy (-DBARTCORE_BLOCK_FUSION
  =0) built from the same tree: every run-* throughput metric within 4% (n=1000
  t75 1.039, t200 1.037, binary 1.032; n=10000 1.001), embedded-offset 1.024,
  setPredictor 0.997/1.017 - NO metric over the 5% gate. NEUTRAL, as the b=1
  anchor predicted (the residual few-% at small n is the irreducible atom
  bookkeeping b=1 pays with no amortization; it vanishes as n grows - n=10000 is
  essentially exact). This de-risks the machinery for Stage B, where b>1 turns
  the same passes into the ~6x DRAM win.

NEXT: Stage B (b>1, the win, the one approved re-record) builds on this map -
the flag becomes the b-dispatch, the A cache and atomOf come back on, `members`
becomes an owned per-block buffer, and the suffstat/roll/draw regroup to
O(atoms).

## 10. Stage B landing note (2026-07-13, LANDED dormant, bench VERDICT: KILL)

Stage B's machinery (plan commits (i)-(iv-a) plus follow-ons,
docs/plans/block-fusion-stage-b.md) is landed behind a runtime blockSize knob.
The shipped default stays 1: the fused b>1 path never engages in production,
Stage A's bitwise b=1 anchor is untouched, and DBARTS_BLOCKSIZE forces b>1 only
for gates and benching.

CORRECTNESS GATES now in place at b>1 (all held): the joint atom-map mutation
fuzzer (patched map matches a from-scratch buildForBlock rebuild after every
move and on rejection, bitwise), plus a run()-path undo fuzzer added with the
move-maintenance rewrite (catches a rejected change/swap leaving a cache leak
outside the touched subtree); affine identity under moves at widths 1-3 (max
|diff| ~4e-13); cross-ISA bitwise (scalar/SSE2/AVX2/NEON) at b>1; the forced
b>1 exact-posterior anchors at small n (via DBARTS_BLOCKSIZE, bypassing the
default so the anchors actually exercise the affine path, the plan's
FORCE-b>1-ON-THE-ANCHORS precondition).

HEADLINE BENCH VERDICT: KILL. dbarts-bench (x86 AVX2, 16 cores, single-chain
single-thread), friedman p=10, m=75, screening protocol (reduced reps of the
biggrid harness; the ~7x margin needs no finer resolution); two recordings
bracket the move-maintenance fix (b=1 essentially unchanged, ~29.1/29.4
ms/iter at n=1e5 and ~498/505 at n=1e6 across both):
- before (tip b74fec9, O(n)-per-move regroup): n=1e5 b=4 488ms (16.8x
  slower), b=8 749ms (25.7x); n=1e6 b=4 7045ms (14.1x), b=8 9945ms (20.0x).
- after (tip 10c59bc, in-place O(affected) move maintenance + targeted
  rejection undo): n=1e5 b=4 229ms (7.8x slower), b=8 285ms (9.7x); n=1e6
  b=4 3482ms (6.9x), b=8 3694ms (7.3x).

The fix bought the modelled ~2.1-2.7x on the fused path but did not change the
sign: b>1 still loses ~7x at n=1e6, the slowdown ratio is roughly flat in n
(no crossover at any feasible size), and b=8 costs about what b=4 does.

ATTRIBUTION: buildForBlock rebuilds the joint atom map from scratch at every
block entry, every sweep -- O(bn) per block, O(mn) per sweep -- because
cross-sweep map persistence (section 4.5, rebuild cadence) was explicitly
deferred past Stage B. Second contributor: the design-admitted near-root
change/swap subtree re-slice (8.1).

PHASE DECOMPOSITION (2026-07-13): env-gated wall timers inside runFusedBlock
(diagnostic only, never committed) split the fused iteration into phases on
dbarts-bench, quiet, m=75, tip 31a36b9, the same screening protocol as above.
ms/iter and share of the fused iteration (timer sums match R-measured wall
within ~1% in every cell):
- n=1e5 b=4 (total 229.6): entry-ex-build 5.3 (2%), buildForBlock 105.1 (45%),
  affine cache writes 0.1, moves 89.2 (39%), draws 9.4 (4%), exit scatter 17.2
  (8%), sweep-end rebuild 4.9 (2%).
- n=1e5 b=8 (total 279.4): buildForBlock 109.2 (39%), moves 127.1 (46%),
  scatter 22.4 (8%).
- n=1e6 b=4 (total 3414.7): buildForBlock 1333.5 (39%), moves 1014.0 (30%),
  scatter 821.4 (24%), draws 135.2, entry 57.5, rebuild 47.4.
- n=1e6 b=8 (total 3691.4): buildForBlock 1341.2 (36%), moves 1295.3 (35%),
  scatter 814.0 (22%).

The decision quantity is the theoretical BEST CASE for persistence (4.5):
delete buildForBlock for free and compare what remains to the same-n b=1
sweep:
- n=1e5 b=4: 124.5 vs 28.6 -> 4.4x slower. b=8: 170.2 vs 28.6 -> 6.0x.
- n=1e6 b=4: 2081.2 vs 497.7 -> 4.2x. b=8: 2350.2 vs 497.7 -> 4.7x.

So persistence alone is not sufficient: buildForBlock is only 36-45% of the
fused iteration. The b>1 move phase alone (atom-map maintenance + affine
cache re-derivation inside metropolisJumpForTree) costs 2-2.6x the ENTIRE
b=1 sweep at n=1e6, and the block-exit scatter adds 22-24% there -- it
scales ~48x for a 10x increase in n, so the treeFits working set at n=1e6
(~600 MB, past any LLC) is DRAM-bound, and that cost is outside anything 4.5
touches.

TRAFFIC-MODEL CORRECTION: re-examining the cost model behind the original
~6x DRAM projection (section 0, section 7 Stage B) finds it over-counted the
amortizable traffic. The n*m treeFits slab must be READ once per sweep (to
form the frozen outside-block field O) and WRITTEN once per sweep (to record
the drawn fits) at ANY block width -- fusion reduces the number of passes,
but each fused pass is b times larger, so those bytes never shrink with b.
Only the running-total roll and the suffstat index-gather genuinely amortize
with b. The bandwidth ceiling attributable to blocking was therefore
~1.3-2x, not ~6x, before any implementation cost was counted. Separately, at
n=1e5 the dominant b=1 costs -- suffstat gather (x86's #1 hotspot, 32%,
section 1) and fit scatter (x86's #3 hotspot, 15%, section 1), together ~47%
of the sweep -- run through the shuffled per-tree index buffer (section 1)
and fit in LLC; they are latency-bound, not DRAM-streaming-bound, so a
bandwidth-amortization design does not address them either. The measured
b=1 throughput (~240 MB per sweep in ~29 ms at n=1e5) implies ~8 GB/s
effective single-core -- already close to the single-core streaming ceiling.

STAGE C: CLOSED, not conditional. Cross-sweep persistence (4.5) is necessary
but, per the phase decomposition above, measured insufficient: its best case
is 4-6x slower than b=1, not a win. The fused approach cannot win on this
architecture at any feasible n or b; Stage B's machinery stays in place,
dormant, behind blockSize, as an exact, gated implementation -- no default
flip, no baseline re-record, no MANIFEST change. The honest summary of the
block-fusion bet: mathematically valid (exact Gauss-Seidel regrouping,
verified Q-drop cancellation), engineering sound, but the underlying
DRAM-amortization premise over-counted the amortizable traffic by ~3x. Any
future revival of b>1 must first wire a forced-b>1 posterior anchor into the
standing CI (e.g. a DBARTS_BLOCKSIZE variant in the sbc workflow) -- the
categorical-exact / logistic-reference / BCF exact-posterior gates at forced
b>1 (section 7 of the plan) were one-time manual validations, and nothing in
CI forces b>1 today.

EXCISED (2026-07-13): the machinery is removed from the tip -- both the b>1
fused path and the Stage A b=1 atom writer (the whole AtomMap, its move kernels,
the blockSize knob, and the DBARTS_BLOCKSIZE override). The complete
implementation is preserved on branch archive/block-fusion at a225e75; the
rationale (KILL verdict, traffic-model correction) stays in sections 9-10 above.
The excision was verified bitwise -- the R equivalence gate reports identical
draws against the same anchor (equivalence-ac6ec2c.rds), because the Stage A
atom writer was bit-identical to the legacy setNodeAverages / computeLeafStats
path it routed around. The legacy writer is again the only suffstat path.
