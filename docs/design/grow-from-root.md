# Grow-from-root: cost and validity of root-down tree sampling

Status: memo, 2026-07-08. Decides whether dbarts adds XBART-style root-down
stochastic tree construction, and in what role. Grounded in a measured
cut-scan kernel (benchmarks/kernels/grow_from_root.c) and in the actual move
dispatch (moves.hpp metropolisJumpForTree, model.hpp CGMTreePrior/DartPrior).
The cut-scan histogram it needs is the same primitive the informed-proposal
construction of parallel-bart-frontier.md 3.1 wants, so the cost table below
reports the numbers that construction's falsifier needs too.

## 1. What grow-from-root is, and is not

Grow-from-root (He, Yalov and Hahn 2019, "XBART") builds a tree top-down in
one shot. Starting at the root with all n observations, it scans every
candidate cut of every variable, forms the integrated-likelihood weight of
each cut under the leaf prior, samples one cut (or the no-split option) from
that full weighted set, partitions, and recurses into the two children. One
sweep re-draws every tree in the forest this way against the current
residual. It is a stochastic tree BUILDER, not a walk: each sweep discards the
previous structure and grows a fresh one.

It is NOT a reversible Metropolis-Hastings kernel targeting the exact BART
posterior. The classic sampler perturbs a tree by one local move per sweep
(birth, death, change or swap) with a Hastings correction that makes the
forest chain pi-stationary for the true posterior. Grow-from-root has no such
correction: the recursive cut draw is a greedy-stochastic construction whose
stationary law is XBART's own approximation, close to the posterior in
predictive fit but not equal to it and not defended by a detailed-balance
argument. This distinction drives the surface decision in section 3.

## 2. The cut-scan kernel: measured cost

The kernel does, for one leaf's member set, one pass per variable over the
members, accumulating per-code sufficient statistics (count, sum w, sum wz,
optionally sum wz^2), then a prefix scan collapsing each per-code histogram
into the integrated-likelihood weight of EVERY candidate cut. Data layout is
faithful to bartcore's ColumnStore: uint16 codes, column-major (column j at
codes + j*N), members addressed through an index segment, so the scan gathers
col[idx[i]], z[idx[i]], w[idx[i]] - the access pattern of
misc_partitionIndices and misc_computeIndexed*SufficientStatisticsFast.

Two access modes are measured. "seq" is contiguous root-order members (the
root expansion of a grow-from-root build, and the best case: pure streaming).
"gather" is a scattered subset of the population (any non-root leaf, and the
regime informed proposals on an existing tree see). The baseline "single-move"
is the current per-move cost model: one gather partition on one column plus the
two-child suffstat recompute Tree::birth performs.

Apple M1 Max, single P-core, single-thread streaming read bandwidth ~52 GB/s
(64 MB out-of-cache sum, eight independent accumulators). Primary table:
255 cuts/variable, sum-wz^2 off. "scan" is one full node-expansion (all p
variables, all cuts); "mult" is scan/single-move; "GB/s" is demand traffic
(code 2 + z 8 + w 8 = 18 bytes per member-variable) over scan time.

  n_leaf   p  access   scan/node-exp  ns/elem  single-move  mult    GB/s
  1e5     10  seq         1.33 ms      1.33      0.51 ms     2.6x    13.6
  1e5     10  gather      5.29 ms      5.29      0.51 ms    10.4x     3.4
  1e6     10  seq        13.3  ms      1.33      4.86 ms     2.7x    13.6
  1e6     10  gather     48.6  ms      4.86      4.86 ms    10.0x     3.7
  1e5     50  seq         6.51 ms      1.30      0.545 ms   11.9x    13.8
  1e5     50  gather     30.3  ms      6.06      0.545 ms   55.6x     3.0
  1e6     50  seq        66.2  ms      1.32      4.63 ms    14.3x    13.6
  1e6     50  gather    247    ms      4.95      4.63 ms    53.4x     3.6

Stable across two runs (bandwidth 51-54 GB/s, scan times within a few
percent). Secondary findings from the full sweep (n_leaf in {1e3, 1e4} and
100 cuts, in the driver output): ns/elem is flat in n_leaf, confirming the
scan is O(p n_leaf) with no size-dependent per-element cost; 100 vs 255 cuts
is within noise (the histogram bins fit in L1 either way); accumulating
sum-wz^2 adds ~20% (an extra multiply and an extra in-cache bin write), so it
is real dead weight in the constant-leaf path (frontier 3.7 flags the same);
and gather at n_leaf <= 1e4 is optimistic because the member set fits in L2,
so the 1e5/1e6 rows are the honest DRAM figures.

DRAM-bound verdict, the two regimes differ:

- seq (root-order) runs at ~13.6 GB/s demand, a quarter of the 52 GB/s
  roofline, so it is NOT bandwidth-bound. It is bound by the histogram update
  - the indexed load-add-store into a bin per member, ~1.3 ns/element (~4
  cycles) regardless of n_leaf. This is the classic scatter-accumulate
  bottleneck, not memory.
- gather (scattered leaf) shows ~3-4 GB/s demand, but each of the three
  gathers touches a fresh 64-byte line, so true line traffic is ~3 * 64 / 5 ns
  ~= 38 GB/s, about 0.7x the roofline. The scattered scan IS effectively
  DRAM/cache-line bound, and column-major storage forces z and w to be
  re-gathered once per variable (they are the same field across variables), so
  a row-major or fused layout is the obvious ceiling-raiser if this ever
  matters.

Cost multiplier scan-vs-single-cut (the number frontier item 4 needs): on a
scattered leaf, scan/single-move ~= p - 10x at p=10, ~53x at p=50 - because one
classic move touches one column's data and the scan touches all p. On
root-order members the multiplier is ~p/4 (the sequential scan is ~4x cheaper
per variable than the gather-mode move it is compared against). Crucially, one
scan EVALUATES the entire p x cuts neighborhood (12750 candidates at p=50,
255 cuts) for that ~53x move cost, i.e. ~240 candidates per unit of classic
single-candidate move cost, so per candidate examined the scan is orders of
magnitude cheaper - the multiplier only looks large against a single blind
candidate.

Threading (gather, p=50, 255 cuts, columns split across workers): 1.80x at 2
threads and 2.93x at 4 threads for n_leaf 1e5; 1.88x / 2.4x for n_leaf 1e6.
Sub-linear and worse at the larger leaf, consistent with the gather scan being
partly bandwidth-bound - workers contend for the shared memory controller. The
scan parallelizes, but a single socket saturates well before its core count.

## 3. Validity story and the surface choice

Three ways to ship, with the exact-posterior gates arbitrating anything
sampled:

(a) Warm-start initializer only. Grow-from-root produces starting trees; the
existing exact MH sweeps then own stationarity, so the initializer needs no
correctness guarantee at all - a warm start may be any distribution. dbarts
already ships installTrees and the warm.start argument, so the consumer exists
and the seam is a producer of that same tree-state object. This is XBART's own
documented usage (grow-from-root to reach a good fit fast, BART sweeps to
refine and quantify uncertainty).

(b) A documented approximate sampler surfaced to users. This ships a
knowingly-biased default and must clear the exact-posterior gates as a
posterior-changing kernel; it would compete with, not complement, the exact
sampler, and dbarts's contract is an exact BART posterior.

(c) A Gibbs-valid variant. The scan can drive an EXACT kernel - but by a
different route than grow-from-root. Frontier 3.1 uses this identical cut-scan
to build a locally-balanced informed birth/death proposal with Zanella (2020)
weights, whose acceptance collapses to a ratio of two scan sums and is
pi-reversible. That keeps MH and the exact posterior while buying the scan's
informed-proposal mixing. So the exactly-valid use of the scan already has a
cheaper, MH-preserving claimant that does not require abandoning detailed
balance.

Recommendation: (a). Grow-from-root's unique, defensible niche is the warm
start - it sidesteps validity entirely, reuses the kernel, and matches the
method's own intended role. If we want the scan's benefit INSIDE the exact
sampler, the vehicle is frontier 3.1's informed proposal, not a standalone
grow-from-root posterior sampler. Ship neither (b) nor (c)-as-grow-from-root.

## 4. The engine seam, honestly

MoveStrategy is design vocabulary (core-generalization.md, gpu-bart.md,
gp-leaves.md); there is no MoveStrategy symbol in the code. What exists is a
per-tree loop in Chain::run (chain.hpp ~600-653) that, for each tree,
rolls the residual, calls setNodeAverages, then calls the free function
metropolisJumpForTree (moves.hpp) for exactly one local move, then
sampleParametersAndSetFits. The proposal draws live on CGMTreePrior
(model.hpp: drawSplitVariable, drawRuleForVariable) and the DART variant on
DartPrior; both are residual-independent (frontier fact 1.1).

A root-down builder would replace the metropolisJumpForTree call for that tree
with a full rebuild, and needs:

- A per-leaf all-cuts scan primitive over the ColumnStore + a node's index
  segment - the kernel prototyped here, promoted from benchmark to engine and
  reading the leaf model's suffstat rather than fixed (sumW, sumWZ, sumWZ2).
- A tree-builder entry parallel to metropolisJumpForTree: a free function,
  growTreeFromRoot(ctx, tree, rng, treeY, sigma), that recursively splits from
  the root, sampling each cut from the scan's integrated-likelihood weights
  and stopping on the no-split draw or the depth/min-node vetoes CGMTreePrior
  already enforces. It writes a fresh Tree in place of the incumbent.
- Interaction with the leaf-model concept: the constant leaf's per-cut
  accumulator is the three scalars measured here, and its integrated
  likelihood is the closed form the prefix scan uses. Linear and GP leaves
  replace the scalars with a small quadratic form (U'WU, U'Wz) per cut - the
  scan generalizes but each bin becomes a tiny matrix, connecting to the U'WU
  cache of linear-leaf-reuse.md. Constant leaf first; non-constant leaves are
  a later generalization, not a blocker.

The scan primitive is shared with frontier 3.1 and 3.3 and is the natural
first resident object for any future device backend (frontier section 4), so
it earns promotion independently of which consumer lands first.

## 5. Kernel verdict and break-even

The scan clears its cost bar for the warm-start role outright: a root-down
build is dominated by the sequential root pass (all n, ~13.6 GB/s, histogram-
bound), each deeper level scans a partition of the same n, and shallow BART
trees are 2-3 levels, so a full tree build costs ~ (levels) x (a root scan)
~= a few seq scans of n. Against a classic sweep that moves each tree by one
MH step, grow-from-root does far more structural work per sweep, but XBART's
payoff is reaching a good fit in far fewer sweeps - exactly the warm-start
economics, where absolute per-sweep cost is amortized by needing few of them.

For the informed proposal (frontier 3.1, gather regime), the break-even is
sharper: it costs ~p classic moves (10x at p=10, ~53x at p=50), so it pays for
itself only if it raises ESS-per-sweep on the tree-structure coordinate by
more than ~p. Informed proposals on structured discrete targets routinely
deliver one to two orders of magnitude ESS gains, so p=10 clears comfortably
and p=50 is at the edge - closed by the free Rao-Blackwellized readout
(frontier 3.3) and the lifted direction bit, which make the whole scanned
neighborhood count rather than the single realized move. The scan is
DRAM-bound in this regime, so wide hardware and the row-major/fused-field
layout are the levers if the multiplier needs shrinking.

## 6. Recommendation

GO on the cut-scan kernel as a shared engine primitive, and on grow-from-root
in role (a), a warm-start tree producer feeding the existing installTrees /
warm.start surface. NO-GO on shipping grow-from-root as a standalone posterior
sampler (roles b and c-as-grow-from-root); the exact-sampler use of the scan
belongs to frontier 3.1's informed proposal, which preserves MH.

Follow-up item for VD to approve: "grow-from-root warm-start producer" -
promote the scan primitive into the engine (constant leaf), add growTreeFromRoot
alongside metropolisJumpForTree, and wire it to emit warm-start tree state.
It shares its one new primitive with the frontier's informed-proposal
prototype (item 4), so the two can be scheduled to land the scan once and split
the consumer work. The implementation plan is a separate item.

## 7. Landing note (2026-07-10)

Shipped in role (a), constant leaf only. The cut-scan is a self-contained
header (src/bartcore/scan.hpp): a per-code (count, sum w, sum wz) histogram
over a node's index segment, prefix-scanned into every ordinal cut's collapsed
left/right suffstats and scored by the leaf's marginal. It is occupancy-aware
(a zero-count side emits a -inf never-selected sentinel, so no empty leaf is
ever built and the MH veto never runs on this path) and omits sum wz^2 by
evaluating the marginal with a zero sumWeightedResponseSq - the dead-weight
term's per-node total cancels between the no-split term and every cut. The
header is leaf-templated so frontier item 4 includes it unchanged and a
matrix-valued leaf substitutes its (U'WU, U'Wz) block without an interface
change.

growTreeFromRoot (src/bartcore/grow.hpp) recurses on the scan: at each node it
draws one outcome from {no-split} U {occupancy-nonempty (var, cut)}, weighted by
the CGM prior's own factors ((1 - growth) L(node); growth P(var) P(cut)
exp(L_left + L_right)) assembled in log space with a max-shift before
exponentiating. Draw count per node is exact: one discrete draw at every
positive-growth node, plus one symmetric missing-direction coin per split on a
column with missing values. Categoricals are ordinal-only in v1 (never scanned,
never split), which keeps that count exact and leaves categorical structure for
the MH sweeps; occupancy on the non-missing counts subsumes the ancestor split
interval and keeps every child non-empty, so the grown forest is a legal chain
state. Chain::growForestFromRoot duplicates run()'s sweep body with the per-tree
MH move replaced by a fresh grow, so the default run path stays byte-identical
(equivalence.R: 21/21 identical draws).

Surface: Sampler::growFromRoot fans across chains on the thread pool
(thread-count-independent, single chain inline on R's stream); the R5 method
dbartsSampler$growFromRoot(n.sweeps = 2L, updateState = FALSE) refuses linear
and gp node priors; bart2 gains n.grow.sweeps = 0L at the initialization fork
(0 keeps the prior-init default byte-identical; > 0 grows then samples as
usual), mutually exclusive with warm.start. The composable cross-sampler
workflow donor$growFromRoot(k) then target$installTrees(donor) reuses the
existing warm-start seam with no new install code. bart() (BayesTree-compat)
does not gain the argument. The exact-sampler use of the scan remains
frontier 3.1's informed proposal, unaffected.
