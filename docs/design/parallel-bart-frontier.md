# Parallel BART: the frontier

Status: research note, 2026-07-07. Companion to gpu-bart.md, which surveys
what exists; this note works out what does not exist yet, from BART's
mathematical shape, and ranks it. Everything here is grounded in this
engine's actual kernels (src/bartcore/{chain,model,moves,tree}.hpp) and is
stated with an exactness class, a cost model, and the cheapest falsifying
experiment. Observation-axis parallelism inside one tree update is treated
as settled (bartz; gpu-bart.md section 2g) - the question is every OTHER
axis, plus what a CRAN package can conditionally ship.

## 1. Three structural facts everything else builds on

1.1 Structure proposals are residual-independent. The node, the split
variable, and the cut are all drawn without consulting the working
response (moves.hpp: drawBirthableNode, drawSplitVariable,
drawRuleForVariable, and the change/swap selections); only the acceptance
ratio and the leaf/sigma draws touch the residual. Consequences: a whole
sweep's proposals can be drawn and scored against a frozen residual and
corrected exactly afterward (section 3.2); two coupled chains propose
identical moves as soon as their structures agree (3.6); and an informed
proposal needs the data only for its weights (3.1).

1.2 The constant leaf's inter-tree coupling is one scalar field. The
constant-leaf sampler consumes only (sumW, sumWZ) per node: the sumWZ2
term cancels algebraically in every move ratio (it is additive across any
partition of a fixed observation set), and sumWZ is affine in the
residual - sumWZ_t(L) = D(L) + mu_t(L) sumW(L), with D the leaf sum of the
single field d_i = w_i (y_i - F_i). So the entire coupling of one tree to
the other m - 1 factors through d. There is no higher-dimensional
interaction to exploit: cross-tabulations, low-rank factors, and cell
algebra over the FULL forest are provably redundant with d and more
expensive to maintain (section 2). This wall is specific to the constant
leaf; linear and GP leaves have quadratic-form marginals, which is exactly
why the U'WU cache (linear-leaf-reuse.md) exists for them.

1.3 Backfitting is already the optimal exact coloring of a small joint
Gaussian. Conditional on structures, the forest is a linear model over
L = sum_k L_k leaf indicators; L = O(m) (shallow trees, 2-5 leaves), so
the joint leaf system is TINY. Same-tree leaves never share observations,
so each tree is an independent set, and the classic leaf stage is
precisely one Gauss-Seidel sweep of this Gaussian with one color per
tree. At default depth the cross-tree leaf graph is nearly complete
multipartite, so no coloring with fewer than m colors exists: exact
cross-tree width at the leaf layer is unavailable, which is the scarcity
blocked-jacobi manufactures around. Meanwhile the leaf prior precision
lambda = k^2 m / nodeScale^2 (= 16m at defaults) GROWS with m, so the
joint system is best-conditioned exactly where backfitting mixes worst -
a gift specific to BART's 1/sqrt(m) prior scaling.

## 2. Conservation laws: what cannot be done

These negative results are load-bearing; each was reached independently
and kills a tempting construction.

- The precision budget is conserved. Any augmentation that (a) renders
  the tree targets conditionally independent given latents, (b) is
  Gaussian, and (c) preserves the marginal N(sum g_k, sigma^2) must
  allocate total variance sigma^2 across the pieces. Per-tree freedom
  summing to sigma^2 cannot equal sigma^2 everywhere; negative
  correlation restores per-piece variance but destroys the conditional
  independence that was the point. Blocked-jacobi's conservative-birth
  tax is therefore a law, not an implementation artifact. What remains is
  ALLOCATION (section 3.5).
- Whole-sweep fused algebra is foreclosed twice over. The m-way joint
  partition's cell count saturates to n by m ~ 16-64 (simulated census;
  the real-forest census is experiment E1), so "dense cell algebra for
  the whole sweep" is O(n^2)-shaped; and even with few cells, fact 1.2
  makes full-forest cross-tabs redundant with the residual field and
  O(m^2 n) to maintain against O(mn) for the field itself. The one-pass
  "score every tree against a frozen residual, then run the sweep as
  algebra" reorganization is exact for exactly one tree; consuming it
  for all m without refresh IS blocked-jacobi (posterior-changing), not
  a free restructuring.
- The determinant lemma solves the wrong half. Full-forest-marginal
  structure moves are rank-2 in the evidence (a birth splits one
  indicator column), so their linear algebra is O(L^2) and free. The
  irreducible cost is the co-occurrence scan - which leaf of every other
  tree each affected observation sits in, O(n_leaf m) per move - which is
  precisely what classic BART's residual compression avoids. Exact
  escapes exist only via delayed acceptance (3.2) or small blocks.
- No exact zero-coupling coloring of structure moves exists: the
  posterior leaf covariance A^-1 is dense even when the Gram is sparse,
  so two moves' capacitance cross-block never exactly vanishes. Weak-
  coupling waves are possible only as a bounded-bias approximation,
  gated like any posterior-changing kernel.
- Plain Jacobi (all trees against the frozen residual, no augmentation)
  targets the wrong stationary distribution - the baseline trap.
- Joint accept/reject of one move per tree dies as rho^m.
- Cross-tree pipelining and speculation on the constant leaf are
  foreclosed by the unconditional per-tree leaf redraw: the residual
  entering tree t+1 always changes, and resolving tree t is already
  cheap, so there is nothing expensive to hide. (Speculative overlap
  becomes real for GP/linear leaves, whose per-tree serial step is an
  O(basis^3) factorization - noted in 3.7.)

## 3. Constructions that survive

Ordered by (exactness class, engine reach, shippability). "Shifting"
means same distribution, different floating-point draws - the e565326
suffstat gate class; "new-kernel" means the posterior is preserved but
the kernel differs, gated against exact posteriors.

3.1 Locally-balanced informed birth/death over the full-cut scan
(new-kernel, exact). One segmented prefix scan per variable over a
leaf's members produces the collapsed marginal for EVERY candidate cut;
propose among the whole birth/death neighborhood with Barker weights
g(u) = sqrt(u), and the acceptance collapses to min(1, Z(T)/Z(T')), a
ratio of two scan sums (Zanella 2020 gives pi-reversibility for any
symmetric neighborhood; birth/death is symmetric including the validity
vetoes). Change/swap stay classic (their neighborhoods do not scan).
Cost rises from one cut to all cuts - O(p n_leaf) per move, flat on wide
hardware - and buys informed proposals whose ESS-per-sweep gain on
structured discrete targets is routinely large. The scan is the same
object grow-from-root needs, the same object a GPU wants resident, and
the same object 3.3 recycles. Falsifier: single-tree enumerable problem;
confirm exactness, then ESS-per-sweep vs the cost multiplier.

3.2 Delayed acceptance, twice - the convergent discovery. Both the
kernel and function-space analyses independently arrived at the same
two-stage exact pattern (Christen-Fox):
  (a) Batched cross-tree scoring: draw all m proposals (residual-
      independent, fact 1.1), score them in one fused pass against the
      frozen start-of-sweep residuals, then run the sequential pass
      where stage 1 uses the precomputed stale score and stage 2 pays a
      true-residual rescan ONLY for stage-1 survivors (~20-40%). The
      serial critical path shrinks 2.5-6x; the bulk scoring is one wide
      launch.
  (b) Collapsed-forest structure moves: stage 1 is the CLASSIC
      conditional ratio (cheap, already computed); stage 2 pays the
      O(n_leaf m) co-occurrence plus a rank-2 determinant-lemma update
      to score the move against the FULL forest marginal - Rao-
      Blackwellizing away the leaf-noise in move decisions - only for
      stage-1 survivors.
  Both are exact compositions. Falsifier for (a) is instrumentation
  only: log stale-vs-true acceptance agreement and survivor rate on a
  stock run. For (b): matched-RNG ESS comparison classic vs collapsed at
  small n.

3.3 Waste recycling / Rao-Blackwellized readout (exact, free once a
scan exists). The informed scan scores the entire neighborhood; posterior
functionals (variable inclusion, DART split counts, pointwise fits) can
be averaged over the neighborhood's selection distribution instead of
the single realized move - unbiased, variance never larger. Plus a
lifted (non-reversible) grow/prune direction bit for the tree-size
coordinate: ballistic instead of diffusive motion on the slowest
coordinate, exact, nearly free to implement.

3.4 Block-fused sub-sweeps over persistent atom maps (shifting-class;
the flagship CPU candidate). Process trees in consecutive blocks of
b ~ 4-8. Within a block, maintain the joint partition of the b trees
("atoms": at b = 8 the simulated census gives hundreds to low thousands
of atoms, vs n as b grows past ~16) with per-atom statics A(c), P(c),
frozen-outside fit Q(c), and the running in-block fit; then every leaf
draw, residual roll, rejected-move bookkeeping, and fit write inside the
block is O(atoms) algebra instead of O(n) streams - the identical
Gauss-Seidel sweep, arithmetic regrouped. Births/changes/swaps still
touch their leaf's data (a new cut slices atom interiors); that is ~7%
of today's traffic, while the field maintenance this removes is ~93%.
Modeled DRAM traffic drops ~7-10x at b in {4, 8} for the confirmed-
common n >= 1e5 single-chain workload; atom maps persist across sweeps
and update in O(changed leaves), generalizing the landed U'WU cache
(residual-independent structure cached, residual-dependent part carried
incrementally). Exactness: same distribution, different summation
grouping - b = 1 must reproduce current draws bitwise, b > 1 gates
statistically plus exact-posterior checks. Composes with within-chain
threading (fewer barriers AND cheaper passes) and with a GPU backend
(fewer resident-to-DRAM round trips). Falsifiers BEFORE any engine work:
E1 atom census on real fitted forests (correlated trees, not simulated
ones) at b in {4, 8, 16}; E2 field-fraction profile confirming the
~90/10 maintenance/scan split.

3.5 Unequal noise splitting as a movement-budget schedule (new-kernel,
exact). Within the conservation law (section 2), the allocation {v_k},
sum v_k = sigma^2, is free: release a rotating handful of trees at
near-full variance and pin the rest, updating the batch independently in
parallel. Exact for any state-independent rotation schedule; the uniform
tax becomes a schedule that spends exploration where movement is wanted.
Falsifier: b = 2 released/pinned vs uniform on a signal-concentrated
synthetic, ESS/sec at equal width. This is the surviving core of the
blocked-jacobi question: not "can the tax be avoided" (no) but "can it
be aimed" (yes, exactly).

3.6 Unbiased BART via couplings (exact/unbiased; width for length). Two
lagged chains sharing proposal randomness meet exactly - discreteness
makes meetings positive-probability events, and fact 1.1 means identical
structures propose identical moves; maximal couplings on accept
indicators, leaf draws, and sigma cascade a structure meeting into a
full meeting. Each met pair yields an unbiased estimator with no burn-in
bias; pairs are embarrassingly parallel. The failure mode is not bias
but heavy-tailed meeting times (infinite variance). Falsifier: meeting-
time census at m = 1, then m = 5, watching the tail index. Uniquely
suited to fleets of small devices rather than one big one.

3.7 Substrate and side findings.
- Joint leaf draw by warm-started preconditioned CG (perturbation
  optimization; backfitting itself is the block-Jacobi preconditioner;
  warm starts work because O(1) leaves change per sweep). Exact to
  tolerance at ~2.5-7x a sweep; its value is mixing (removing leaf-layer
  sloshing) and being 3.2b's substrate, never FLOPs. Decisive number:
  does the preconditioned iteration count stay single-digit?
- sumWZ2 is dead weight in the constant-leaf hot path (fact 1.2): it
  cancels in every decision. Removing it saves a stream but changes
  rounding - a shifting-class micro-item for the next baseline
  re-record, not a free win.
- Speculative overlap of the GP/linear leaf's O(basis^3) resolve with
  the next tree's residual-independent work is bitwise-exact latency
  hiding where the serial step is genuinely long.

## 4. The CRAN axis: shipping a conditionally-enabled GPU path

Three mechanisms, in increasing order of ambition, none requiring CRAN
to build device code:
- Configure-conditional OpenCL (precedent: the OpenCL package) - builds
  everywhere, enables where a runtime exists; aging ecosystem.
- A backend-seam ABI: dbarts ships CPU-only and exposes a registered
  kernel-provider interface in dbarts.h - the same additive registration
  pattern as dbarts_sampler_setCallback - and a companion package
  (CRAN + OpenCL, or R-universe/GitHub + CUDA/Metal) registers device
  implementations at runtime. Conditional enabling becomes a runtime
  fact, not a build fact.
- The torch precedent: a CRAN package that downloads its compiled device
  runtime at first use. CRAN-accepted at scale; the strongest existing
  answer for CUDA specifically.

The design consequence that binds the math to the seam: per-call
offload of stateless kernels loses to transfer costs (gpu-bart.md's
falsifying experiment), so a viable seam gives the DEVICE resident
ownership of x (and its sorted orders), the residual field, and the
scan/atom workspaces, with the host driving tree logic through a command
stream. Ranked by fit to that shape: the full-cut scan trio (3.1 + 3.2a
+ 3.3) is the natural first resident object - one scan services the
informed proposal, the batched stale scoring, and the recycled readout;
block fusion (3.4) reduces command traffic for whatever backend exists;
the CG leaf solve (3.7) is matvec-shaped. The seam should be DEFINED
only after the CPU-side experiments pick winners - freezing a kernel ABI
before knowing which kernels matter is how the wrong interface ships.

## 5. Ranked next actions

Instrumentation before prototypes; each item names its kill condition.

1. E2 field-fraction profile of the current sweep (hours; kills or
   confirms 3.4's ceiling). DONE 2026-07-08. SURVIVES (does not kill
   3.4). Measured per-sweep DRAM-byte split, constant-leaf sweep,
   m = 75, Friedman-1 (p = 10), n in {1e4, 1e5}, default (base 0.95,
   power 2) and deep (power 1) configs; driver
   benchmarks/R/parallel-falsifiers.R. Field maintenance (residual roll
   + whole-tree suffstat recompute + fit scatter + totalFits rebuild +
   change/swap index-segment snapshot/restore) = 82-88% of traffic;
   per-move scan work (partition + affected-subtree suffstat) = 12-18%
   (scan reported as a [lower, upper] band on the index-swap writes).
   Dominant and clear of the ~70% collapse floor, but below the note's
   ~93% model: the biggest correction is that shallow-tree change/swap
   moves re-partition near the root, so scan is ~2x the modeled 7%. The
   removable field is ~85%, so the atom map's DRAM drop is ~5.5-6.6x
   (scan-upper band), i.e. the 7-10x claim lands at its low end / just
   under. Stable across n and config (field share drifts +2pt from
   n = 1e4 to 1e5 as trees settle). Field is dominant: E2 clears.
2. E1 atom census on real fitted forests, b in {4, 8, 16} (hours; kills
   3.4 if real atoms saturate at small b). DONE 2026-07-08. SURVIVES
   decisively. Distinct per-observation leaf-assignment tuples over
   consecutive tree blocks, 10 post-burn-in kept forests (nskip 200,
   keepevery 20), two block positions, same configs/data as E2, seed 99.
   Mean distinct atoms (atoms/n; occupancy n/atoms):
     b = 4 : default n1e4 43 (0.004; 233), n1e5 74 (0.0007; 1355);
             deep n1e4 51, n1e5 113 (0.001; 886).
     b = 8 : default n1e4 312 (0.031; 32), n1e5 1120 (0.011; 89);
             deep n1e4 318, n1e5 1596 (0.016; 63).
     b = 16: default n1e4 2196 (0.220), n1e5 12618 (0.126);
             deep n1e4 2293, n1e5 20755 (0.208).
   At b in {4, 8} atoms sit 2-4 orders of magnitude below n and atoms/n
   FALLS as n grows (sub-linear atom growth: 10x n gives ~1.7x atoms at
   b = 4, ~3.6x at b = 8), so the confirmed-common n >= 1e5 regime is the
   best case. Saturation toward n only begins past b = 16 (atoms/n ~0.13-
   0.23, occupancy still 4-8), exactly the ceiling the simulated census
   predicted. Real correlated forests do NOT saturate at small b: E1
   clears. The deep config raises atoms only marginally (same order).
3. Stale-residual agreement logging on a stock run (hours; kills 3.2a if
   the surrogate disagrees often or survivors are near-universal). DONE
   2026-07-08. SURVIVES. Per proposed move, the accept decision under the
   frozen start-of-sweep residual (batched-scoring surrogate,
   y - totalFits_start + oldFit_t) vs the true rolled residual, SAME
   proposal randomness (rng state saved/restored around a snapshot
   excursion; the true decision drives the stock chain). Same configs;
   nskip 100, ndpost 100. Frozen-vs-true agreement 0.966-0.978 overall
   (birth lowest at 0.92-0.98, change highest at 0.98-0.996); stage-1
   survivor (frozen-accept) rate 0.068-0.125 overall, ~= the true accept
   rate. High agreement (surrogate rarely misdirects) AND survivors far
   below universal - in fact BELOW the note's modeled 20-40%, so the
   serial-path shrink is if anything larger than 2.5-6x. Neither kill
   trips: 3.2a clears. (Caveat: agreement is a sweep average; it is worst
   for late-in-sweep trees carrying the most accumulated drift - a
   per-position refinement was not run, but the 0.97 average bounds it.)
4. Single-tree informed-kernel prototype vs an enumerated posterior,
   then ESS-per-sweep vs cost (days; decides 3.1 and by extension the
   scan trio).
5. Preconditioned-CG iteration-count instrumentation (day; decides 3.7's
   substrate and 3.2b's affordability).
6. Coupling meeting-time census at m in {1, 5} (day; decides 3.6).
7. b = 2 released/pinned noise-splitting ESS/sec (days; decides 3.5 and
   with it the remaining value of blocked-jacobi-trees).
If 3.4 survives E1/E2, it is the flagship engine item (CPU-first, exact
sweep, DRAM-bound workload) and proceeds independently of any GPU
decision. The backend seam is designed after items 1-5 report.

Verdict (2026-07-08, items 1-3 measured; drivers under benchmarks/R and
raw records in docs/plans/parallel-falsifiers.md). 3.4 needed BOTH E1 and
E2: both cleared, so 3.4 SURVIVES and is the flagship engine candidate,
with one correction to the cost model - measured field share is ~85% not
~93%, so the realized DRAM drop is ~6x (7-10x low end) rather than the top
of that range; the gap is shallow-tree change/swap moves re-partitioning
near the root, which enlarges the surviving scan term. 3.2a needed the
stale-residual logging to clear (high agreement AND sub-universal
survivors): it did, comfortably, so 3.2a SURVIVES and the batched-scoring
delayed-acceptance path earns a prototype. Both flagship CPU candidates
live; neither was killed. Next: items 4-5 (informed-kernel and CG-leaf
prototypes) before the backend seam.

## References

Zanella (2020), "Informed proposals for local MCMC in discrete spaces,"
JASA. Christen and Fox (2005), "MCMC using an approximation," JCGS.
Jacob, O'Leary, Atchade (2020), "Unbiased MCMC with couplings," JRSS-B.
Hastie and Tibshirani (2000), "Bayesian backfitting," Statistical
Science. Papandreou and Yuille (2010); Parker and Fox (2012) -
perturbation-optimization / CG Gaussian sampling. Frenkel (2004), waste
recycling. Entezari, Craiu, Rosenthal (2018), LISA. Petrillo (2024),
bartz, arXiv:2410.23244. See gpu-bart.md for the surveyed landscape
these build on.
