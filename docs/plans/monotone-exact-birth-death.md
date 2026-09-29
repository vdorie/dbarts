# monotone-exact-birth-death: the monotone birth/death move targets the documented prior

Status: PLANNED 2026-09-29 (dec-B144); derivation and gate verified numerically on an R prototype, not implemented

agent: opus (engine numerics: move seam, order counting, gate)
rng: posterior-changing for every fit with an active monotone constraint (all of its draws move); unconstrained
fits byte-identical, since every engine change lives in the monotone instantiation
window: before 1.0-0 (TODO monotone-exact-birth-death)
budget: ~800 lines (engine ~250, bridge and R ~60, tests/cpp ~200, tinytest ~60, gate ~180, docs ~50). Plan
estimates have run 1.5-2x low: expect up to ~1,600.

## Goal

A monotone birth or death is accepted with the collapsed ratio of the documented prior: the CGM tree prior, and
given the tree, iid (c-inflated) normal leaves restricted to the cone C(T) and normalized per tree by
Z_T = P(unconstrained leaves lie in C(T)). Every path that installs a monotone state leaves a feasible one, no
path reaches a structural move other than birth/death, and an exact gate over multi-split trees, which the
current move fails decisively, guards it.

## Context

- The defect. [`MonotoneConstantGaussianLeaf::oneLeafLogMarginal`](../../src/bartcore/model.hpp) and
  [`MonotoneConstantGaussianLeaf::twoLeafCoupledLogMarginal`](../../src/bartcore/model.hpp) divide the touched
  leaves' constrained marginal by d, their prior cone mass given the frozen neighbours. d equals the whole-tree
  ratio only when the touched leaves have no frozen constrained neighbour, which covers every move of part (a)
  of the existing gate ([9. Gates](../design/monotone.md#9-gates)). That is why the gate passed.
- Derivation. pi(T, M) is proportional to p(T) prod_k phi_k(mu_k) 1{M in C(T)} lik / Z_T, with
  phi_k = N(0, s_k^2) and s_k = c scale / k when leaf k has a constrained neighbour, as
  [`MonotoneConstantGaussianLeaf::priorSd`](../../src/bartcore/model.hpp) sets it. Integrate the touched leaves
  given the rest (same): pi(T, same) = p(T) / Z_T * prod_same phi lik * I_T(same), where I_T is today's score
  numerator without the division. So the birth ratio is
  [p(T*) / p(T0)] [Z_T0 / Z_T*] [I_T*(same) / I_T0(same)] [q(T* -> T0) / q(T0 -> T*)], death its reverse,
  and in log terms today's value + log d0 - log d* + log Z_T0 - log Z_T*. The frozen phi_k cancel because a
  birth or death never changes whether a frozen leaf has a constrained neighbour: a neighbour of the split leaf
  stays adjacent to at least one child. This held on all 101 births of the 3x2 enumeration. The proposal
  terms in [`birthOrDeathMove`](../../src/bartcore/moves.hpp) are unchanged. The touched-leaf redraws
  ([`MonotoneConstantGaussianLeaf::redrawAfterBirth`](../../src/bartcore/model.hpp),
  [`MonotoneConstantGaussianLeaf::redrawAfterDeath`](../../src/bartcore/model.hpp)) are already exact
  conditionals, and the leaf Gibbs sweep is exact given T. An empty cone (a death whose merged leaf's frozen
  bounds cross) keeps its -HUGE_VAL sentinel, which is correct: pi(T*, same) = 0 there.
- Which pairs are constrained. [`monotoneNeighborBounds`](../../src/bartcore/model.hpp) relates j < k when,
  along a constrained axis, j's code box ends one code below k's start (direction -1 flips it), and the two
  boxes share at least one code on every other axis on either leaf's path. Boxes that touch only at an edge or
  corner share no code and are not a pair; their order follows by transitivity through the leaves between.
  The relation is acyclic by recursion on the tree. At a split on axis v, the only cross pairs run along v, in
  one direction. A split on an unconstrained axis leaves no cross pairs, so its subtrees fall into separate
  components. Related leaves all carry sd c scale / k and isolated leaves drop out, so
  Z_T = e(P_T) / L! = prod over components C of e(C) / |C|!, with e the number of linear extensions. This was
  checked against a brute-force oracle (iid draws, fitted function tested on the cell grid) on 32 trees with
  1-3 constrained axes and mixed directions: every |z| <= 1.6.
- No closed form. A tree that splits only on one constrained axis is a chain (e = 1). One constrained and one
  free axis already produce non-series-parallel components: the 3x3 enumeration has an N (e = 5 of 4!).
  Counting is exact by a DP over down-sets per component, costing O(down-sets x |C|). Only the component
  holding the touched leaf, or its children, differs between T0 and T*.
- Measured sizes (old sampler, n 1000, p 5, 75 trees, 1/2/3 constrained predictors, 7,500 trees each):
  at most 9 leaves, largest component 9, at most 43 down-sets. The C++ count averages 0.14 us.
- Stress (random guillotine partitions, worst of 20): two constrained axes, 116 us at L = 40. One constrained
  plus one free axis, 1.8 ms at L = 24, 52 ms at L = 32 and 0.1 s at L = 40 (850k down-sets).
- A monotone tree step costs ~120 us today (2.4 ms per 20-tree sweep, quadrature-bound), so the count adds
  under 5% at realistic sizes.
- Evidence, one-tree enumerations on cell grids with 10 rows per cell and batch-means chi-square against the
  exact law. The R prototype of the corrected move passes: 17/14 df (1-D, 4 cells), 40/60 (x1, x2
  constrained, 3x2), and 16-48 on 49 df at weak signal. The same prototype with the old ratio reproduces the
  dbarts draws (11/14, 34/51, 28/60), so it mirrors the code. The current dbarts fails with 767/14, 641/51 and
  429/60. On one slow-mixing design (root mass 0.6%) the corrected prototype sat 2.8 regenerative SEs off in
  the root variable. The same geometry at root mass 3-5% passes (|z| <= 0.6), so gate designs must mix.
- Reachability. At creation R forces birth/death ([`resolveSamplerSpec`](../../R/spec.R)).
  [`xbart()`](../../R/xbart.R) refuses `monotone`. [`ruleGibbsMove`](../../src/bartcore/moves.hpp) compiles out
  for a param-scoring leaf. But:
  - [`dbartsSampler$setControl`](../../man/dbartsSampler-class.Rd) installs any mixture on a monotone sampler.
    Probed: the default 0.6/0.4 change mix is accepted, so change and swap run with a score that is exact only
    for birth/death.
  - [`dbartsSampler$installTrees`](../../man/dbartsSampler-class.Rd) and
    [`dbartsSampler$setState`](../../man/dbartsSampler-class.Rd) accept an unconstrained donor's leaves.
    Probed: the fitted surface falls by up to 1.4 along the constrained axis, and by 0.48 after one sweep.
  - [`Chain::growForestFromRoot`](../../src/bartcore/chain.hpp) draws leaves against the previous tree's mu at
    reused node ids, so neighbour bounds can cross and hit the clamp fallback of
    [`MonotoneConstantGaussianLeaf::drawTruncatedNormal`](../../src/bartcore/model.hpp) with a > b. No
    violation was seen in one probe.
- Unaffected. The prior draw ([`MonotoneConstantGaussianLeaf::drawFromPriorForTree`](../../src/bartcore/model.hpp))
  already samples the documented prior by rejection. The level-fibre shift is a conditional given T and
  prediction reads leaves. The help pages describe the prior being kept, and monotone is new in 1.0-0, so NEWS
  gets no entry.
- Claims to reword: [4. Decision - marginal likelihood for the structure moves](../design/monotone.md#4-decision---marginal-likelihood-for-the-structure-moves)
  (B "targets the EXACT constrained posterior", B' "exact posterior up to 1-D quadrature error"),
  [11. Costs, risks, and confidence](../design/monotone.md#11-costs-risks-and-confidence) (confidence HIGH),
  [Plan-vs-code note](../design/monotone.md#plan-vs-code-note), and dec-B16 in [decisions.md](../decisions.md).

## Decision

Pathological sizes. Is the count uncapped, or does it run under a work budget?

Recommendation: a budget of 2^22 down-sets per component, under ~0.5 s and ~200 MB at worst. A birth whose tree
exceeds it is rejected, and the prior draw rejects such trees, so the target is exactly the documented prior
restricted to trees within the budget, stated as part of the prior. This removes no approximation (dec-B14).
Only a component of 23 or more leaves can exceed the budget, and the CGM(0.95, 2) prior puts 1.4e-18 on 23 or
more leaves; no measured fit came within 5 orders of it.

The alternative, no budget, keeps the prior verbatim, but one pathological tree can stall a sweep for minutes
and exhaust memory.

Evidence that would change this: a real fit whose count exceeds ~1 ms.

## Constraints

- Exact for the stated prior: no approximation (dec-B14). Counts accumulate in double: exact to 2^53, with
  relative error under 1e-12 beyond that, like the quadrature tolerance already in the score.
- Unconstrained samplers byte-identical: the new seam compiles out, as the three existing monotone seams do.
  No dbarts.h change.
- Out of scope:
  - an exact linear-extension draw for the prior leaves, which would replace rejection and its 1e6 cap
    (a 9-leaf chain fails the cap ~6% of the time);
  - change moves under the constraint;
  - quadrature speed (TODO monotone-leaf-quadrature);
  - empty leaves pinned at 0, which are transient under the veto.

## Steps

1. Order counter (engine, beside the monotone geometry in model.hpp). The input is a tree and a leaf. The
   counter collects that leaf's component by BFS over the same adjacency test as `monotoneNeighborBounds`,
   then runs the layered down-set DP on a 64-bit mask (a component over 64 leaves is over budget). It returns
   log e(C) - lgamma(|C| + 1), or an over-budget flag.
2. Seam. Add an optional leaf concept, `logTreeNormalizer(tree, leaves)`, which the monotone leaf declares.
   `birthOrDeathMove` evaluates it on the touched leaves in each state: before and after `tree.birth`, and
   before and after `orphanChildren`. It multiplies the prior ratio by Z_T0 / Z_T*, and the census logs that
   ratio in the prior column. An over-budget proposal rejects.
3. Drop the d terms: `priorMass` in `oneLeafLogMarginal`, `denom` in `twoLeafCoupledLogMarginal`. Leave
   `coneProbability` in place.
4. Budget in the prior draw. [`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp) adds "within budget"
   to its whole-tree rejection predicate, for monotone only.
5. Reachability:
   - A bridge refusal, called after [`parseProposalProbs`](../../src/R_interface_bartcore.cpp) at creation and
     in [`bartcore_setModel`](../../src/R_interface_bartcore.cpp): an active monotone constraint with any
     nonzero swap, change, perturb or rule_gibbs probability errors (the frozen all-zero mixture stays
     allowed). setControl names the refusal in R.
   - `growForestFromRoot` reseeds a monotone tree's mu to the all-zero feasible seed before its draw.
   - [`Chain::installForest`](../../src/bartcore/chain.hpp) reseeds any installed tree failing
     [`monotoneTreeIsFeasible`](../../src/bartcore/model.hpp): a warm start keeps the donor's structure, not
     values the constraint forbids.
   - [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) refuses such a tree on the setState path, so a
     state restores exactly or not at all.
6. tests/cpp:
   - The count against brute-force permutation counts on hand-built trees (1-3 axes, mixed directions, one N)
     and 200 random trees.
   - Log Z against e / L! on the enumeration's trees.
   - An RNG-free ratio test: on hand-built A < B -> A < B1 < B2 and A < B -> {A < B1, B2} trees with pinned
     mu_A, the move's log ratio against an independent closed form. It must differ from the current code's
     value.
   - [`testMonotoneMarginal`](../../tests/cpp/test_model.cpp) loses its d_* = 1/2 normalizer check.
   - The budget flag fires on a star component above 2^22.
7. tinytest ([test-monotone.R](../../inst/tinytest/test-monotone.R)): setControl with a change mix errors;
   installTrees from an unconstrained donor leaves the fit monotone at once; setState of that state is refused;
   growFromRoot then one sweep is monotone. A statistical check does not fit here: the most sensitive cheap
   functional (the root cut on the 1-D design) sat at |z| 0.6-0.8 against the old move at 50k-100k draws and
   15-30 s.
8. Gate: benchmarks/R/monotone-reference.R gains part (c), reusing
   [`runSampler`](../../benchmarks/R/monotone-reference.R): an exact enumeration over multi-split one-tree
   structures. Port the check's enumeration: CGM prior, cell sums, per-tree
   exp(sum base) P_post(C) / Z_T, with P_post(C) summed over linear extensions. Map sampled trees to canonical
   keys and test with a batch-means chi-square: fail above qchisq(1 - 1e-4, df), and report max |z|. Designs
   (10 rows per cell):
   - c1: x1 constrained, 4 cells, sigma 0.6 (15 structures);
   - c2: x1 and x2 constrained, 3x2, sigma 1.5 (62);
   - c3: x1 constrained and x2 free, 3x2, sigma 1.0 (62).
   Quick mode runs 300k draws each, adding about 6 min (block-wise getTrees reads) to the 40-min
   exact-gates job; full mode runs 900k and two seeds.
   Power at quick size: the current move scores 244 (bound 43), 213 (95) and 203 (95). The corrected
   prototype scores 11, 48 and 47.
9. A monotone scenario in benchmarks/R/equivalence.R (x1, x2 constrained, 20 trees). The gaussian baseline
   re-records with that scenario added and the other 53 bitwise; its MANIFEST row names part (c) as the ORACLE
   (P17). No monotone scenario exists today, so nothing else re-records.
10. Docs:
    - monotone.md sections 4, 9 and 11 and the Plan-vs-code note restate B' with the whole-tree normalizer, and
      section 9 gains part (c). Also record that mBART's d-normalized eq. 4.11 targets neither this prior nor
      the software's.
    - dec-B16's Record line: superseded in part by dec-B144.
    - Status lines and INDEX at landing, and the TODO item removed.

## Verification

- `cd tests/cpp && make && ./test_bartcore`: the new count, ratio and budget checks pass.
- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'`:
  all pass.
- `R_LIBS=<lib> Rscript benchmarks/R/monotone-reference.R quick`: parts (a), (b) and (c) PASS. Mutation run:
  restore the two d divisions and drop the Z term, then `touch` the header and reinstall; each of c1-c3 must
  FAIL at 2x its bound or more.
- The full exact-gates.yaml list with `quick` on the slice's library: every gate PASSES.
- `Rscript benchmarks/R/equivalence.R compare <current>`: 53 scenarios "identical draws (same RNG stream)", no
  "max |z|". The BCF and multinomial harnesses are identical. The four test-reproducibility snapshot files
  carry no monotone fit and pass unchanged on the reference build.
- Release level: the SBC monotone arm, whose design section in [sbc-family-tiers.md](sbc-family-tiers.md) is
  still in review. Its `monotone-1` arm (1 tree, 0.4 s per replicate, a few minutes at R 400) flags the current
  move and must pass. Then the 20-tree arm runs (~85 min at R 200, 4 h at 50 trees) before admission to the
  matrix.
- Speed: on a quiet machine, monotone sweep time before and after at 20 trees with 1 and 2 constrained
  predictors must stay within 5%, and bench-sampler compare must show the unconstrained paths unchanged.
- `Rscript tools/check-doc-freshness.R .` passes.
