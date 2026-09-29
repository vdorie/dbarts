# monotone-exact-birth-death: the monotone birth/death move targets the documented prior

Status: PLANNED 2026-09-29 (dec-B144), revised after a blind critique; the size budget below is a decision for
the maintainer. Derivation and gate verified on an R prototype; not implemented.

agent: opus (engine numerics: move seam, order counting, exact pair redraw, gate)
rng: posterior-changing for every fit with an active monotone constraint (all of its draws move, prior draws
included); unconstrained fits byte-identical, since every engine change lives in the monotone instantiation
window: before 1.0-0 (TODO monotone-exact-birth-death)
budget: ~1,000 lines (engine ~420, bridge and R ~90, tests/cpp ~320, tinytest ~90, gate wiring ~10, docs ~70;
the gate script itself is already in the tree). Plan estimates have run 1.5-2x low: expect up to ~2,000.

## Goal

The monotone chain targets the documented prior exactly. That prior is the CGM tree prior and, given the tree,
iid (c-inflated) normal leaves restricted to the cone C(T), normalized per tree by
Z_T = P(unconstrained leaves lie in C(T)). This covers four pieces:

- the birth/death acceptance;
- the redraw of a birth's two children;
- the treatment of leaves left empty;
- the prior leaf draw.

No path reaches a structural move other than birth/death. Every state the sampler holds lies in the cone. An
exact gate over multi-split trees, N-shaped leaf orders included, fails the current engine decisively and
guards the fix.

## Context

- The defect. [`MonotoneConstantGaussianLeaf::oneLeafLogMarginal`](../../src/bartcore/model.hpp) and
  [`MonotoneConstantGaussianLeaf::twoLeafCoupledLogMarginal`](../../src/bartcore/model.hpp) divide the touched
  leaves' constrained marginal by d, their prior cone mass given the frozen neighbours. d agrees with the
  whole-tree normalizer only when no touched leaf has a frozen constrained neighbour. That covers every move of
  part (a) of the existing gate ([9. Gates](../design/monotone.md#9-gates)), which is why it passed.
- Derivation. pi(T, M) is proportional to p(T) prod_k phi_k(mu_k) 1{M in C(T)} lik / Z_T, with
  phi_k = N(0, s_k^2) and s_k = c scale / k for a leaf with a constrained neighbour
  ([`MonotoneConstantGaussianLeaf::priorSd`](../../src/bartcore/model.hpp)). Integrating the touched leaves
  given the rest ("same") gives pi(T, same) = p(T) / Z_T * prod_same phi lik * I_T(same), where I_T is today's
  score numerator with the division removed.
  - Birth, exactly:
    log alpha = log[p(T*) / p(T0)] + log[q(T* -> T0) / q(T0 -> T*)] + log I_T*(same) - log I_T0(same)
    + log Z_T0 - log Z_T*. Death is the reverse.
  - Against today's value this adds + log d* - log d0 + log Z_T0 - log Z_T*.
  - The frozen phi_k cancel because a birth or death never changes whether a frozen leaf has a constrained
    neighbour: 0 of 101 births in a 3x2 enumeration.
  - The proposal terms in [`birthOrDeathMove`](../../src/bartcore/moves.hpp) are unchanged.
  - An empty merged cone keeps its -HUGE_VAL sentinel, which is correct: there pi(T*, same) = 0.
- Constrained pairs. [`monotoneNeighborBounds`](../../src/bartcore/model.hpp) relates j < k when, along a
  constrained axis, j's code box ends one code below k's start (direction -1 flips it), and the boxes share a
  code on every other path axis.
  - Edge or corner contact is not a pair; its order follows by transitivity through the leaves between.
  - The order is acyclic: at a split on axis v the only cross pairs run along v, one way. A split on a free
    axis creates no cross pairs.
  - Related leaves share the sd c scale / k, so Z_T = e(P_T) / L!, with e the number of linear extensions.
  - The brute-force oracle agrees (zcheck mode below): 32 trees, 1-3 constrained axes, mixed directions, worst
    |z| 1.64.
- No closed form. A tree splitting on one constrained axis only is a chain (e = 1), but one constrained plus
  one free axis already yields N-shaped orders (e = 5 of 4!). The exact count is a DP over down-sets,
  O(down-sets x L).
- A merged leaf's relations are the union of its children's, so a death maps the new tree's down-sets
  injectively into the old tree's, and the whole-tree down-set count never grows under a death (0 of 6,575
  deaths across three enumerations).
- Measured sizes (current sampler, n 1000, p 5, 75 trees, 1/2/3 constrained predictors, 7,500 trees each):
  at most 9 leaves and 43 down-sets. The C++ count averages 0.14 us, against ~120 us per monotone tree step
  today (2.4 ms per 20-tree sweep, quadrature-bound).
- Stress (random guillotine partitions): two constrained axes, 116 us at 40 leaves. One constrained plus one
  free axis, 52 ms at 32 leaves and 0.1 s at 40 (850k down-sets).
- Redraw defect (blocking). [`MonotoneConstantGaussianLeaf::redrawAfterBirth`](../../src/bartcore/model.hpp)
  gives up after 100 tries and keeps mu[upper]. That slot is stale: freed node slots are recycled, and the
  resize zero-fills only new ones. It then draws the lower leaf on [aL, min(bL, stale)], which may be empty,
  and the clamp in [`MonotoneConstantGaussianLeaf::drawTruncatedNormal`](../../src/bartcore/model.hpp) then
  runs with lo > hi (undefined behaviour), leaving an infeasible state. The critique measured the fallback on
  about 0.5% of accepted births on flat data, 41% under a mildly decreasing partial residual, and 97% on a
  steeper one.
- Empty leaves. [`dbartsSampler$setPredictor`](../../man/dbartsSampler-class.Rd) with forceUpdate,
  updatePredictor, [`dbartsSampler$setData`](../../man/dbartsSampler-class.Rd) and
  [`dbartsSampler$setCutPoints`](../../man/dbartsSampler-class.Rd) can strand a member-empty leaf.
  [`MonotoneConstantGaussianLeaf::drawOneLeaf`](../../src/bartcore/model.hpp) pins it at 0 with no bound
  check, so the cone breaks and Z_T is no longer e / L!. This recurs in embedded use.
- Other paths. At creation R forces birth/death ([`resolveSamplerSpec`](../../R/spec.R)).
  [`xbart()`](../../R/xbart.R) refuses `monotone`, and [`ruleGibbsMove`](../../src/bartcore/moves.hpp)
  compiles out. Everything below was probed.
  - [`dbartsSampler$setControl`](../../man/dbartsSampler-class.Rd) installs a change/swap mixture.
    [`dbartsSampler$setModel`](../../man/dbartsSampler-class.Rd) with a model that lacks the monotone
    attribute does the same.
  - [`dbartsSampler$installTrees`](../../man/dbartsSampler-class.Rd) and
    [`dbartsSampler$setState`](../../man/dbartsSampler-class.Rd) accept an unconstrained donor's leaves. The
    fit then decreases by up to 1.4 along the constrained axis, and still by 0.48 after one sweep.
  - [`Chain::growForestFromRoot`](../../src/bartcore/chain.hpp) draws against the old tree's mu at reused
    node ids. It can also grow, as can installTrees and setState install, a tree of any size.
- Prior draw. [`MonotoneConstantGaussianLeaf::drawFromPriorForTree`](../../src/bartcore/model.hpp) samples the
  right law by rejection, but its 1e6-attempt cap fails ~6% of the time on a 9-leaf chain (1/9! = 2.8e-6).
- Unaffected: the leaf Gibbs sweep and the level-fibre shift (conditionals given T), and prediction.
  Monotone is new in 1.0-0, so NEWS gets no entry.
- Claims to reword: [4. Decision - marginal likelihood for the structure moves](../design/monotone.md#4-decision---marginal-likelihood-for-the-structure-moves),
  [11. Costs, risks, and confidence](../design/monotone.md#11-costs-risks-and-confidence),
  [Plan-vs-code note](../design/monotone.md#plan-vs-code-note), and dec-B16 in [decisions.md](../decisions.md).

## Decision

Question: should the exact count run under a work budget, and if so, measured how?

1. Whole-tree down-set budget, B = 2^22 (recommended). Trees whose down-set count exceeds B get zero prior
   mass, and that restriction is stated as part of the prior. The target stays exact for the stated prior, so
   dec-B14 holds.
   - Closure: deaths never raise the count, so the within-budget set is closed under deaths. Every
     within-budget tree still reaches the root through within-budget trees, and the restricted chain stays
     irreducible.
   - Prior mass removed: exceeding B takes 23 or more leaves, which carry 1.4e-18 under CGM(0.95, 2). No
     measured fit came within 5 orders of B.
   - Cost: at most ~0.5 s and ~200 MB for one proposal.
   - Handling: a birth over budget rejects. A current tree can never be over budget, because every entrance
     prunes or refuses (step 7), so no death needs an uncountable Z_T0.
   - Down-set keys are multi-word bitsets, so a long chain of any length is counted, not refused.
2. No budget. The documented prior stays verbatim, but one pathological tree can stall a sweep (exponential
   beyond ~40 leaves on a free axis) and exhaust memory. Every entrance must count whatever it installs.
3. A per-component leaf cap. Simpler to state, but free-axis deaths merge components, so it is not closed
   under deaths; it refuses cheap long chains and passes expensive wide stars. Not recommended.

Option 1 costs ~50 lines over option 2. Evidence that would change the call: a real fit whose count exceeds
~1 ms.

## Constraints

- Exact for the stated prior (dec-B14). Counts are held in double: exact to 2^53, with relative error under
  1e-12 beyond.
- Unconstrained samplers byte-identical: the new seam compiles out, like the three existing monotone seams.
  No dbarts.h change.
- Out of scope: change moves under the constraint, quadrature speed (TODO monotone-leaf-quadrature), and
  reconciling a chi k hyperprior with the truncated law.

## Steps

1. Order counter (engine, beside the monotone geometry). Down-set keys are multi-word bitsets.
   - Relation: build the tree's relation matrix with the adjacency test of `monotoneNeighborBounds`, and split
     it into components.
   - Count: run the layered down-set DP per component, keeping each layer's counts, and return the whole-tree
     log Z_T = sum over components of log e(C) - lgamma(|C| + 1), or over-budget (Decision).
   - A move changes only the component of the touched leaf in T0 and the component(s) holding the two children
     in T*. After a free-axis birth the children may share a component: take the distinct components.
2. Seam. Add an optional leaf concept, `logTreeNormalizer`, declared by the monotone leaf.
   [`birthOrDeathMove`](../../src/bartcore/moves.hpp) evaluates it in each state: before and after
   `tree.birth`, before and after `orphanChildren`. It adds log Z_T0 - log Z_T* to the log prior ratio (logged
   in the census's prior column), and an over-budget proposal rejects.
3. Drop the d terms: `priorMass` in `oneLeafLogMarginal`, `denom` in `twoLeafCoupledLogMarginal`.
4. Exact pair redraw in `redrawAfterBirth`:
   - Keep the capped rejection for the upper leaf; it is exact whenever it accepts.
   - When the cap is reached or acceptMax underflows, draw the upper leaf by inverting the CDF of its marginal,
     phi_U(u) [Phi_L(min(bL, u)) - Phi_L(aL)] on [max(aL, aU), bU]. Use the `coneProbability` quadrature's
     cumulative with safeguarded Newton, evaluate the density exactly, and work on the log scale in the tail.
   - Then draw the lower leaf on [aL, min(bL, mu_U)], which is never empty.
   - Mixing a capped exact rejection with an exact fallback is exact. Never read mu[upper] before writing it,
     and never call the clamp with lo > hi.
5. Empty leaves are leaves with no data. `drawOneLeaf`, both redraws and the prior draw give an empty leaf its
   prior truncated to its neighbour bounds instead of pinning it at 0, and
   [`monotoneTreeIsFeasible`](../../src/bartcore/model.hpp) checks empty leaves too.
   - The cone and Z_T = e / L! then hold in every state. Trees with an empty leaf still carry zero posterior
     mass through the veto, and the veto ranks move the chain out.
   - The cost is that a test row routed to an empty leaf predicts that leaf's draw instead of 0.
   - The alternative, refusing mutations that strand a leaf, breaks the embedded use the sampler exists for.
6. Exact prior leaf draw. Per component, draw a uniform linear extension by backward sampling on the kept DP
   counts: remove a maximal element x with probability e(D - x) / e(D). Draw |C| iid N(0, c scale / k),
   sort them, and assign them in extension order; isolated leaves draw alone. This replaces the rejection loop
   and its 1e6 cap in `drawFromPriorForTree`.
7. Reachability:
   - A bridge refusal after [`parseProposalProbs`](../../src/R_interface_bartcore.cpp), at creation and in
     [`bartcore_setModel`](../../src/R_interface_bartcore.cpp). It is keyed on the engine's leaf kind
     ([`LeafModelKind`](../../src/bartcore/model.hpp)), not on the incoming model, and refuses any nonzero swap,
     change, perturb or rule_gibbs share. The all-zero frozen mixture stays allowed.
   - setControl mirrors creation: a defaulted mixture is rewritten to birth/death silently, and an explicit
     non-birth/death one is refused.
   - growForestFromRoot reseeds mu to the all-zero feasible seed before its draw.
   - [`Chain::installForest`](../../src/bartcore/chain.hpp) reseeds infeasible leaves.
   - growForestFromRoot and installForest both prune an over-budget tree by deaths of its deepest nodes until it
     is within budget. This terminates because deaths never raise the count.
   - [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) (setState, copy, reload) refuses an infeasible
     or over-budget tree. With steps 4 and 5 every state the sampler produces passes this check.
   - [`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp) adds "within budget" to its rejection
     predicate.
8. tests/cpp:
   - The count against brute-force permutation counts on hand-built trees (1-3 axes, mixed directions, N, a
     star, a chain over 64 leaves) and 200 random trees.
   - Log Z against the enumeration's e / L!.
   - RNG-free ratio tests on A < B -> A < B1 < B2 (constrained split) and A < B -> {A < B1, B2} (free split)
     with pinned mu_A. Compare the move's log alpha with the closed form of the corrected statement in Context
     (it must differ from the current code's value).
   - The pair redraw on a design with P_post(lower <= upper) ~ 1e-4, where the old loop exhausts almost every
     time: every draw is feasible, and the upper leaf's draws match the quadrature CDF (KS).
   - The linear-extension draw is uniform over extensions (chi-square on a 5-leaf N-plus-chain), and the budget
     flag fires.
   - [`testMonotoneMarginal`](../../tests/cpp/test_model.cpp) loses its d_* = 1/2 check.
9. tinytest ([test-monotone.R](../../inst/tinytest/test-monotone.R)):
   - a setControl change mix errors, and a defaulted one is rewritten;
   - setModel without the monotone attribute cannot install a change mix;
   - installTrees from an unconstrained donor leaves the fit monotone at once, and setState of that state is
     refused;
   - growFromRoot plus one sweep is monotone;
   - after 2,000 sweeps on data decreasing along the constrained axis, copy() and setState round-trip;
   - setPredictor(forceUpdate = TRUE) stranding a leaf keeps the fit monotone.
   A statistical check does not fit: the most sensitive cheap functional sat at |z| 0.6-0.8 against the current
   move at 50k-100k draws.
10. Gate: [`buildDesign`](../../benchmarks/R/monotone-exact-enumeration.R), [`prototypeKeys`](../../benchmarks/R/monotone-exact-enumeration.R) (details under Verification).
    - The script enumerates every one-tree structure on cell grids with the CGM prior, weights each by
      exp(sum base) P_post(C) / Z_T, and computes P_post(C) by the same down-set recursion carrying
      Richardson-extrapolated trapezoid quadrature, self-checked to 2e-12 against the closed form.
    - It compares draws per root rule, which only the root-only tree can change, against each group's
      conditional law, using a Hotelling T^2 on batch-mean cell frequencies with its F reference. A group fails
      at p < 1e-4.
    - Designs, 10 rows per cell: c1, x1 constrained, 4 cells; c2, x1 and x2 constrained, 3x2; c3, x1
      constrained and x2 free, 3x2; cN, x1 constrained and x2 free, 2x3, where the two halves change level at
      different cuts, giving 26% N mass in the tested x1-rooted group.
    - The script replaces the planned part (c) of monotone-reference.R. It joins exact-gates.yaml's list in the
      fix commit, since it fails the current engine by design.
    - The general DP beyond these sizes rests on step 8's brute-force checks.
11. A monotone scenario in benchmarks/R/equivalence.R (x1 and x2 constrained, 20 trees), with the other 53
    scenarios bitwise. Its MANIFEST row names the enumeration gate as the ORACLE (P17).
12. Docs:
    - monotone.md sections 4, 9 and 11 and the Plan-vs-code note restate B' with the whole-tree normalizer, and
      record that mBART's d-normalized eq. 4.11 targets neither this prior nor the software's.
    - The prior statement names the budget if option 1 is taken.
    - dec-B16 is marked superseded in part by dec-B144.
    - Status lines and INDEX at landing, and the TODO item removed.

## Verification

- `cd tests/cpp && make && ./test_bartcore`: the count, ratio, redraw, extension-draw and budget checks pass.
- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'`.
- `R_LIBS=<lib> Rscript benchmarks/R/monotone-exact-enumeration.R quick`: every group passes. This takes about
  6 min; full mode runs 900k draws. The mutation run restores the d divisions and drops the Z term (then
  `touch` the header and reinstall), and must fail like the current engine.
  - Current engine, quick: c1 p 1.6e-24 and 6.5e-35 (two root rules); c2 3.8e-10 and 4.3e-5; c3 1.4e-15;
    cN 1.9e-29.
  - `prototype` mode (the corrected move in R, with an exact pair redraw): every group p >= 0.07.
  - `prototype-old` reproduces the current engine (c1 T^2 196 and 246, against the engine's 163 and 253).
  - Against an e(N) miscount of 10%, cN's noncentrality is ~80, so power is about 1.
  - `zcheck` mode verifies Z_T = e / L!.
- `R_LIBS=<lib> Rscript benchmarks/R/monotone-reference.R quick`: parts (a) and (b) still pass. Then run the
  whole exact-gates.yaml list with `quick`.
- `Rscript benchmarks/R/equivalence.R compare <current>`: 53 scenarios "identical draws (same RNG stream)" and
  no "max |z|". BCF and multinomial compare identical. The snapshot files carry no monotone fit.
- Release level: the SBC arm in [Monotone arm: design](sbc-family-tiers.md#monotone-arm-design). monotone-1
  (0.4 s per replicate) flags the current move and must pass; then the 20-tree arm (~85 min at R 200) must pass
  before admission.
- Speed: on a quiet machine, monotone sweep time at 20 trees, 1 and 2 constrained predictors, within 5% of
  today. bench-sampler compare unchanged on the unconstrained paths.
- `Rscript tools/check-doc-freshness.R .` passes.
