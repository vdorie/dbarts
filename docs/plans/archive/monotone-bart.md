# monotone-bart

agent: opus (all three commits: MCMC-kernel, constrained numerics, and the
  exact-posterior gate; routing is per-step below).
rng: C1 neutral; C2 neutral for the default path, posterior-changing only on
  the not-yet-reachable monotone leaf; C3 posterior-changing (monotone becomes
  user-reachable).
budget: 600-1000 engine lines across the three commits (design section 3), plus
  the R surface and gates.

## Goal

Per-variable monotone-increasing/decreasing constraints, `monotone = c(x1 =
"+", x3 = "-")` on dbarts()/bart2, enforced through a constrained constant-leaf
model: mBART's exact conditional structure-move target with honest low-
dimensional quadrature (approach B'), the default tree prior kept clean.
Unconstrained fits stay byte-identical at every commit.

## Context

docs/design/monotone.md is the spec (all VD resolutions inline; B', monotone-
first, base-speed a non-issue). Load-bearing anchors, re-verified 2026-07-18:
- Draw: `sampleParametersAndSetFits` constant branch, [[chain.hpp:2299-2319@bd098397]]; `mu`
  is a reference into the persistent `forest.muByTree[t]`, zeroed+refilled AFTER
  the moves ([[chain.hpp:2303-2304@bd098397]]), so it survives the prior sweep into this
  sweep's moves - the new work is keeping it valid THROUGH in-sweep births/deaths.
- Move seam: `metropolisJumpForTree` [[chain.hpp:702@bd098397]]; `logLikelihoodForBranch`
  [[moves.hpp:47-66@bd098397]] reads node suffstats and ZERO leaf params; birth/death score at
  [[moves.hpp:172-191@bd098397]], ratio shape [[moves.hpp:198-204@bd098397]]/258-264.
- Leaf concept: `ScalarLeafModel` per-node draw only, [[model.hpp:47-55@bd098397]]; posterior
  N(m_k,s_k^2) at [[model.hpp:127@bd098397]]; factory/instantiation seam [[facade.hpp:469@bd098397]].
- Neighbor bounds: `Tree::splitInterval` [[tree.hpp:313@bd098397]]. Truncated-normal primitive
  `ext_rng_simulateTruncatedNormalScale1` already exists, [[random.h:111@bd098397]] (REUSE).

## Constraints

- Unconstrained (all-zero direction) fits BYTE-IDENTICAL at every commit: the
  monotone leaf is a separate construction-time instantiation, if-constexpr dead-
  branch elided; NO runtime flag on the hot path (a truncated-normal-infinite-
  bounds path would perturb the Ziggurat stream). Equivalence trio bitwise +
  suite at every commit; bench-sampler neutral on the default path.
- v1 scope: constant leaves only; categorical predictors refuse; birth/death +
  single-site leaf Gibbs moves only (no change/swap); per-variable only (no
  convexity/multivariate shape). Out of scope: linear/gp leaf constraints.
- No dbarts.h change (internal leaf-model + bridge args only).

## Steps

1. C1 - the seam (2a), rng-neutral, opus. A tree-granularity draw concept method
   `drawParametersForTree` the constant leaf declares (default stays the per-node
   loop, byte-identical); the generic MoveStrategy M-access hook (score may read
   the current leaf vector M) that heteroscedastic later specializes; and the
   muByTree in-sweep lifetime change gated to the new path so the constant path
   keeps its post-move zero-and-refill. No monotone leaf yet - inert.
   Gate: equivalence trio bitwise, suite, tests/cpp seam-inert component test.
2. C2 - the monotone leaf (2b), opus. The constrained truncated-normal single-
   site Gibbs draw (via the C1 hook), `ConstrainedConjugateMove` (via the M-access
   hook) with neighbor geometry (splitInterval-derived [a,b] bounds), the B'
   constrained birth/death marginal (1-D/2-D honest quadrature, no grid), c^chi
   inflation (section 6), and empty-cone -> 0 feasibility. Reachable only from
   C++ tests, not R.
   Gate: equivalence trio STILL bitwise (default unchanged); tests/cpp - neighbor-
   geometry brute-force oracle, feasibility invariants under fuzzed moves,
   truncated marginal/draw vs quadrature, c-inflation variance property.
3. C3 - R surface + gates, opus. `monotone` arg parse on dbarts()/bart2 (accepts
   +/-, increasing/decreasing, +/-1, length-p); categorical refusal; k default ->
   fixed 2 (explicit chi refuses); proposal.probs error on explicit non-default;
   bridge wiring to pass directions. The two-part exact-posterior gate
   (benchmarks/R/monotone-reference.R: (a) one-cut enumerable <= 2-D, (b) fixed
   3-leaf double-bounded companion at 3-D); recovery test; Rd + NEWS.
   Gate: exact-posterior both parts to MC+quadrature error; equivalence trio
   bitwise; air format; full suite.

## Verification

Per commit: R CMD INSTALL --preclean; tests/cpp from make clean; equivalence
trio bitwise (equivalence-f494156, bcf-99205ee, multinomial-8c2b5fc); full
tinytest. C3 additionally: the exact-posterior gate passes both parts; a bench-
sampler compare confirms the default path is speed-neutral (quiet-machine grant).

## Landing

All three commits landed 2026-07-19, each independently gate-verified by the
orchestrator (equivalence trio bitwise at every commit - unconstrained fits are
byte-identical throughout):

1. df1b6e0 - the seams (C1). Two C++20 detection concepts (TreeDrawLeafModel,
   ParamScoringLeafModel) no shipped leaf declares, so both if-constexpr branches
   compile out; M-access threaded through MoveContext::leafParams (no move-signature
   churn) - the generic hook heteroscedastic reuses. Byte-neutral scaffolding.
2. 1944966 - the constrained monotone leaf (C2). Sequential single-site
   truncated-normal Gibbs draw over the surviving mu block, B' constrained
   birth/death marginal (closed-form + one adaptive quadrature, honest normalizer,
   empty-cone sentinel), splitInterval neighbor geometry, c-inflation. Four
   component-test gates (neighbor oracle, feasibility fuzz, marginal/draw vs
   quadrature, c-inflation variance). Reachable only from tests/cpp.
3. 3862dd0 - R surface + the exact-posterior gate (C3). monotone arg
   (named/worded/signed/positional), categorical + linear/gp refusal, fixed-k rule,
   birth/death-only forcing, bridge route with NO dbarts.h change. The two-part gate
   (one-cut enumerable + fixed 3-leaf double-bounded companion). Suite 3291.

Landing story, recorded because it is the methodology working: building the C3
exact-posterior gate exposed that C2's deterministic feasible-init redraw BIASED the
stationary posterior (the orchestrator had accepted it on an MCMC-validity argument -
the gate proved the argument wrong). C3 replaced it with the design's exact eq-4.17
conditional redraw (rejection cone draw on birth, truncated draw on death), which the
gate confirms. The bias was never shipped (the monotone path was unreachable from R
until C3). Then the gate's own part (b) reference showed E[mu2] at z=3.18: diagnosed
(reference grid sweep + 9-seed sampler pool) as O(d) quadrature discretization in the
COARSE reference, not sampler bias - the sampler is exact, and 90b6e41 tightened the
reference grid (0.0015 -> 1e-4) so every part (b) quantity sits at |z| <= 2.2.

Bench-sampler default-path speed-neutrality CONFIRMED 2026-07-20 (same-machine A/B,
since a naive compare against the cold b9d53c7 CSV was swamped by ~7% cross-session
machine drift). Method: rebuild b9d53c7 and HEAD (which carries monotone +
heteroscedastic + the data-store-residuals batch, all landed after b9d53c7) under
identical conditions; HEAD vs a freshly-recorded b9d53c7 baseline geomean 0.981 (HEAD
marginally faster), every arm inside the +/- 6% same-code noise floor. The control
that proves it was drift, not code: the IDENTICAL b9d53c7 binary re-measured now ran
+7.4% geomean (up to +18% on the tiny setPredictor-reject arm) versus its own 07-17
cold CSV. The equivalence trio already proved the draws bitwise-identical; this rules
out an incidental (layout/branch) slowdown, and none appeared. v2 doors: change/swap
moves under the
constraint (a >2-leaf branch marginal beyond the product form), linear/gp leaf
constraints, convexity/multivariate shape - all out of v1 scope per the design.
