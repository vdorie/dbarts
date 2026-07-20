# Per-variable monotone constraints (mBART): design

Status: LANDED 2026-07-19 (ee4ca79). Plan: docs/plans/monotone-bart.md (this is its
step 1). Users declare a monotone-increasing or -decreasing relationship for any
subset of predictors; each tree of the forest is constrained to be monotone in
those predictors, so the sum-of-trees fit is monotone (Chipman, George,
McCulloch, Shively 2022, Bayesian Analysis 17(2):515-544; arXiv:1612.01619 -
"mBART"). The constraint rides as a per-column sign flag on the model spec,
enforced through the tree-granularity leaf-parameter draw the core provisioned
(docs/design/core-generalization.md:141) plus a constrained branch score for the
birth/death moves. Constant leaves only in v1; categorical predictors refuse the
constraint. rng class: posterior-changing when any constraint is declared, and
BYTE-neutral (no rng-stream perturbation, no hot-path cost) when none is
(section 8).

There is no maintained implementation to inherit: mBART never shipped a CRAN
package (the reference C++ is bitbucket.org/remcc/mbart), and stochtree states it
lacks shape constraints. The cross-ecosystem user expectation is set by XGBoost /
LightGBM `monotone_constraints` and sklearn's `monotonic_cst` - a per-feature
{-1, 0, +1} vector; follow-on BART literature is active through 2025 (Probit
Monotone BART, arXiv:2509.00263). This note tracks the paper's algorithm and
flags where v1 deviates (honest normalizers, section 4; exact-posterior gate,
section 9).

## 1. The mBART model and its reduction to our engine

**Per-tree monotonicity (paper Section 2).** A sum-of-trees function is monotone
in a predictor set S whenever each component tree is monotone in S, so the
constraint is imposed tree by tree. A single tree partitions x into rectangular
terminal regions R_k = {x : x_i in [L_ik, U_ik)}. Along a constrained coordinate
i, region R_k is an ABOVE-neighbor of R_k* when their boxes adjoin with R_k on
the larger-x side (L_ik = U_ik*) and their projections overlap on every other
coordinate; a BELOW-neighbor adjoins on the smaller-x side (U_ik = L_ik*). A tree
is nondecreasing in x_i iff every terminal mu is

    mu_k <= min over above-neighbors of mu,   and
    mu_k >= max over below-neighbors of mu

(paper Section 2, conditions (a)/(b)); nonincreasing flips the inequality, which
is equivalent to negating x_i. This coupling among neighboring leaf values is the
whole of the constraint - it is a cone C(T) = {M : g(.;T,M) monotone} in leaf
space, an intersection of pairwise half-spaces mu_below <= mu_above.

**The engine already has the two phases mBART needs.** mBART is MH-within-Gibbs,
backfitting one tree at a time against partial residuals (paper Sections 4.1-4.2),
with exactly two steps per tree that the constraint touches:

1. a birth/death Metropolis move on the tree structure, and
2. a draw of the leaf parameters given the structure.

Our chain runs these as `metropolisJumpForTree` (chain.hpp:702) then
`sampleParametersAndSetFits` (chain.hpp:2270). The constraint modifies each: the
move's branch score becomes a constrained (truncated) local marginal (section 4),
and the parameter draw becomes a sequential truncated-normal Gibbs sweep coupling
neighboring leaves (section 4). No new sweep structure is introduced; the two
seams that exist are re-pointed for the constrained forest.

**Box geometry is already computed.** Each leaf's [L_ik, U_ik] along an ordinal
column is `Tree::splitInterval` (tree.hpp:313), the ancestor-constrained cut
interval used for availability (`hasAnyAvailableVariable`, tree.hpp:437,
moves.hpp:81). Neighbor determination for a constrained axis walks the tree's
bottom nodes and, for each ordered pair, tests boundary adjacency in the
constrained coordinate and interval overlap in every OTHER coordinate the two
leaves split on - so the overlap test spans the distinct split variables on the
two leaves' ancestor paths, not just |S|: O(b^2 * d) per tree for b leaves and d
distinct ancestor-path split variables. Still negligible because BART trees are
small (the paper leans on this same fact, Section 3.3). This is the one genuinely
new piece of tree geometry and the component-test target of section 9.

## 2. Decision - constraint representation and surface

**Per-column sign flag on the model spec.** Resolve `monotone` into a borrowed
per-predictor direction vector, values in {-1, 0, +1}, carried to the engine as a
new `SamplerOptions` field beside the other per-column designations
(`columnTypes` chain.hpp:51, `leafCovariateColumns` chain.hpp:82) and consumed
once at construction. A nonzero entry on any column selects the constrained
instantiation at the factory (`createSampler`, facade.hpp:469); an all-zero
vector is treated as null and selects the existing constant-leaf path unchanged
(section 8). This mirrors how the linear-leaf designation rides and dispatches,
and keeps the constraint out of the data container - the same data can be fit
constrained or not.

**Surface (recommend `monotone`, a named vector).**

    dbarts(y ~ ., data, monotone = c(x1 = "+", x3 = "-"))
    bart2(x, y, monotone = c(x1 = "+", x3 = "-"))

Accept the sign glyphs "+"/"-", the words "increasing"/"decreasing", and the
integers +1/-1; unnamed integer/character vectors of length p are accepted as a
positional {-1,0,+1} spec (the XGBoost/LightGBM form). Names resolve against the
model-matrix column names after expansion; a name that does not resolve, or a
sign on a categorical column, is a hard error at spec time (see below). The
argument also rides `dbartsData`/`dbartsSampler` and `setModel` refuses to change
it after construction (the constraint is structural, like the tree/chain counts
setControl freezes, core-generalization.md:739).

Naming rationale (fork 6, a collision trap): NEVER `mbart`. CRAN's BART package
exports `mbart`/`mbart2` for MULTINOMIAL BART; a `dbarts::mbart` or a
`family = "mbart"` would read as multinomial to anyone who knows that package.
There is no new top-level function and no new `family` value - monotonicity is an
orthogonal argument that composes with gaussian and (probit/logistic) binary
responses, not a response family. `monotone` is short, matches the ecosystem's
`monotone_constraints`/`monotonic_cst`, and cannot be confused with a BART
function name. Alternatives considered: `monotone.increasing = c("x1")` +
`monotone.decreasing = c("x3")` (two args, no sign vocabulary to parse, but two
places to keep consistent and verbose for many columns) - rejected as less dense;
`monotonic_cst`/`monotone_constraints` (verbatim sklearn/XGBoost) - rejected for
not matching dbarts's `lower.case.dotted` R style.

**Categorical predictors refuse.** Monotonicity is undefined on an unordered
factor (category codes 0..K-1 carry no order; core-generalization.md:621). A
nonzero direction on a categorical column is refused at spec time in R and,
defensively, at the factory (the categorical-refusal already guarding leaf
covariates, facade.hpp:456). ORDERED factors and numeric columns are stored as
ordinal `ColumnType` (core-generalization.md:473) and DO accept the constraint -
their codes are ordered, so `splitInterval` and the neighbor test are meaningful.
So "categorical refuses" is precisely "unordered factors refuse; ordinal columns
and numerics are eligible," the clean `is.ordered()`-adjacent split.

## 3. Decision - the constrained leaf-parameter draw

RESOLVE as: sequential single-site Gibbs of truncated normals over the tree's
leaves, one node at a time, each conditional on the current values of its
neighbors - mBART's own scheme (paper Section 4.3, "draws of a single mu
component given T and all the remaining elements of M").

For leaf k with unconstrained full-conditional N(m_k, s_k^2) (the standard
constant-leaf posterior, model.hpp:127, from the node's (sum w, sum wz)
sufficient statistic), the constrained draw is that normal truncated to

    [ max over below-neighbors of current mu,  min over above-neighbors of current mu ]

along each constrained axis, intersected across axes. Because the previous sweep
left every mu feasible and single-site updates keep the others fixed, this
interval is always non-empty (it contains the leaf's own current value), so the
draw never fails. Sweeping all leaves in a fixed order is a valid Gibbs kernel on
the truncated stationary law; a random scan is the documented alternative if the
fixed order mixes poorly.

**Seam - what actually has to change (the "hook" is aspirational).** The
provision (core-generalization.md:141) reads as though a tree-granularity leaf
draw seam already exists; it does NOT. Today the draw is a hardcoded per-node
loop that rebuilds mu from zero every sweep (`mu.assign(tree.nodes.size(), 0.0)`,
chain.hpp:2304, then an independent `drawFromPosteriorForNode` per bottom node,
chain.hpp:2305-2319), the leaf-model concept exposes ONLY that per-node scalar
draw (`ScalarLeafModel`, model.hpp:47-55), and the moves are leaf-templated free
functions that read the node sufficient statistics and ZERO leaf parameters
(`logLikelihoodForBranch`, moves.hpp:47-66). Three real changes follow, none a
mere re-pointing:

1. A new concept method - a tree-granularity draw
   `drawParametersForTree(rng, tree, k, sigma2, muOut)` the constant monotone leaf
   provides and `sampleParametersAndSetFits` calls once per tree, in place of the
   per-node loop at chain.hpp:2305, when the leaf model declares it. Independent
   per-leaf stays the default for every other model, byte-identical (section 8).
2. Leaf values must stay VALID through the move phase. The move score needs
   mu_same, the frozen neighbor values, DURING `metropolisJumpForTree`
   (chain.hpp:702). Those values already survive there: `sampleParametersAndSetFits`
   runs AFTER the moves and only then zeroes and refills `muByTree[t]`
   (`mu.assign(tree.nodes.size(), 0.0)`, chain.hpp:2303-2304 - `mu` is a reference
   INTO the persistent `forest.muByTree[t]`), so during the moves the vector still
   holds the previous sweep's draw. The new work is keeping it consistent ACROSS
   in-sweep structural changes: an accepted birth adds two nodes and a death removes
   one, so `muByTree[t]` (sized to the pre-move tree) must be resized/patched at each
   acceptance rather than rebuilt once at the end, with the touched leaves' redraw
   written in (step 3) so later moves in the same sweep read consistent neighbors.
   The constant path keeps its post-move zero-and-refill unchanged; only the
   monotone path maintains muByTree incrementally through the moves - a change to
   the leaf store's in-sweep lifetime, not the creation of persistence that already
   exists.
3. The moves gain a leaf read/write responsibility they entirely lack today.
   `ConstrainedConjugateMove` must READ muByTree for the frozen neighbors when it
   scores a proposal, and on an accepted birth WRITE the eq.-4.17 redraw of the
   touched leaves back into muByTree so later moves in the same sweep see
   consistent neighbors - new plumbing between the move machinery and the
   leaf-parameter store, not a tweak to the existing conjugate score.

The alternative seat - reusing `FunctionLeafModel::beginTreeDraw` (model.hpp:82-85)
and carrying mu as a degenerate function leaf - is rejected: it drags the
function-leaf test-cache and fits-are-parameters machinery (chain.hpp:2277-2298)
onto a constant leaf for no benefit and muddies the chi-k accounting.

**Budget.** The plan front-matter's "~500 lines" (docs/plans/monotone-bart.md) no
longer holds once (2) and (3) are counted. Persisting leaf state into the move
phase, a new constrained MoveStrategy that reads and writes it, the neighbor
geometry (section 1), the tree-granularity draw, and c-inflation (section 6) are
plausibly 600-1000 lines of engine C++ before the R surface, bridge, and gates.
Recommend the plan split step 2 into 2a (leaf-state persistence + the
tree-granularity draw refactor - a mechanical but wide change) and 2b (the
constrained move + neighbor geometry + c-inflation), and re-budget each.

**Primitive (already in tree).** The Gibbs draw needs a doubly-truncated normal
with arbitrary mean and sd. `ext_rng_simulateTruncatedNormalScale1(rng, mean,
lower, upper)` already exists (random.h:111, landed by the ordinal arc:
inverse-CDF in the bulk, Robert 1995 rejection when the tail-probability gap
underflows). Monotone REUSES it, scaling the interval and the draw by the
posterior sd s_k, so the numerically hard tail-interval work is shared and was
tested once (deep-tail coverage already lives in ordinal's primitive test). No new
rng primitive is needed.

## 4. Decision - marginal likelihood for the structure moves

This is the load-bearing fork. Under the constraint the leaf values are DEPENDENT
given T, so the birth/death score is no longer a sum of independent per-leaf
Gaussian marginals (the current `logLikelihoodForBranch`, moves.hpp:47-66). The
options, and what mBART actually does:

- **A. Fully-collapsed exact marginal.** Integrate ALL leaf values over the full
  cone C(T): p(r | T) = integral over C(T) of prod_i N(r_i; mu_{node(i)}) *
  truncated-product-normal dM. This is a Gaussian orthant/box probability in b
  dimensions - no closed form and expensive to approximate honestly for every
  proposal. Rejected; also NOT what mBART does.

- **B. mBART's conditional (partially-collapsed) marginal (paper eq. 4.11).**
  Condition on mu_same, the leaf values NOT created or destroyed by the move, and
  integrate out only the one or two leaves the move touches, over their cone
  given the frozen neighbors. For a birth splitting leaf mu_0 into (mu_L, mu_R),

        p(r_new | T*, mu_same) = integral over C of
            N(r_L | mu_L) N(r_R | mu_R) * phi(mu_L) phi(mu_R) / d_*  d mu,

    with C = {a <= mu_L <= mu_R <= b} (paper eq. 4.12-4.13), a,b the outer bounds
    from the frozen neighbors. The acceptance is (paper eq. 4.11)

        alpha = min{ [q(T*->T0)/q(T0->T*)] *
                     [p(r_new|T*,mu_same) / p(r_old|T0,mu_same)] *
                     [p(T*)/p(T0)], 1 },

    the SAME shape as the unconstrained ratio the engine computes today
    (moves.hpp:198-204 birth, moves.hpp:258-264 death). This is a valid
    MH-within-Gibbs block update of (T, touched-mu) with the other leaves held
    fixed, then a redraw of the touched mu from the truncated posterior over C
    (eq. 4.17). It targets the EXACT constrained posterior - the conditioning is
    not an approximation. That exactness is the paper's own derivation (Section
    4.2, going from eq. 4.10 to eq. 4.11): the standard CGM98 reverse-move
    probability q(T*->T0) supplies the transition ratio, and rsame gives "the
    same multiplicative contribution to the top and bottom of the acceptance
    ratio so that it cancels out" conditional on mu_same. Taken here on the
    paper's word (the reverse-move-normalizer step is not re-derived in this
    note); the exact-posterior gate (section 9) is what checks it empirically.

- **B'. B with HONEST normalizers (recommend).** mBART's IMPLEMENTATION (paper
  Section 4.3) then makes two approximations for speed that we should NOT inherit,
  both stated verbatim there (Greek and tildes transliterated to ASCII):
    (i) it grids each integral - "we approximate them by summing over a grid of
    (mu_L, mu_R) values" (Section 4.3, near eq. 4.14), reporting "a 20x20 grid
    size (which yields excellent results for each two-dimensional numerical
    approximation)";
    (ii) it drops the constrained normalizers d_* (birth, bivariate) and d_0
    (death, univariate) - "we reduce the computational burden by letting d_* and
    d_0 equal one and then compensating for this omission with an adjustment of
    our T prior", then "We compensate for the omission of d_* and d_0 by letting
    alpha = .25 and beta = .8 rather than using standard BART default values of
    alpha = .95 and beta = 2" (Section 4.3; alpha/beta are the tree prior's
    base/power).
  That couples the constraint to the tree prior and makes the sampler target an
  approximation, which the exact-posterior gate (section 9) would then fail. We
  keep the normalizers honest instead:
    * The touched-leaf marginal factorizes cleanly. A birth NOT on a constrained
      axis leaves mu_L, mu_R each independently truncated to their own [a,b]
      (their only shared neighbor bound), so the marginal is a PRODUCT of two 1-D
      Gaussian-times-truncated-Gaussian integrals - closed form in Phi. A birth
      ON a constrained axis couples them only through mu_L <= mu_R, a single
      linear inequality; the marginal is one bivariate-normal box probability,
      one 1-D adaptive quadrature (not a 20x20 grid). Death is always 1-D
      (eq. 4.18), closed form.
    * The redraw of the touched leaves is EXACT, not gridded: the constrained-axis
      birth draws (mu_R, then mu_L | mu_R) as two sequential 1-D truncated normals
      honoring mu_L <= mu_R; every other case is independent 1-D truncated
      normals. No mesh, no grid artifacts.
  Honest d keeps the DEFAULT tree prior (base=0.95, power=2, chain.hpp:39)
  unchanged, so the constraint does not smuggle a different structural prior, and
  the sampler targets the exact posterior up to 1-D quadrature error.

- **C. Fully non-conjugate MH (propose-and-draw).** Draw the touched leaves from
  the constrained prior as part of the proposal and accept on the joint
  likelihood without integrating them out - the reversible-jump / prior-grown
  machinery gp-leaves.md and the open Phase 6 (core-generalization.md:780) build
  for non-integrable leaves. Correct, and it needs no truncated marginal, but it
  gives up the Rao-Blackwellization B/B' keep (the touched leaves are integrated,
  not sampled, in the accept decision), so it mixes worse - the standard argument
  for BART's collapsed moves. Rejected for v1 as the primary path.

Recommendation: **B'**. It is mBART's exact target with our numerics substituted
for its grid, maps onto the existing move seam (replace the two
`logIntegratedLikelihoodForNode` evaluations for the touched leaves at
moves.hpp:172-191 with the constrained joint marginal that reads the frozen
neighbor mu), and keeps the tree prior clean for the gate.

RESOLVED (VD, 2026-07-18): **B'**, and proceed with the monotone arc now -
before heteroscedastic - building the generic MoveStrategy M-access hook that
heteroscedastic later specializes. The design was blind-critiqued (REJECT: gate
closure, k-default, the aspirational-hook framing), amended, and re-critiqued
(ACCEPT). Base-speed confirmed a non-issue: the unconstrained constant-leaf path
is byte-identical by construction-time instantiation + if-constexpr (section 8),
so the abstraction costs no default-path speed - the work is the 600-1000 line
engine leg (2a/2b/C3), not a tradeoff.

**Shared-machinery sequencing (flagged per the plan).** The current conjugate
MoveStrategy reads NO leaf parameters - it integrates them all out. Both B' and
option C need a MoveStrategy that (i) reads the current leaf-parameter vector M
during scoring and (ii) computes a bespoke local marginal rather than the
closed-form sum. That "score reads M + custom local marginal" seam is ALSO what
the queued heteroscedastic arc's non-conjugate MoveStrategy needs
(core-generalization.md:780, gp-leaves.md). Consideration: monotone does NOT have
to wait for the full non-conjugate arc - B' still integrates conjugately, just
over a truncated region, so it is a distinct, simpler MoveStrategy
(`ConstrainedConjugateMove`, compile-time-selected for the monotone leaf like the
plain conjugate move requires `IntegrableLeaf`). But whichever of {monotone,
heteroscedastic} lands first should introduce the M-access hook GENERICALLY so
the other specializes it (constrained-conjugate marginal here; Laplace /
pseudo-marginal there) rather than forking two parallel move-scoring paths. If
sequencing is free, landing the non-conjugate seam first and making monotone a
specialization is the lower-duplication order.

## 5. Decision - proposal feasibility (a birth can violate monotonicity)

The paper is silent on empty cones (it never states a hard-rejection rule), but
the structure makes the answer clean, and it differs by move:

- **Birth NEVER produces an empty cone.** Splitting a feasible leaf mu_0 into
  (mu_L, mu_R) yields C = {a <= mu_L <= mu_R <= b} where a,b are mu_0's own frozen
  outer bounds; since a <= mu_0 <= b held, (mu_L, mu_R) = (mu_0, mu_0) is in C, so
  C is non-empty. Feasibility is therefore NOT a hard gate on births - the
  constraint only reshapes the marginal (B') and the redraw. This is an invariant
  worth asserting in a component test: every accepted birth leaves a feasible
  tree.
- **Death** re-merges two ordered children into one leaf whose cone is the parent
  bound [a,b], always non-empty. Never infeasible.
- **Change / swap CAN empty a cone** (reordering splits can force some leaf's
  max-below above its min-above). Handled with NO special-casing: the constrained
  marginal integrates the truncated-normal product over an empty region to 0, so
  p(r_new | T*, mu_same) = 0 and alpha = 0 - the move is rejected by the ordinary
  acceptance test, the same way the empty-leaf veto (-HUGE_VAL, moves.hpp:62)
  already rejects unoccupied leaves. No infeasible state can be accepted.

**v1 move set (scope decision).** Restrict the constrained forest to birth/death
plus the single-site leaf Gibbs - exactly mBART's examples ("in all our examples
we use birth/death moves and draws of a single mu component," paper Section 4.3).
Set the constrained forest's move mix to birth/death only
(birthOrDeathProbability = 1, chain.hpp:40), a legitimate difference from the
unconstrained default mix (0.5/0.4/0.1 birth-death/change/swap). This overrides
the user-facing `proposal.probs` (default c(birth_death = 0.5, swap = 0.1,
change = 0.4, birth = 0.5), R/dbarts.R:339): resolve the clash by ERRORING at spec
time when `monotone` meets an EXPLICIT non-default `proposal.probs` (the user
asked for swap/change the constrained sampler cannot honor), and forcing
birth/death only, silently, when `proposal.probs` is left at its default.
Rationale: a
constrained CHANGE move whose children are both terminal is still <= 2-D and
supportable (the paper says so), but the general change/swap touches more than
two leaves and needs a > 2-D cone integral, which is where mBART itself stopped.
Cost/risk: change moves aid mixing (docs/design/change-move-balance.md); dropping
them may slow the constrained chain. v1 accepts this; the children-terminal
change move is the recorded v2 extension (still 1-2 D, so B' covers it).

## 6. Decision - prior calibration under truncation

mBART keeps the CGM10 leaf prior mu ~ N(0, sigma_mu^2) with
sigma_mu = 0.5/(k sqrt(m)) (Y scaled to [-0.5, 0.5], k=2 default, m trees) - which
is EXACTLY dbarts's constant-leaf prior sd = scale/k with scale =
nodeScale/sqrt(numTrees), nodeScale = 0.5 (model.hpp:103, chain.hpp:38,372). The
one change (paper Section 3.3, eq. 3.6): a leaf that IS constrained uses an
inflated variance c^2 sigma_mu^2 with

    c^2 = pi / (pi - 1) ~= 1.4669   (c ~= 1.2112),

so that the post-truncation (skew-normal) MARGINAL variance of a singly-
constrained leaf equals the unconstrained sigma_mu^2. The variance-matching
derivation lives in paper eq. 3.7-3.9 (the constrained pair's marginals are
skew-normal with the stated means and common variance sigma_mu^2 at c^2 =
pi/(pi-1)) and Azzalini (1985); it is taken here on the paper's word, not
re-derived, and its consequence is checked by the c-inflation component test
(section 9). This balances the prior's effect across constrained and
unconstrained predictors and preserves the induced N(m*mu_mu, m*sigma_mu^2) prior
on E(Y|x) (paper Section 3.3, CLT argument: constrained means carry pairwise
offsetting biases).

Recommendation: adopt the c-inflation. It is per-leaf and tree-dependent - a leaf
uses c*scale iff it currently has at least one neighbor along a constrained axis,
determined during the same tree walk that finds the truncation bounds, so it is
free. Following mBART, IGNORE multiplicity (a leaf constrained on two sides gets c
once, not c^2); the paper justifies this because trees are small and the
correction is second-order (Section 3.3). Both the constrained marginal (section
4) and the Gibbs draw (section 3) use c*scale for constrained leaves.

**Sub-decision - the chi-k hyperprior.** dbarts optionally samples k from the
accumulated sum of standardized squared leaf values (`forest.updateK`,
chain.hpp:2287,2315). Under truncation and a per-leaf-variable prior scale, the
"standardized" square is no longer param/(scale/k), so feeding truncated draws
into the chi-k update biases k. v1 uses FIXED k under `monotone` (mBART itself
uses fixed k=2), but must not turn that into a default-fit error: the binary
default k IS chi(1.5, 2.0), not a user choice (resolveNodeHyperprior k=NULL ->
chi(1.5, 2.0) for binary, R/model.R:367; .kDefault, R/bart.R:283; bart2 k=NULL
resolves the same), so a plain bart2(..., monotone=..., family="probit") supplies
a chi hyperprior by default. Resolution: under `monotone`, an UNSUPPLIED k
resolves to fixed k = 2 (the continuous default and mBART's value) for BOTH
continuous and binary responses - silently, since the user requested no
hyperprior; only an EXPLICIT chi request (k = chi(...)) refuses, with an
informative message pointing at fixed k. Reconciling the truncated leaf law with
the chi-k posterior is a documented follow-up, not v1 scope.

## 7. Seams summary

- **R:** `monotone` arg on dbarts/bart2/dbartsData; resolves to a per-column
  {-1,0,+1} vector on the model spec; refuses categorical columns and post-hoc
  change (`setModel`). Probability transforms and prediction are unaffected -
  monotonicity constrains f, not the link.
- **Bridge (src/R_interface_bartcore.cpp):** carry the direction vector into a new
  borrowed `SamplerOptions::monotoneDirections` (length numPredictors), consumed
  at construction like `columnTypes`. No dbarts.h / ABI change in v1 (the
  robust-errors/ordinal precedent) - monotone is reachable only through the R
  surface and internal bridge, so no LinkingTo consumer sees it.
- **Factory (facade.hpp:469):** any nonzero direction selects the monotone
  constant-leaf instantiation (`MonotoneConstantGaussianLeaf` +
  `ConstrainedConjugateMove`); all-zero selects the existing constant-leaf path
  verbatim (section 8).
- **Leaf store (chain.hpp):** muByTree must PERSIST across the move phase - stop
  the zero-rebuild at chain.hpp:2304 so the moves can read frozen neighbors - and
  the tree-granularity draw replaces the per-node loop at chain.hpp:2305 (section
  3, changes 1-2). This is the widest mechanical change and the reason for the
  re-budget.
- **Leaf model (model.hpp):** a constant leaf with a tree-granularity Gibbs draw
  (section 3) and per-leaf c-inflation (section 6); holds a pointer to the
  direction vector and computes neighbor bounds via `Tree::splitInterval`
  (tree.hpp:313).
- **Move (moves.hpp):** `ConstrainedConjugateMove` replaces the per-leaf marginal
  sum (moves.hpp:54-65) with the conditional truncated joint marginal of the
  touched leaves (section 4), READING frozen neighbor mu from muByTree and WRITING
  the accepted eq.-4.17 redraw back (section 3, change 3).

## 8. Bitwise neutrality for unconstrained fits

Binding requirement (rng class: neutral when `monotone` is absent). When no
constraint is declared the direction vector is null and the factory selects the
UNCHANGED `ConstantGaussianLeaf` + conjugate MoveStrategy instantiation - not a
constrained type with the constraint switched off. This is the neutrality
guarantee: the unconstrained code path is literally the same object code, so it
draws the identical rng stream (a truncated-normal-with-infinite-bounds does NOT
reproduce `ext_rng_simulateStandardNormal`'s Ziggurat stream, so a runtime "off"
flag inside one leaf type would perturb draws - construction-time type selection
avoids that entirely). No per-observation cost is added anywhere: the neighbor
walk, truncation, and constrained marginal live only in the monotone
instantiation. Binary size grows by one instantiation (constant x monotone),
bounded like the existing leaf-model matrix (core-generalization.md:118).

Gates that prove neutrality (section 9): the equivalence baselines reproduce
byte-identically with NO re-record (monotone adds no draw to any existing family
or leaf model), and bench-sampler shows no regression on the unconstrained
constant-leaf paths.

## 9. Gates

**Exact-posterior gate (the plan's required gate).** TWO pieces, because one toy
cannot both enumerate tree space cheaply AND exercise a leaf bounded on both
sides - and the dimension must be stated honestly. The established pattern is
benchmarks/R/logistic-reference.R, categorical-exact.R.

(a) Enumerable ONE-cut gate. One constrained predictor with EXACTLY ONE available
cut, Gaussian response, sigma FIXED. The predictor then admits only two structures
- root, or the single split into an ordered pair {mu_L <= mu_R} - and the sampler
itself never grows past two leaves, so every move marginal AND the reference stay
at most 2-D. (With more than one cut the chain does grow past two leaves - depth-1
growth is live, moves.hpp:186-188 - and a b-leaf structure's reference marginal
becomes a b-D constrained integral, NOT the 1-D/2-D quadrature; the one-cut
restriction is what removes that closure gap.) For each structure the marginal
likelihood is the constrained integral

    integral over C of  prod_i N(r_i; mu_{node(i)}, sigma^2)
                        * prod_leaf N(mu_leaf; 0, (c^chi * scale/k)^2)  dM,

C the monotone cone (root: none, 1-D over mu_0; split: {mu_L <= mu_R}, 2-D), c^chi
inflating constrained leaves (section 6), an empty leaf contributing a factor of 1
(its prior integrates to unity). This is mBART's eq. 4.12 integral computed
HONESTLY (including the d normalizers, section 4) - only the honest version (B')
matches the quadrature, so this doubles as the gate on B' vs mBART's d=1
approximation. Each structure's posterior weight is its tree prior times this
marginal, renormalized; match the sampler's posterior mean of the identified
quantities (fitted values, structure probabilities) to the quadrature, to Monte
Carlo error. This closes exactly and validates the tree-structure posterior, the
move score, and the coupled TWO-leaf draw - but its two leaves are each
single-bounded (mu_L only above by mu_R, mu_R only below by mu_L), so it does NOT
touch a leaf bounded on both sides, the highest-risk geometry (section 1).

(b) Fixed-structure multi-leaf companion. A double-bounded interior needs >= 3
ordered leaves (mu_1 <= mu_2 <= mu_3, mu_2 bounded below by mu_1 AND above by
mu_3), whose joint leaf posterior is a genuine 3-D constrained integral. Rather
than enumerate 3-leaf tree space, PIN a hand-built 3-leaf monotone tree (two cuts
on the constrained axis), disable structure moves, and run ONLY the single-site
truncated-normal leaf Gibbs (section 3); match its stationary leaf-parameter
posterior (means and pairwise covariances) against a direct 3-D
constrained-quadrature of the joint N(.; posterior) x monotone-cone - the
posterior here carries the same c^chi-inflated leaf prior (section 6) the draw
uses, entering through each leaf's full conditional. This isolates and validates
the double-bounded interior draw at its true 3-D dimension without paying the
tree-space enumeration.

Tolerances bound MC plus quadrature error and are never widened to pass. Together
(a) and (b) cover the two axes one toy could not: (a) the full tree-structure
posterior and move score at <= 2-D, (b) the double-bounded interior leaf draw at
3-D. The alternative - a single gate enumerating the full chain with up to
(C+1)-D constrained quadrature for C cuts - is recorded but rejected as costlier
for no extra coverage: its 3-leaf structures need exactly the 3-D quadrature (b)
already does, atop the tree-space sum (a) already does.

**Component tests (RNG-free where possible).**
- Neighbor geometry: on hand-built trees over 2-3 constrained axes, the
  above/below-neighbor sets and per-leaf [a,b] bounds match a brute-force box
  adjacency oracle. This is the highest-risk new code (section 1).
- Feasibility invariants under every move: after any accepted birth/death/leaf
  draw the tree is monotone (every leaf within its neighbor bounds); a change/swap
  that empties a cone scores 0 and is rejected (section 5). Assert on a fuzzed
  move sequence.
- Truncated marginal and draw: the 1-D and 2-D constrained marginals match direct
  quadrature; the constrained-axis birth redraw (mu_R then mu_L|mu_R) has support
  {a<=mu_L<=mu_R<=b} and matches the analytic truncated moments, including a
  deep-tail interval exercising the Robert-rejection fallback.
- c-inflation: a singly-constrained leaf's post-truncation marginal variance
  equals the unconstrained sigma_mu^2 (the eq. 3.7-3.9 property).

**Recovery test (the plan's step 3).** Simulated data from a monotone truth at
moderate n: the constrained fit recovers it with tighter intervals than
unconstrained BART, and a non-monotone truth fit under the constraint shows the
EXPECTED bias (flattening against the constraint), not a crash or a
degenerate fit - the family-level smoke beyond the exact gate.

**Neutrality trail.** Re-run benchmarks/R/equivalence.R compare and expect
IDENTICAL draws for every existing scenario (no re-record); bench-sampler compare
on a quiet machine shows the unconstrained constant-leaf paths unchanged (the
hook must not slow the default, plan Verification). A NEW equivalence scenario
records the constrained channels (fitted values, a monotone-recovery summary) as
the posterior-changing baseline for this arc.

## 10. Out of scope, and the doors

- **Non-constant leaves.** Constant leaves only in v1. Linear or GP leaves under a
  monotone constraint are a different problem (a linear leaf is monotone iff its
  slope has the right sign - a coefficient constraint, not a neighbor-ordering
  one; linear-leaves.md, gp-leaves.md). Refused at spec time when `monotone`
  meets a leaf-covariate designation.
- **Convexity / general shape constraints.** Out of scope; monotonicity is a
  first-difference sign constraint, convexity a second-difference one, needing a
  different neighbor calculus. Not provisioned here.
- **Joint / multivariate monotonicity** beyond per-variable (e.g. monotone in
  x1+x2 jointly). Out of scope; the constraint vocabulary is per-column signs.
- **General change/swap under constraint** (section 5): v2, starting with the
  children-terminal change move.
- **chi-k under truncation** (section 6): v1 fixes k; reconciling the hyperprior
  with the truncated leaf law is deferred.
- **dbarts.h exposure:** none in v1 (the ordinal/robust-errors precedent). A
  future flat-C entry could expose the direction vector for embedded use once
  there is demand.

## 11. Costs, risks, and confidence

- Per-constrained-predictor cost. mBART reports ~5x per-iteration slowdown for
  each constrained predictor at a 20x20 grid; our honest-but-lower-dimensional
  integrals (mostly 1-D closed form, 2-D only for a constrained-axis birth) should
  beat that, but constrained forests are inherently slower than unconstrained -
  acceptable because the cost falls only on declared columns and unconstrained
  fits are untouched (section 8).
- Mixing. Birth/death-only (section 5) plus single-site leaf Gibbs can mix slowly
  when many leaves are mutually constrained (tight truncation intervals). The
  recovery test watches for it; random-scan Gibbs and the children-terminal change
  move are the tuning knobs held in reserve.
- Neighbor geometry is the correctness-critical new code; a wrong adjacency test
  silently fits the wrong constrained model. The brute-force oracle test (section
  9) is the guard.
- Confidence in the constrained-draw algorithm AS SPECIFIED: HIGH for
  birth/death plus single-site truncated-normal Gibbs. The paper is explicit and
  the two structural claims that make it clean are verified here from the text -
  the acceptance conditions on mu_same and reduces to the engine's existing ratio
  shape (eq. 4.11), and a birth can never empty the cone (section 5). The
  deviation from mBART (honest normalizers instead of d=1 plus a retuned tree
  prior) STRENGTHENS correctness - it is what lets the exact-posterior gate pass -
  at a modest, columns-only compute cost. The residual risk is engineering, not
  algorithmic: the neighbor-adjacency geometry and the tree-granularity leaf-draw
  seam, both covered by the component tests above.

## Plan-vs-code note

The plan stub (docs/plans/monotone-bart.md:22-25) frames fork 3 as "exact
constrained marginals vs mBART's approach," implying mBART's approach is not
exact. Finding: mBART's TARGET is exact (the conditional-on-mu_same marginal,
eq. 4.11, hits the true constrained posterior); only its IMPLEMENTATION
approximates (grid + d=1 + retuned tree prior, Section 4.3, quoted verbatim in
section 4 B'). The real fork is
therefore B (mBART's approximated implementation) vs B' (mBART's exact target with
honest, lower-dimensional numerics), not exact-vs-mBART - and B' is both exact and
what the plan's own exact-posterior gate requires. The plan's "sequential
conditional truncated normals" description (line 23) is correct for the leaf draw
and for the constrained-axis birth redraw; it is NOT how the structure-move
MARGINAL is scored (that is a truncated integral, not a draw), a distinction
section 4 makes explicit.
