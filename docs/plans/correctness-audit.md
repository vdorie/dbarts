# correctness-audit

agent: opus derivers, paired and independent; fable adjudicates
rng: neutral (findings only; any fix becomes its own gated item)
budget: findings per block recorded here; no code changes

Review 1 of the retrospective program
(docs/plans/archive/retrospective-reviews.md).

## Disposition

Every finding this audit filed as a fix item was confirmed and has
since been fixed. The blocks below are the derivations and the
evidence, not a list of live defects; each finding's own paragraph
names its fix.

- Block 1, the degenerate-root guard gap: birthOrDeathMove returns a
  no-op for a root-only tree with no birthable node
  (src/bartcore/moves.hpp).
- Block 2, the CHANGE move's missing proposal density:
  logProposalCorrection is computed and consumed in changeMove
  (src/bartcore/moves.hpp).
- Block 3, the chi hyperprior's mislabeled degrees of freedom:
  ChiKHyperprior::draw uses shape 0.5 (M + nu)
  (src/bartcore/model.hpp), with the binary default relabeled
  chi(1.5, Inf) so default draws are unchanged.
- Block 4, the sigma degrees of freedom over zero-weight rows:
  numPositiveWeights_ feeds the sigma draw (src/bartcore/model.hpp).
- Block 5, the empty-leaf chi-k count: constant and linear leaves now
  match GP, as that finding's own paragraph records.
- Blocks 6 and 7, test fits under a multi-forest coupling:
  refuseUndefinedTestFits guards bartcore_setTestPredictor,
  bartcore_setTestOffset and bartcore_predict
  (src/R_interface_bartcore.cpp).

Each fix carries its own archived plan under docs/plans/archive/.

## Goal

Every acceptance ratio and conjugate update in the engine is either
independently re-derived and CONFIRMED, or flagged with the derivation
that contradicts it. The exact-posterior gates verify a few enumerable
configurations end-to-end; this audit checks the formulas themselves,
term by term, where a compensating pair of errors or an untested
configuration could hide.

## Method

Per block, two independent Opus derivers with different orientations:
deriver A derives the correct expression from the model statement
first, then diffs against the code; deriver B translates the code into
math and checks self-consistency (detailed balance, normalization,
domain edges). Neither sees the other's work; neither trusts comments.
Verdicts per term: CONFIRMED (matching derivation) or DISCREPANCY
(derivation + the exact code location). Fable adjudicates; a surviving
DISCREPANCY gets numerical verification (quadrature or simulation at
fixed inputs) before becoming a fix item. Two derivers at a time,
blocks sequential (the program's fan-out constraint).

## Blocks, in priority order

1. Birth/death: acceptance ratio incl. p_birth/p_death step
   probabilities at boundary trees, birthable/prunable node selection
   (drawBirthableNode, probabilityOfSelectingNodeForBirth,
   birthableNodeExists, probabilityOfBirthStep), rule-draw vs
   tree-prior cancellation on the cut grid, DART-weighted variable
   selection in both directions, growthProbability depth terms, the
   empty-leaf veto's reversibility interaction (moves.hpp; tree prior
   in model.hpp).
2. Change/swap: good-rule rejection samplers and their forward/reverse
   cancellation claims, categorical gauge preservation, swap validity
   walk, likelihood-only acceptance (moves.hpp, tree.hpp).
3. Hyperpriors and DART: ChiKHyperprior k^2 gamma posterior (the
   0.5/scale^2 rate term), DART normalized-gamma Dirichlet update and
   the alpha grid draw, tau priors and the grouped tau posterior
   (model.hpp; grouped pieces in chain.hpp).
4. Response models: gaussian sigma draw and its qchisq-calibrated
   prior, fixed-sigma, probit truncated-normal latents, logistic
   Polya-Gamma draw + working weights/response construction, weights
   placement in every family (model.hpp, chain.hpp).
5. Leaf marginals and draws: constant (integrated likelihood +
   posterior draw + k/sigma placement), linear (U'WU marginal,
   posterior slope draw, calibration inheritance), GP (nugget,
   over-cap fallback delegation) (model.hpp).
6. BCF and grouped: calibration map, a/b glue draws, two-forest
   residual accounting; grouped-intercept Gibbs and its offset
   plumbing (chain.hpp).
7. Cross-cutting algebra: residual roll and totalFits maintenance,
   scale transforms (range scaling, original-scale sigma in state),
   warm-start/state-restore rebuild math (chain.hpp, sampler.hpp).

Blocks 4-6 have partial end-to-end cover (logistic-reference,
categorical-exact, bcf-exact); 1-3 and 7 are the least directly gated
and run first.

## Constraints

- Findings only; the tree does not change under this item.
- Derivers read the code and the model statement; they do not read
  docs/design 'derivation' prose as authority (it may share a wrong
  source with the code).
- Every DISCREPANCY that survives adjudication gets a numerical check
  before it is called real.

## Verification

- Per-block verdict tables recorded in this file's Status; real
  findings become fix items with the exact-posterior gates as
  arbiter.

## Status

Block 1, birth/death (2026-07-08): CONFIRMED by both derivers
independently - every factor of the implemented acceptance ratio is
exactly pi(T')q(T|T') / pi(T)q(T'|T) for the proposal the code draws.
Jointly confirmed: likelihood ratio over the changed branch; growth
and no-split tree-prior factors at the correct depths; the rule
prior/proposal cancellation (ordinal interval, categorical gauge
pattern, missing-direction coin - proposal IS the live prior, both
omitted); node-selection counts recomputed on the post-move tree in
both directions; p_birth/p_death boundary handling (root forces
birth, reverse density uses P_birth(T) = 1 not the constant);
move-class probability cancels; DART weights enter proposal and prior
identically over the same available set and cancel; the empty-leaf
veto defines one consistent surrogate target across all moves.

Findings routed onward:
1. Degenerate-root guard gap (deriver A #17, adjudicated REAL but
   latent): drawBirthableNode's single-node branch skipped the
   availability check, so a root-only tree with NO available variable
   forced a birth and drawRuleForVariable indexed
   data.types[(size_t)-1] - out-of-bounds. The all-constant-column
   case did not trigger it empirically (a constant column still
   quantizes to >= 1 cut; births propose, hit the empty-leaf veto,
   reject - verified by running bart on constant x), but zero-cut
   columns were reachable at least via setCutPoints on a root-only
   sampler (invalidCutPoints only protects existing splits), and the
   death branch was equally unguarded for single-node trees. Filed as
   fix item moves-degenerate-root-guard, and FIXED there:
   birthOrDeathMove now returns a no-op when a root-only tree has no
   birthable node.
2. Cross-move categorical prior flag (deriver B): for > 54 reachable
   pooled categories, ruleForVariableLogProbability uses an
   approximate closed form while the draw density is exactly
   1/(2^R - 2); birth/death never evaluates the prior (cancellation)
   but change/swap may - routed to block 2 as a targeted question.
3. Deriver A fragility note (no action): the entire rule prior incl.
   DART rests on the exact proposal-equals-prior identity with no
   backstop term; any future proposal that deviates from the live
   prior silently breaks birth/death. Worth a comment or an assert
   when that surface is next touched.

Block 2, change/swap (2026-07-08): SWAP CONFIRMED by both derivers
(deterministic involution, symmetric proposal - the swappable set and
child coin are topology-only - so the subtree pi-ratio is exactly the
MH acceptance; validity walk and mask-pool restore clean). Veto
composition and rejected-move state restoration CONFIRMED. Block-1's
targeted >54-category question RESOLVED as a false alarm by both: the
"approximate" closed form is an exact algebraic rewrite of
-log(2^R - 2) (error O(2^(1-R)), sub-epsilon), so no cross-move
target divergence exists there.

MAJOR FINDING, since FIXED (both derivers independently, opposite
orientations; orchestrator re-derived and concurs): the CHANGE move's
acceptance (moves.hpp changeMove, the exp(yLogPi + yLogL - xLogPi -
xLogL) with no transition term) omitted the proposal-density ratio.
The proposal drew the new variable from the prior but the new RULE
uniformly over the descendant-valid good set (ordinal:
findGoodOrdinalRules; categorical: reject-until-valid), while the
acceptance retained the node's local rule prior normalized over the
ancestor-only interval. Those cancel only for same-variable redraws.
Cross-variable changes were mis-weighted by
[p_var(v')/p_var(v)] * [|Valid(v)|/|Valid(v')|]: the chain's
stationary distribution carried an effectively SQUARED rule prior at
changed nodes, biased toward low-cardinality /
descendant-constrained variables. Invisible when all variables have
equal cut counts and trees are shallow - which is exactly the
existing exact-posterior gates' regime. INHERITED: the deleted
classic engine's changeRule.cpp computes the identical
pure-pi-ratio acceptance (verified in git history at b354f3a~1), so
this was a CGM-lineage defect dbarts had carried since its origin,
not a rewrite regression. changeMove now carries the missing term as
logProposalCorrection. Deriver B's toy check shows detailed balance
failing by exactly 2x in a 4-vs-2-cut config; deriver A's analysis
predicts stationary rule-prior mass ~ 1/a_v^2 rather than 1/a_v.
Fragility notes (no action): swap's correctness rests on the unstated
lemma that a child can never carry its parent's exact rule; the
64-attempt categorical abort is variable-dependent and compounds the
asymmetry.

Verification (2026-07-08): CONFIRMED empirically. Engine-level
exact-enumeration test, benchmarks/R/change-balance.R: n = 100,
single tree, x1 with 19 cuts vs x2 with 2 cuts, fixed sigma, 1M kept
draws vs a depth-6 memoized region DP for the exact posterior and for
the predicted wrong target (rule prior squared); truncation moves
root marginals < 3e-7. Result: P(root = x2 | split) engine 0.2988 vs
exact 0.0774 vs wrong-target 0.3704 - the engine fails the exact arm
at z = +479 and sits 76% of the way to the wrong target (between the
arms as predicted, since birth/death still target the correct pi).
Within-variable cut distributions match the exact posterior (max gap
0.02), so the corruption is purely between-variable - exactly the
predicted channel. Control with equal root cut counts: root margins
match (z = -2.2) validating the enumeration, with a small honest
residual (z = -27, 19x smaller) from depth >= 1 changes whose
descendant-constrained intervals still differ - the same mechanism.

Block 3, hyperpriors and DART (2026-07-08): CONFIRMED throughout by
both derivers except one adjudicated labeling defect. Jointly
confirmed: the k draw's rate (0.5 * sum(param^2)/leafScale^2 + the
finite-scale 0.5/scale^2 term; infinite-scale limit clean; the leaf
sum accumulated raw on the internal scale, no k^2-vs-k^4 error, reset
per sweep); the full DART machinery (Dirichlet alpha/p + m_j with
counts recounted over all trees each sweep; normalized-gamma with a
1e-300 floor and uniform fallback; the alpha grid living in lambda
space with the Beta prior expressed in lambda so NO Jacobian is owed,
normalizer constants exact; alpha conditions on s not counts - the
correct Linero step; updateDelay holds s uniform then consumes fresh
counts); tau priors (cauchy and gamma parameterizations match R
exactly, scale-family internal conversion owes no Jacobian), the tau
posterior conditioning only on the group effects (correct block), the
weighted per-group conjugate intercept update (reduces to R's
unweighted form at unit weights), offset plumbing and recording
de-scaling coherent end to end.

ADJUDICATED FINDING, since FIXED (both derivers, same math, framing
reconciled): the k posterior shape was 0.5(M + 2 nu - 1), the exact
Gibbs step for a chi(2 nu - 1) prior on k - NOT the chi(nu) the
degreesOfFreedom argument and docs described. chi(1.25) delivered
chi(1.5); the two coincide only at nu = 1. Bit-identical to classic
parameterPrior.cpp (a Jacobian slip in the original k^2 derivation,
inherited). The sampler was internally valid for its own prior; the
defect was the LABEL. Filed as chi-hyperprior-df, where VD picked the
code-fix: ChiKHyperprior::draw now uses shape 0.5 (M + nu), and the
binary default is relabeled chi(1.5, Inf), which leaves default
draws bit-identical.

Targeted-test note for the SBC/blindspot reviews: the 1e-300
normalized-gamma floor feeds log(1e-300) into the alpha grid's
sum-log-s under extreme sparsity (alpha/p << 1, many zero-count
variables), plausibly biasing alpha low; worth a high-sparsity
calibration check.

CONSEQUENCE, as it stood before the fix: dbarts' tree-structure
posterior over-weighted low-cardinality and descendant-constrained
split variables in every configuration with unequal effective cut
counts (mixed continuous/categorical predictors especially), and had
done so since the package's origin. The fix decision was VD's:
posterior-changing class, changing 15-year-old semantics. The
candidate shapes were (a) propose change rules from the unrestricted
prior and let invalid descendants reject via pi = 0 (restores the
exact cancellation with no valid-set counting - counting is
infeasible for wide categorical masks - and removes the 64-try
asymmetry), or (b) restricted proposals with explicit |Valid| ratios
where countable. (Resolved 2026-07-08: VD picked the hybrid after the
change-move-fix stage-2 measurements; see that plan. The landed term
is logProposalCorrection in changeMove, src/bartcore/moves.hpp, and
change-balance.R is its regression gate.)

Block 4, response models (2026-07-08): CONFIRMED on every mainline
path by both derivers except one adjudicated zero-weight defect.
Jointly confirmed: the gaussian sigma^2 Gibbs step is the exact
conjugate update for the documented scaled-inverse-chi-squared prior
- lambda = sigest^2 qchisq(1-quant, df)/df computed at the bridge,
posterior draw (df lambda + S)/chisq(df + n) with the weighted
S = sum w_i r_i^2, internal range scaling coherent in and out, n = 0
reproducing the prior; fixed-sigma reaches every consumer (moves,
leaf draws, latents, BCF glue, recording) through the single chain
sigma_ with the draw gated off; probit latents are exact Albert-Chib
(sign-folded lower-truncation gives the correct upper/lower
truncated N(fits + offset, 1); working response strips the offset;
sigma pinned to 1); logistic is the exact binomial Polya-Gamma
update (omega ~ PG(w, psi) by integer replication, working response
kappa/omega - offset with kappa = w(y - 0.5); the leaf update's
omega * working recovers kappa exactly, so user weights compose
through PG without double-counting); weights appear in the correct
place in every family (constant/linear/GP suffstats, sigma SSR,
grouped-intercept precision; probit correctly weightless with the
R-level guard).

ADJUDICATED FINDING, since FIXED (deriver B; deriver A had graded the
df term only as a fragility note; orchestrator verified numerically):
the sigma posterior's degrees of freedom counted ALL n observations
(model.hpp drawSigmaSqFromPosterior, df + numObservations) while
zero-weight rows contribute nothing to S, nothing to any leaf
conditional, and are documented as "ignored" (the weights validator
warns exactly that). The draw was not the conjugate update of any
coherent model once w_i = 0 rows existed: df was over-counted by
their number, deflating sigma. Verification: paired fits, 50 real
rows vs the same 50 plus 50 EXACT DUPLICATE rows at w = 0 (duplicates
cannot alter cut points, the y scaling, the leaf prior, or
leaf-emptiness patterns, so the df term is the only open channel):
mean sigma 0.365 -> 0.048, z = -264 against the promised equality.
The first-order df ratio predicts 0.72; the observed collapse is the
stationary feedback loop (deflated sigma lets trees absorb residual
as structure, shrinking S, deflating sigma further). Filed as fix
item sigma-df-zero-weights and FIXED there: numPositiveWeights_ is
the count the posterior df now uses. Posterior-changing only for fits
with zero weights (the zeroweights equivalence scenario shifted; no
default touched).
Gate-blindspot note: the zeroweights equivalence scenario compares
the code to its own baseline and no exact-posterior gate covers
zero weights, which is how this survived.

Fragility notes (no action): the base GaussianResponse::drawSigma
always draws - only the chain's sigmaIsFixed_ gate keeps fixed sigma
fixed; logistic omega uses lround(w) while kappa uses raw w, an
integer-weights invariant enforced only by the R validator (direct
dbarts.h consumers can pass non-integer w and get a silent
omega/kappa mismatch); the probit truncated-normal sampler's NaN
fallback substitutes a sign-correct DBL_EPSILON latent on extreme
tail failure; logistic setOffset reshifts the working response
against omega drawn at the old psi (repaired by the next sweep's
refreshLatents - the standard embedded-Gibbs pattern).

Block 5, leaf marginals and draws (2026-07-08): CONFIRMED throughout
by both derivers; no surviving DISCREPANCY. Jointly confirmed, term
by term including constants and determinants: the constant leaf's
integrated likelihood and conjugate mu draw (shared tau0 = (k/scale)^2,
the retained prior normalizer that birth/death needs, empty/single-
observation/all-zero-weight edges all coherent); the linear leaf's
U'WU marginal (ridge = tau0 sigma^2, Cholesky determinant and
half-solve quadratic form, exact q = 0 reduction to the constant
form on the same dropped-constant convention) and its
N(M^-1 b, sigma^2 M^-1) draw, with calibration inherited from the
constant leaf's k for every coefficient; the GP leaf's nugget
marginal (observation noise sigma^2/w_i on V's diagonal, distinct
from the 1e-6 conditioning jitter inside C), its Matheron-rule draw
(covariance verified equal to the exact posterior), and the correct
chi-k sufficient statistic f'C^-1 f. The 60a116c U'WU cache
preserves the right invariant: the crossproduct depends only on
(ordered member list, covariates, weights); the memcmp member-list
key auto-misses on any structural change or rollback, covariate
mutations clear via regather on every setPredictor/setData path,
weight mutations via invalidateStatistics at setWeights/PG-refresh/
latent setResponse, and the residual-dependent pieces always rescan
- no path serves a stale crossproduct. The over-cap fallback and
zero-weight routing branch on identical predicates in score and
draw, and all three leaf models drop the same additive observation
constant, so it cancels in every MH ratio including size-switched
births across the GP cap - score and draw target one coherent model
per node everywhere.

Findings routed onward:
1. Empty-leaf chi-k count inconsistency (deriver A, latent):
   forced-to-zero empty leaves contribute to the k hyperprior's leaf
   count inconsistently across leaf models - constant adds 1, linear
   adds q+1, GP adds 0 - while all three add nothing to the sum of
   squares. Counting a forced zero as a genuine N(0, sigma_mu^2)
   draw inflates the shape without a matching rate term, deflating
   k. GP's not-counting is the coherent choice given forced-zero
   params. Latent (the empty-leaf veto keeps empty leaves out of
   normal fits; reachable only via mid-sampler data/weight/state
   mutations that strand one). Filed as fix item
   chi-k-empty-leaf-count: make constant/linear match GP;
   draw-neutral whenever no empty leaf exists. (Landing refinement:
   the mutation paths all collapse empty leaves - applyNewData and
   forceRefreshTrees collapse, revalidateTrees rejects - so no
   public path strands one; the defect was pure-latent, reachable
   only by fabrication. Fixed and landed with a fabrication-based
   tests/cpp regression; see chi-k-empty-leaf-count.md.)

Fragility notes (no action): the constant leaf's centered
sum-of-squares is a catastrophic-cancellation form under a large
shared residual offset (score-side rounding only; flagged for the
numerical-robustness review); the cache invariant (U'WU/kernel
depend only on members/covariates/weights and every mutation path
clears) is enforced by convention with no backstop assert - a
future response family with per-sweep weights that misreports
workingWeightsVaryPerSweep, or a covariate path skipping regather,
would silently serve stale statistics; zero-weight rows count
toward the GP max.leaf.size cap (consistent between score and draw,
a calibration surprise only); the 1e-6 jitter makes the realized GP
prior s^2(C + 1e-6 I), identically in score and draw.

Block 6, BCF and grouped (2026-07-08): CONFIRMED throughout by both
derivers; no surviving DISCREPANCY. Jointly confirmed: the
calibration map (s = n-1 sample sd of the range-scaled net-offset
working response; mu forest total N(0, s^2) and tau forest total
N(0, (sdModerate s / 0.674)^2) via the engine's nodeScale/sqrt(T)
leaf convention with k fixed at 1); the glue draws (a conjugate
against residual y* - b_z tau with precision 1/aVariance +
sum w mu^2/sigma^2; b0/b1 per-arm conjugates against y* - a mu with
an empty arm falling back to its prior; the half-Cauchy auxiliary
aVariance ~ IG(1, (a^2 + aPriorScale^2)/2), Monte-Carlo verified to
marginalize a to Cauchy(0, aPriorScale)); the two-forest residual
accounting (per-forest response resid/m and weights w m^2 with
per-observation m = a or b_z, whose crossproduct sufficient
statistics are exactly the untransformed problem's - deriver B
verified numerically that the 1e-9 multiplier floor reduces a
zero-coefficient arm's contribution to the informative-arm-only
posterior to machine precision); the sigma feed (combined
a mu + b_z tau fit with ORIGINAL user weights and the
positive-weight df); updateA/updateB gating; the grouped-intercept
Gibbs (precision sum-of-weights/sigma^2 + 1/tauG^2 - a weight sum,
not a count: the block-4 sigma-df bug class is absent here);
grouped offset plumbing (forests backfit z - b, shiftFits feeds
F + b + offset to the base latent refresh, so probit/logistic
latents and PG weights compose correctly); zero-weight rows inert
in every posterior block; recording de-scales sigma/tau/group
effects by sigmaScale and blends trainingFits with the live glue;
state save/restore round-trips a, aVariance, b0, b1, groupEffects,
groupTau on the internal scale with scratch rebuilt per sweep.
Sweep order (mu -> tau -> latents -> sigma -> glue -> intercepts)
is a coherent blocked Gibbs scan on one joint. BCF and grouping are
mutually exclusive, but by an explicit bridge refusal rather than by
a hardwired Gaussian response: the chain selects its response off the
spec's family, and gaussian, probit and logistic all build. That
interaction still cannot arise, but this block's derivations were run
against the Gaussian arm and have not been re-run for the probit or
logistic multi-forest chain.

Findings routed onward:
1. BCF testFits single-forest (deriver B, adjudicated REAL but
   latent; verified by orchestrator in chain.hpp):
   results.trainingFits had an explicit BCF blend but
   results.testFits unconditionally reported scale *
   forests_[0].totalTestFits + shift (+ testOffset) - the bare mu
   forest with the combined-response centering, no a, no b_z tau.
   A correct test-row blend is ill-defined without a test treatment
   vector (the API carries none; R-side prediction recombines
   through the facade's forestTotalFits/bcfGlue channels), so the
   channel cannot silently emit mu-only numbers dressed as fits.
   Filed as fix item bcf-testfits-guard and FIXED there: the testFits
   channel comes back NaN under BCF, refuseUndefinedTestFits guards
   the entries that would otherwise fill it, and
   results.k/variableCounts/splitProbabilities are documented as
   mu-forest diagnostics under BCF.

Fragility notes (no action): kHalfNormalMedian is the truncated
literal 0.674 (exact value 0.67449; ~0.07 percent tau-scale
inflation, and the design note itself states 0.674); the 0.674
constant is decoupled from bPriorVariance, so a host setting
bPriorVariance != 1/2 silently breaks the median-effect
calibration; BCFSpec carries no initial a/b0/b1, so updateA/B =
false can only pin the defaults 1.0/0.0/1.0; the 1e-9 multiplier
floor is a magic constant unrelated to data scale (verified benign
for the sufficient statistics); scaledResponseSd includes
zero-weight rows in the prior anchor s - the single place a
zero-weight row participates anywhere (prior calibration only);
with b0 initialized to 0 the sweep-0 tau forest sees control rows
at floored weight (warm-up transient only); recorded sigma
conditions on the previous sweep's glue while trainingFits uses the
new glue (valid Gibbs, half-step phase when eyeballing single
draws).

Block 7, cross-cutting algebra (2026-07-08): CONFIRMED throughout by
both derivers; the sole DISCREPANCY both found independently was the
BCF testFits gap already filed from block 6 (bcf-testfits-guard -
upgraded there: deriver B traced a concrete route, the then unguarded
bartcore_setTestPredictor entry accepting test predictors on a BCF
pointer, so the item's fix grew a guard on that entry.
refuseUndefinedTestFits, in src/R_interface_bartcore.cpp, is that
guard, and it covers setTestOffset and predict as well). Jointly
confirmed by full symbolic trace: the fused residual roll (per-tree
treeY is exactly the response net of all other trees' current fits,
accept/reject-agnostic because parameters are always redrawn and
rejected moves restore node sufficient statistics exactly -
birth/death/change/swap all verified); totalFits maintenance
(telescoped through the roll in a pure run, re-anchored by direct
summation on every mutation and restore path); the two-forest roll
(each forest's transformed response/weights recomputed from the
other's fresh totals every sweep); every mutation entry point leaves
membership, params, fits, totals, residuals, scaling, and latents
mutually consistent before the next sweep or read (setData recovers
params before the store moves; setPredictor validates all trees across
all chains two-phase before any rebuild, and its per-observation
session commits all-or-none; setCutPoints force-refreshes
unconditionally by design); the affine scale map and its inverse are
exact, sigma lives internally and is original-scale at every boundary
(draw, record, serialize, restore), and the prior re-anchors on every
range change; state store/restore is symmetric and validates cuts and
occupancy all-or-none before mutating (semantic, documented
non-bitwise); keepTrees flatten/unflatten indexing is uniform across
store, read, and predict including circular wrap, every read walking
the ring from the write cursor over the recorded draws (before 1.0-0
the readers walked raw slot order, which rotated a store two recorded
runs had written); the per-sweep callback observes the previous
sweep's fully-coherent end state.

Fragility notes (no action): totalFits is never exactly re-anchored
inside a pure run - the telescoped rounding error random-walks
(empirically ~1e-12 relative at the 1e5-sweep workload; treeFits
and state export stay exact; never assert bitwise totalFits ==
sum(treeFits) in the pure-run path), and rebuildFitsFromParameters
updates incrementally so it preserves entry drift where the other
mutation paths zero-and-resum; BCF correctness is silently coupled
to constant leaves (vector/function branches draw against
response workingWeights, not the forest-transformed weights - the
static_assert is the only guard); setState accepts an
observation-count-mismatched gaussian state when occupancy happens
to validate (caller contract, unguarded); range-changing mutations
leave leaf params on the old internal scale until the next sweep
re-adapts (documented approximate continuation; no read occurs
before a sweep on the normal path); installForest applies the
donor's scale anchor to the destination's response (correct for
the multi-start use case only); zero-weight rows influence tree
STRUCTURE though not the numeric posterior (the empty-leaf veto
keys on numObservations, so zero-weight-only leaves are legal
split outcomes and enter varcounts/DART probabilities - routed to
the SBC and gate-blindspot reviews); an under-filled savedTrees
buffer (numSamples < capacity) leaves default slots that predict
iterates over; the callback's sweepIndex is the raw iteration
index including thinned sub-iterations.
