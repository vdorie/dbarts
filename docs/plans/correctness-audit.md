# correctness-audit

agent: opus derivers, paired and independent; fable adjudicates
rng: neutral (findings only; any fix becomes its own gated item)
budget: findings per block recorded here; no code changes

Review 1 of the retrospective program (retrospective-reviews.md).

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
   latent): drawBirthableNode's single-node branch skips the
   availability check, so a root-only tree with NO available variable
   forces a birth and drawRuleForVariable indexes
   data.types[(size_t)-1] - out-of-bounds. The all-constant-column
   case does not trigger it empirically (a constant column still
   quantizes to >= 1 cut; births propose, hit the empty-leaf veto,
   reject - verified by running bart on constant x), but zero-cut
   columns are reachable at least via setCutPoints on a root-only
   sampler (invalidCutPoints only protects existing splits), and the
   death branch is equally unguarded for single-node trees. Filed as
   fix item moves-degenerate-root-guard.
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

MAJOR FINDING (both derivers independently, opposite orientations;
orchestrator re-derived and concurs): the CHANGE move's acceptance
(moves.hpp changeMove, the exp(yLogPi + yLogL - xLogPi - xLogL) with
no transition term) omits the proposal-density ratio. The proposal
draws the new variable from the prior but the new RULE uniformly over
the descendant-valid good set (ordinal: findGoodOrdinalRules;
categorical: reject-until-valid), while the acceptance retains the
node's local rule prior normalized over the ancestor-only interval.
These cancel only for same-variable redraws. Cross-variable changes
are mis-weighted by [p_var(v')/p_var(v)] * [|Valid(v)|/|Valid(v')|]:
the chain's stationary distribution carries an effectively SQUARED
rule prior at changed nodes, biased toward low-cardinality /
descendant-constrained variables. Invisible when all variables have
equal cut counts and trees are shallow - which is exactly the
existing exact-posterior gates' regime. INHERITED: the deleted
classic engine's changeRule.cpp computes the identical
pure-pi-ratio acceptance (verified in git history at b354f3a~1), so
this is a CGM-lineage defect dbarts has carried since its origin, not
a rewrite regression. Deriver B's toy check shows detailed balance
failing by exactly 2x in a 4-vs-2-cut config; deriver A's analysis
predicts stationary rule-prior mass ~ 1/a_v^2 rather than 1/a_v.
Fragility notes (no action): swap's correctness rests on the unstated
lemma that a child can never carry its parent's exact rule; the
64-attempt categorical abort is variable-dependent and compounds the
asymmetry.

Verification in flight: engine-level exact-enumeration test
(single tree, two ordinal predictors with unequal cut counts;
engine long-run root-split frequencies vs the exact posterior AND vs
the predicted wrong target) - the protocol's numerical step before
the finding becomes a fix decision for VD. Fix would be
posterior-changing (strictest gates) and changes 15-year-old
sampler semantics; VD decides with the evidence.
