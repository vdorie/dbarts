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
