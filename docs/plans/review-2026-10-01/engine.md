# Review 3 - lens: engine

Tree: 01dee4b4 (detached), library r3-lib. Diff base 7ad0bbea.

Covered (by probe or by reading changed code whole):
- Monotone leaf (model.hpp MonotoneConstantGaussianLeaf: branch score, cone quadrature, pair and
  single-leaf truncated draws, prior draw by rejection under "joint"); fitted-function monotonicity on
  grids, two constrained axes, NA in constrained and free columns, ordered-factor axis incl. an unused
  level, both priors.
- Geweke successive-conditional probes (my own, fixed sigma, 5 trees, n = 40, pivots plus forest
  functionals vs the engine's prior): plain; level-fibre step on, auto at a frozen mixture, off at a
  frozen mixture; swap, perturb, rule_gibbs mixtures; numeric + unordered factor + ordered factor with
  NA under all five moves and birth/death only; monotone "joint"/"leaf" with and without NA; monotone +
  level step; probit and logistic through setResponse. All |z| < 3 at 1000-3000 chains (one perturb
  z = 2.9 did not reproduce at two other seeds with 3000 chains).
- drawLevelShift algebra (conditioning, totalFits staleness, declines), chi-k exact draw and k = Inf
  limits, multi-forest leaf-prior writer (setForestFixedK, setForestMapSd, mapLeafScale), multinomial
  zero-trial composition, AFT censored redraw against the variance surface incl. setResponse at
  updateScale, variance-forest re-anchor, linear-leaf statistics-cache invalidation sites.
- Thread-count invariance (n = 20000, test rows, gaussian and logistic, chains 1 and 2, threads 1 vs
  2): bitwise identical. Worker exception/cancel path in Sampler::run read.
- Mutation invariants: setPredictor rollback (forceUpdate = FALSE) with columns named twice and in
  both orders - continuation bitwise identical to an untouched twin; forced setPredictor, whole-matrix
  setPredictor, setCutPoints, setData under plain/leaf/joint - training fits equal predict() to 1e-15.
- storeState + copy() continuation across 15 configurations (all families, variance forest, aft +
  variance, monotone both priors, monotone probit, chi-k, level step, all moves): same trajectory,
  differences at 1e-14 rounding only.
- Edge cases: n = 1, n = 2, constant predictors, all-but-one zero weights, weights 1e-300 / 1e300,
  y scaled 1e+-200, constant y, nbinom all zeros, ordinal with empty middle category, aft all censored,
  variance forest at n = 2 and constant y: no NaN, no crash.

Not covered: the "leaf" prior's order counter and pair ratio beyond what the Geweke and one-cut probes
exercise (gated by monotone-exact-enumeration.R and tests/cpp); nbinom dispersion grid, hazard,
hurdle, bcf/amplitude glue exactness (each has its own exact gate); GP and linear leaves; grow-from-root;
sparse predictor replacement; Windows threading.

## engine-01 - MAJOR - monotone cone probability is inaccurate below 1e-12 and zero beyond ~38 sd

Location: src/bartcore/model.hpp, MonotoneConstantGaussianLeaf::coneProbability (via
monotoneIntegrate / monotoneAdaptiveSimpson), twoLeafCoupledLogMarginal; also normalMass in
oneLeafLogMarginal.

Claim: the constrained-axis birth/death score takes log of a linear-space adaptive Simpson integral
with an ABSOLUTE tolerance of 1e-12, so whenever the touched pair's posterior cone mass is below
~1e-12 (the pair's data contradict the constraint by more than ~7 joint posterior sd) the refinement
never runs and the 16-panel estimate is used, off by up to ~19 nats; past ~38 sd the integral (and
normalMass) underflows to 0 and the move's -HUGE_VAL "infeasible" sentinel fires, so such births are
never accepted and such deaths always are. The sampler then follows the wrong law under both priors,
contradicting the design's "each targeted exactly".

Probe 1 (r3-engine-cone.R: line-for-line R transcription of coneProbability, monotoneIntegrate and
monotoneAdaptiveSimpson; root birth, no frozen neighbours, so exact = pnorm(-gap, log.p = TRUE)):

    sR/sL = 1.2:  gap 6 err 0.000, 7 +0.077, 8 +0.119, 9 -0.239, 15 -0.185, 40 -Inf (exact -804.6)
    worst over sR/sL in {0.05..20}:  gap 8 -0.59, 10 -4.28, 15 -3.42, 30 -14.13, 35 -18.97 nats

Probe 2 (r3-engine-onecut.R n0 n1 delta sigma nSamp): one constrained binary predictor, increasing,
one tree, fixed sigma; n0 rows at x = 0, n1 at x = 1 sitting delta below. Exact law of {root, split}
from closed forms (the split tree's restricted CGM prior includes the (1 - 0.95/4)^2 its vetoed
children owe, validated by the unconstrained twin r3-engine-onecut-free.R and the gap-2 control):

    gap 2.0 (control) leaf : P(split) exact 0.3519, sampler 0.3518
    gap 7.2  sU/sL 19.7 leaf : logCone exact -29.02 engine -30.94 | P(split) exact 0.1482, sampler 0.0249
    gap 7.2             joint: exact 0.0800, sampler 0.0124
    gap 9.9             leaf : logCone exact -51.78 engine -57.51 | exact 0.1126, sampler 0.0003
    gap 13.2            leaf : exact 0.0730, sampler 0.0037
    gap 19.6 sU/sL 5.0  leaf : exact 0.0301, sampler 0.0217   (100k draws; 5-se band <= 0.0025)

In every row the sampler matches the law implied by the engine's quadrature, not the exact one.

Reach: any constrained-axis birth whose two children are of unequal size and whose data run against
the constraint by more than ~7 posterior sd - common with large leaves, a misspecified or
nearly-flat constraint, or partial residuals other trees overshoot. Fitted values move little
(the pair sits near equality either way); split counts, structure posterior and mixing do not.

Why gates missed it: monotone-exact-enumeration.R designs keep gaps under ~7 sd, and its own
reference (orderProbability, a linear trapezoid over m +- 10 s) returns -Inf in the same regime (a
cX design with mu = c(1.2, 0.8, 0.4, 0), sigma 0.05 gives logPostCone = -Inf for every split tree),
so it cannot see it; its quadrature self-check stops at gap 3; SBC draws data from the
(monotone-consistent) prior; tests/cpp has no deep-tail cone case.

Fix: integrate relative to the integrand's peak in log space (the device monotoneInvertLogConcave
already uses) or make the tolerance relative to the running estimate; use logStandardNormalMass in
oneLeafLogMarginal; the no-frozen-bound case has the closed form log Phi((mU - mL) / sqrt(sL^2 + sU^2)).

## KNOWN (TODO truncated-normal-upper-tail) - new evidence of bias, not only resolution

ext_rng_simulateTruncatedNormalScale1's bulk branch (gap > 0) on an interval far above the mean, used
unreflected by OrdinalResponse::drawLatents and the bridge's ordinal latent draw. R emulation of the
same Rf_pnorm5/Rf_qnorm5 arithmetic, 1e5 draws, mean 0: (8.2, 9]: 2 distinct values, mean 8.606 vs
exact 8.318 (0.29 sd biased); (7.5, 8]: 282 values, mean 7.6188 vs 7.6192. The record says
"lose precision"; it is also biased. Reach is a latent 8 sd from its category, so MINOR.

## Checked and found correct

- Level-fibre step (drawLevelShift): conditional c = u - v (1'u)/(1'v) is the exact Gaussian
  conditioning; totalFits/roll bookkeeping consistent; Geweke passes on, auto-frozen, and monotone.
- swap, perturb, rule_gibbs (incl. categorical, ordered factor, NA routing): Geweke passes.
- Monotone fitted function is monotone on every grid probed (NA in constrained and free columns,
  ordered-factor axis incl. unused level, after setPredictor/setCutPoints/setData reseeds).
- Monotone pair/leaf truncated draws reflect upper-tail intervals (no bulk-branch defect there).
- chi-k: shape 0.5 (M + nu) correct; empty leaves excluded from M; k = Inf limits take 0 marginals.
- Multi-forest writer: setForestMapSd re-derives through the constructor's expression; the
  half-Cauchy median equals the scale; setForestFixedK refuses drawn-k and map forests.
- Multinomial zero-trial rows: PG(0) point mass skipped, zero precision and veto weight.
- AFT censored redraw reads the variance surface in its entry units at updateScale; re-anchor scales
  the surface by (prev/new)^2 and each tree factor by its m-th root.
- Thread count does not change any draw; chain RNG is per chain; worker exceptions are rethrown after
  join with the cancel flag set; SIGINT masked for workers.
- setPredictor rollback restores bitwise; storeState/copy() resumes the same trajectory everywhere.
