# check-shape-fixes: the count mixing gate's shape, a calibration arm for k, a weekly full run

Status: PLANNED

agent: sonnet implementer, one (R scripts and CI workflows only); no reviewer needed beyond the coordinator's read.
rng: NEUTRAL. No package code changes (nothing under R/, src/, inst/ or man/).
window: pre-release; the by-hand calibration run is a step BEFORE the merge to main.
budget: ~400 lines.

## Goal

The count model's mixing gate measures agreement on a shape whose posterior spreads, starts its chains apart, and
has no rule for constant chains. The calibration suite gains one arm that draws k under its default
hyperprior. The monotone enumeration gate's full mode, and the other exact gates' full modes that fit, run weekly on
a longer limit. The calibration suite can be started by hand on this branch.

## Context

- dec-B301 (the maintainer, 2026-10-07: "Fix the shape.") ruled the first two items and the by-hand run; dec-A140
  recorded the three choices it reshapes; dec-A132 (the maintainer: "If we have weekly tests or somesuch we can do
  exhaustive on an infrequent basis.") ruled the weekly run. Root TODO: check-shape-fixes and
  weekly-full-exact-gates, one slice.
- [negbin-mixing.R](../../benchmarks/R/negbin-mixing.R) is the mixing gate; `splitRhat` held the constant-chain
  rule, `runCell` the cases.
- [sbc.R](../../benchmarks/R/sbc.R) is the calibration harness; `sbcReplication` is the plain arms' replication,
  `sbcMatrixConfigs` and `sbcMatrixFunctionals` the matrix's admission level. [sbc.yaml](../../.github/workflows/sbc.yaml)
  runs the matrix weekly. No arm there ranked k: nbinom holds it at 3 / 0.68 and the plain arms at 2.
- [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) runs every exact gate in quick mode per push and the
  monotone enumeration gate in its own job; [equivalence.yaml](../../.github/workflows/equivalence.yaml) and
  [rchk.yaml](../../.github/workflows/rchk.yaml) are the weekly (Monday) runs.

Ran, on the tip with the installed package, one R process at a time, two cores at most:

- the mixing gate as it stood (quick: both cases' chains started at one value; the r0 = 2 case held r = 2 in every
  draw of both chains, the r0 = 5, n = 2000 case 5 or 6 in 99%);
- exploratory fits of the shape's posterior over six candidate cells, and the sampler's interface for the start: no
  setter, but the stored state carries a shape per chain and `setState` installs it (read back by `getShape`);
- the k step's interface: `setLeafPrior(normal(k = ))` installs a fixed k, and writing the hyperprior again keeps
  that k as the chain's current value (probed, and made a self-check of the arm);
- k's autocorrelation under the arm's own prior draws (ACF 0.1 at lags 113 to 619 over five datasets);
- every exact gate's full mode, timed (below); the enumeration gate in quick mode for one prior.

Only read: the engine's k update (`ChiKHyperprior`: k^2 | leaves is gamma with shape (M + nu) / 2, which makes the
prior k = scale * sqrt(chisq(nu)), the median 1.906 for chi(1.5, 2) the help states), and GitHub's rule that
`schedule` and `workflow_dispatch` bind to the default branch.

## The rules

1. Mixing gate (dec-B301). Two cells whose posterior on r spreads: r0 = 8 at n = 400 and r0 = 10 at n = 500, in
   place of r0 = 5 at n = 2000 and r0 = 2 at n = 500 (the latter a single value; the former 5 or 6). The rule that
   chains constant at one shared value agree is gone: split-Rhat with no spread in any half is undefined and fails.
   Chain 1 is set to the shape 1 and chain 2 to 50 (the grid's ends) through the state before the run. Each chain
   keeps 2000 draws, not 500: at 500 a correct sampler read split-Rhat 1.059 on one seed, the chains' block means
   drifting together by 4000 draws. Check (i) is now that a chain leaves its own start value.
2. Calibration arm (dec-B301). `probit-k`: the plain probit arm with k drawn under chi(1.5, 2), the default for
   binary fits. theta0's k is drawn by the harness, installed as a fixed k, the leaves are drawn at it, and the
   fit starts from a second independent draw handed back to the hyperprior; k is ranked beside the existing
   functionals. It joins the matrix (M = 77 + 7 = 84). Settings follow the measurement: R = 100, L = 100, thin 1000
   and a 30000-sweep burn (sbc.R supplies it), since k and its leaves are a funnel.
3. Weekly (dec-A132). A new workflow, `exact-gates-weekly.yaml`, Mondays, runs the full mode of every gate in
   exact-gates.yaml's main list (read from that file) except bcf-latent-exact, and the monotone enumeration gate
   and the successive-conditional check in full mode, one job per prior, on limits sized from the timings.
4. By hand, before the merge. sbc.yaml already carries `workflow_dispatch` and a push trigger on its own file;
   this slice edits the file, so the push that carries it starts the suite on bartcore. The command and the
   reading are in the report.

## Steps

1. negbin-mixing.R: cells, `startShapes` set through `setState` (with a read-back refusal), the run driven on the
   sampler (`samplerOnly`), the posterior predictive drawn from the run's test log means and shapes, `splitRhat`
   returning NA for no spread and the cell failing on it, header and README text.
2. sbc.R: `sbcConfigProbitK`, `sbcKDraw`, the k install and hand-back in `sbcReplication`, the k rank, the
   `sbcCheckKHandBack` self-check, the matrix entry and count, the burn default. sbc.yaml: the matrix row, the
   count and the timeout note.
3. exact-gates-weekly.yaml (new); exact-gates.yaml header and one job comment say where full mode runs.

## Proof

- Mixing gate. Tip, quick: passes (13 s); full (seeds 1 to 3): passes (40 s); six seeds: the r0 = 10 cell passes
  all six, the r0 = 8 cell misses its 95% set on seed 4 (a 95% set missing r0 once in twelve, not a mixing finding;
  full mode runs seeds 1 to 3). A mutated build whose shape draw is discarded (a scratch copy, private library)
  fails every check but the band: both chains stay at their starts, split-Rhat undefined, coverage 0.98.
- k arm. Tip at the shipped settings: all seven functionals pass, k's chi-square p 0.90 and ecdf difference 0.074
  of a 0.195 band; at thin 300 and a 10000 burn k's ranks piled at the top (31 of 200 in the last bin) and
  shrank away at thin 1000, so the setting is the evidence's, not a guess. A mutated build with every k draw
  multiplied by 1.4 flags all seven functionals, k's ecdf difference 0.99.
- Weekly. Each gate's full time (laptop, one core, seconds): bd-balance 3, swap-balance 13, perturb-balance 21,
  rule-gibbs-balance 29, change-balance 113, aft-exact 1, aft-hetero-pit 0, backfit-exact 11, bcf-exact 30,
  bcf-exact-weak 1, bcf-exact-restricted 1, categorical-exact 4, heteroscedastic-exact 3, linear-exact 2,
  hazard-exact 11, hurdle-exact 9, mask-redraw-exact 402, multinomial-exact 156, negbin-exact 16, negbin-mixing 40,
  ordinal-exact 12, t-exact 2, monotone-reference 9, hazard-reduction and hurdle-reduction 0, bd-balance zeroweight 2,
  monotone successive-conditional 25 for both priors; all passed. logistic-reference did not finish in 580 s: its
  one-seed ensemble comparison took 32 s and runs 20, so about 13 minutes. The enumeration gate in quick mode took
  272 s for the leaf prior, so about 14 minutes full per prior. The list sums to about 28 minutes; the hosted
  runner took about five times the laptop's time on the enumeration gate (23 minutes against 4.5).

## Gates

lintr with zero lints, air format check, doc-freshness, `verify-anchors`; each changed script in the mode CI runs
it in; the three workflows parsed.

## Out of scope

The package; the carried reference fits; a count-arm SBC that draws k (dec-B301 weighed it: run time for what
probit shows); bcf-latent-exact in the weekly run (about 40 minutes on the laptop alone, not timed here; a matrix
row once someone has); a setter for the shape.

## Calls made in planning

- The two mixing cells: the maintainer asked for one case replaced; the r0 = 5, n = 2000 posterior is 5 or 6 in
  99% and goes with it, since a run with no 6 would be undefined.
- 2000 draws a chain in the mixing gate, and the failure of an undefined Rhat, over keeping 500 and widening.
- Chains set to the grid's ends through the state, since no setter exists; the first sweep redraws r, so the start
  tests that r moves, not that it stays far.
- thin 1000 and a 30000 burn for probit-k, by the chain-length evidence; R = 100 as the gaussian arm.
- A separate weekly workflow reading the gate list from exact-gates.yaml, over a schedule inside exact-gates.yaml
  (its concurrency group would let a push to main cancel the weekly run).
- bcf-latent-exact left out of the weekly run.
