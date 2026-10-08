# check-shape-fixes: the count mixing gate's shape, a calibration arm for k, a weekly full run

Status: PLANNED

agent: sonnet implementer, one (R scripts and CI workflows only), over the build and two fix rounds; two independent
reviews, the second "land after fixes" with one blocking finding (the mixing gate's false-failure claim), fixed in round 2.
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
prior k = scale * sqrt(chisq(nu)), the median 1.906 for chi(1.5, 2) the help states), and, in the first
draft, a rule that `schedule` and `workflow_dispatch` bind to the default branch: true of `schedule` only
(revdep-smoke.yaml records a dispatch from bartcore with the file absent from main).

## The rules

1. Mixing gate (dec-B301). Two cells whose posterior on r spreads: r0 = 8 at n = 400 and r0 = 10 at n = 500, in
   place of r0 = 5 at n = 2000 and r0 = 2 at n = 500 (the latter a single value; the former 5 or 6). The rule that
   chains constant at one shared value agree is gone: split-Rhat with no spread in any half is undefined and fails.
   Chain 1 is set to the shape 1 and chain 2 to 50 (the grid's ends) through the state before the run. A chain
   keeps 2000 draws in quick mode and 4000 in full, not 500. Per run: a chain leaves its own start value, split-Rhat under 1.05, r0 inside the
   99.9% set, coverage of fresh counts in [0.84, 0.95]; in full mode also the 95% set holds r0 in all but at most
   2 of the 6 runs. A correct sampler failed the first draft's per-run 95% set about one dataset in twenty, so the
   checks are sized from 45 fresh datasets a cell and engine streams on the pinned datasets; full mode takes 4000
   draws because at 2000 the pinned r0 = 10 seed 3 dataset failed Rhat in 1 stream in 100 (the header of the script
   states the measured rate and what is caught and missed).
2. Calibration arm (dec-B301). `probit-k`: the plain probit arm with k drawn under chi(1.5, 2), the default for
   binary fits. theta0's k is drawn by the harness, installed as a fixed k, the leaves are drawn at it, and the
   fit starts from a second independent draw handed back to the hyperprior; k is ranked beside the existing
   functionals. It joins the matrix (M = 77 + 7 = 84). Settings follow the measurement: R = 600, L = 99, thin 1000
   and a 30000-sweep burn (sbc.R supplies it), since k and its leaves are a funnel. R = 600 flags an error of one
   half in the shape of the k^2 conditional about 80% of the time (at 100 it does not); k is flagged on its ecdf
   band or a chi-square p below the band's alpha.
3. Weekly (dec-A132). A new workflow, `exact-gates-weekly.yaml`, Mondays, runs the full mode of every gate in
   exact-gates.yaml's main list (read from that file) except bcf-latent-exact, and the monotone enumeration gate
   and the successive-conditional check in full mode, one job per prior. Each gate has its own limit, so one that
   hangs does not take the rest; the limits are sized from the hosted runner's 1 to 2 times the laptop's time.
4. By hand, before the merge. sbc.yaml carries `workflow_dispatch` and a push trigger on its own file, and a
   dispatch runs from a branch the file is on: `gh workflow run sbc.yaml --ref bartcore`. Only `schedule` waits for
   main. The command and the reading (thresholds, how many small p-values to expect) are in benchmarks/README.md.

## Steps

1. negbin-mixing.R: cells, `startShapes` set through `setState` (with a read-back refusal), the run driven on the
   sampler (`samplerOnly`), the posterior predictive drawn from the run's test log means and shapes, `splitRhat`
   returning NA for no spread and the cell failing on it, header and README text.
2. sbc.R: `sbcConfigProbitK`, `sbcKDraw`, the k install and hand-back in `sbcReplication`, the k rank, the
   `sbcCheckKHandBack` self-check, the matrix entry and count, the burn default. sbc.yaml: the matrix row, the
   count and the timeout note.
3. exact-gates-weekly.yaml (new); exact-gates.yaml header and one job comment say where full mode runs.

## Proof

- Mixing gate, on the correct sampler, full mode (4000 draws): 160 runs on fresh engine streams (100 on the r0 = 10
  seed 3 dataset, 12 on each of the other five), none failed (95% interval for a run 0 to 2.3%); seed 3's largest
  Rhat 1.027 and coverage 0.849 to 0.870, the others' largest Rhat 1.007; the 95% set missed r0 in 10 of 160 runs, all
  on seed 3; the script's own full run exited 0 on 14 engine streams (0 of 14). Estimated false-failure rate of a
  full run: about 0.2% (the 95% set count) plus 2 in 10000 (Rhat). Quick mode (2000 draws): 0 of 48 runs failed; 45
  fresh datasets a cell, 0 of 90 (coverage 0.858 to 0.929, largest Rhat 1.018). Caught, over 12 engine streams on
  the six pinned datasets: never moves, quick and full; one sweep in 100: full 12 of 12, quick 11 of 12 (2000); one
  sweep in 20: full 1 of 12 (7 of 12 at 2000 draws), quick 3 of 12 (2000); every draw one step up: full 12 of 12,
  quick 0 of 12 (2000). Not caught: every draw one step down, half the draws one step up (2000): the gate is
  one-sided. Full mode takes 73 s on the laptop (two threads), 40 s at 2000 draws.
- k arm, R = 600, L = 99, thin 1000, burn 30000: 67 minutes on the laptop (6.7 s a replication, two runs at once on
  a loaded machine), so 134 at twice that, inside the 180-minute limit. Correct sampler: all seven functionals pass,
  k's chi-square p 0.70 and ecdf difference 0.041 of a 0.082 band. Shape of the k^2 conditional + 1/2: k flags, ecdf
  difference 0.142 of 0.082 and chi-square p 0.000, 86 of 600 ranks in the lowest bin against 30. Power: the band
  at R = 600 is 0.082 and the ecdf's noise at its worst point about 0.02, so a true gap g is flagged with
  probability about pnorm((g - 0.082) / 0.02), 80% at g = 0.10; the +1/2 error's gap is 0.10 to 0.2 (mean rank 45
  against 54 on the same 100 seeds, paired ecdf gap 0.20). A tree's leaves left out of the sum of squares and every
  k draw 1.05 high flagged at R = 85 and 50. The thin-300 run, whose k ranks piled in the top bin with chi-square
  p 0.000, now flags.
- Weekly. Each gate's full time (laptop, one core, seconds): bd-balance 3, swap-balance 13, perturb-balance 21,
  rule-gibbs-balance 29, change-balance 113, aft-exact 1, aft-hetero-pit 0, backfit-exact 11, bcf-exact 30,
  bcf-exact-weak 1, bcf-exact-restricted 1, categorical-exact 4, heteroscedastic-exact 3, linear-exact 2,
  hazard-exact 11, hurdle-exact 9, mask-redraw-exact 402, multinomial-exact 156, negbin-exact 16, negbin-mixing 73,
  ordinal-exact 12, t-exact 2, monotone-reference 9, hazard-reduction and hurdle-reduction 0, bd-balance zeroweight 2,
  monotone successive-conditional 25 for both priors; all passed. logistic-reference about 11.5 minutes (the probit
  half 107 s). The enumeration gate in full mode took 13.9 minutes (leaf) and 12.6 (joint). The list sums to about
  32 minutes. The hosted runner took 1.4 to 1.8 times the laptop's time on the quick enumeration gate (6.1 and 7.9
  minutes against 4.5) and 1.0 to 1.9 times on the SBC arms.

## Gates

lintr with zero lints, air format check, doc-freshness, `verify-anchors`; each changed script in the mode CI runs
it in; the three workflows parsed.

## Out of scope

The package; the carried reference fits; a count-arm SBC that draws k (dec-B301 weighed it: run time for what
probit shows); bcf-latent-exact in the weekly run (about 40 minutes on the laptop alone, not timed here; a matrix
row once someone has); a setter for the shape.

## Calls made in planning

Made by the orchestrator, in the fix round after review:

- The mixing gate's checks are sized to fail a correct sampler under 1% of the time over engine streams and over
  datasets: a 99.9% set per run, the 95% set read over the full mode's runs, a coverage range from the measured
  spread across datasets. Stuck and one-in-100 stay caught.
- The probit-k arm runs at R = 600 inside a 180-minute limit, and k is flagged on its chi-square p as well as its
  band; no other design for the k step.
- The by-hand run, its command and its reading go in benchmarks/README.md; a dispatch is not bound to the default
  branch.
- Each weekly gate has its own limit, sized from the runner's 1 to 2 times the laptop.

Made in planning:

- The two mixing cells: the maintainer asked for one case replaced; the r0 = 5, n = 2000 posterior is 5 or 6 in
  99% and goes with it, since a run with no 6 would be undefined.
- 2000 draws a chain in quick and 4000 in full for the mixing gate, and the failure of an undefined Rhat, over keeping 500 and widening.
- Chains set to the grid's ends through the state, since no setter exists; the first sweep redraws r, so the start
  tests that r moves, not that it stays far.
- thin 1000 and a 30000 burn for probit-k, by the chain-length evidence; R = 100 as the gaussian arm (raised to 600 in the fix round).
- A separate weekly workflow reading the gate list from exact-gates.yaml, over a schedule inside exact-gates.yaml
  (its concurrency group would let a push to main cancel the weekly run).
- bcf-latent-exact left out of the weekly run.
