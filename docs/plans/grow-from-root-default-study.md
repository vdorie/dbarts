# grow-from-root-default-study

agent: opus (harness + analysis; no engine change)
rng: neutral - measurement only. A GREEN verdict AUTHORIZES a separate
  posterior-changing item; it does not perform it.
budget: one new file, benchmarks/R/grow-init.R, ~450-550 lines.
window: none - runs on any load (see Constraints).

## Goal

Decide whether bart2's init can default to n.grow.sweeps > 0, and measure
what a data-fitted init costs in LONG-RUN behaviour. Two stages: Stage 0
calibrates every threshold against measured sampling variance and gates
the study's own premise; Stage 1 is confirmatory and changes nothing.

## Context

- Item: TODO `grow-from-root-default`. Surface: `grow-from-root-warm-start.md`.
  Validity: `docs/design/grow-from-root.md` section 3(a) - posterior
  invariant, "a warm start may be any distribution".
- Corrected literature and why it settles nothing: an untracked
  companion doc; design
  decisions, overturns and evidence: `synthesis.md` beside it. This file
  restates neither and supersedes `memo.md` section 2.
- `Chain::growForestFromRoot` (chain.hpp) grows per chain, each on its own
  Mersenne Twister (sampler.hpp), so draw sequences are
  thread-count-independent. Harness idiom: `benchmarks/R/grouped-mixing.R`.
- Sizing priors below are from an unpreserved scratch probe (8 chains,
  1000 draws, 500-iterate tail, warm k=4 vs cold, matched seeds; recorded
  in `synthesis.md`). Stage 0 re-measures all of them in-harness.

## Decision

**D1. Categorical defect.** The grow.hpp arithmetic is true - categoricals
count in `numAvailable` and `availableSplitProbability`, then the
candidate loop skips them - but the harm is unestablished. Normalization
is proportional over the whole candidate vector, so relative
probabilities among split candidates are unchanged and only the
split-vs-no-split ODDS are scaled by A = prior mass on scannable
available columns, bounded in [1/p, 1] (0.5 at the root of a 5-categorical
/ 5-ordinal design), against integrated-likelihood ratios of exp(O(n)).
The floor is "grows no better than cold", not "worse than cold", and
section 3(a) documents the ordinal-only intent, so "correctness fix" is a
category error.

Fork. (i) MEASURED ARM: patch the renormalization, install a second
variant, three arms on S7 and S5. Costs an engine arc plus
re-verification of `test-bart2-grow-from-root.R` (its
`earlyRMSE(growFit) < 0.9 * earlyRMSE(coldFit)` threshold is
draw-dependent) and `test_grow.cpp`'s seed-pinned replay - to learn a
quantity whose likeliest value is zero and whose renorm direction pushes
reserved mass onto ordinal NOISE columns in exactly the designs where it
fires. (ii) DOCUMENTED BEHAVIOUR: state the A-factor and the
graceful-degradation floor in `grow-from-root.md` and the `n.grow.sweeps`
Rd; carry S7 here as a DESCRIPTIVE cell with a null prediction.
**Recommendation: (ii)**; the renormalization is prescribed nowhere in
this plan. Pre-registered reopening trigger: S7's contrast landing beyond
the COST kill margin.

**D2. Prerequisite (1) scope.** TODO lists categorical scan support as a
prerequisite before the flip is decidable (authorized 2026-08-08). D1(ii)
demotes it to a documented scope note - VD's call, not the study's;
`memo.md` settled it in a rider without saying so. Signing D1(ii) signs
the demotion; declining keeps the flip blocked behind a categorical-scan
code arc whatever this study returns. Family consistency (four refusing
bart2 branches, rbart_vi's absent surface) is a surface decision, out of
scope.

## Constraints

- Read-only on the package: no engine edit, no default change, no
  baseline re-record inside this item.
- LOAD-INSENSITIVE: every gated quantity is a function of the draw
  sequence. `n.threads` is pinned only so the ungated runtime footer is
  comparable; wall clock enters no criterion.
- `combineChains = FALSE` in every fit (shapes become chains x draws x
  pars). bart2 defaults to TRUE, and every gated statistic needs chains
  separated.
- `n.burn = 0L` everywhere; the whole trajectory is the datum.
- Seeds per `grouped-mixing.R`: data `set.seed(BASE_SEED + s)`, sampler
  `seed = s`; both arms of a pair share `s`. Two values, not one.
- Thresholds FREEZE at the end of Stage 0, written into this file before
  Stage 1 runs.

## Metrics

| id | quantity | contrast (positive = warm worse) |
|---|---|---|
| M0a | between-chain SD of the t=0 test fit over the gold posterior SD: rho0 | premise gate |
| M0b | mean pairwise L1 between per-chain inclusion vectors at t=0 | reported |
| M0c | modal-root-variable share across chains, per tree, meaned | premise gate |
| B1 | per-iterate test RMSE vs true mu at t in {1,25,100,250} | cold - warm, relative |
| B2 | per-chain first-hit B: least t with R(s) <= band_hi on [t, t+50] | log2(B_warm/B_cold) |
| C1 | mean per-iterate test RMSE over the last 500 iterates | relative gap |
| C2 | test RMSE of the pooled posterior mean, last 500 iterates | relative gap |
| C3 | pointwise 95% CI coverage of true mu, 200 held-out points | cold - warm, absolute |
| C4 | between-chain SD of per-chain time-averaged inclusion, meaned over columns | relative gap |
| S1m | between-chain SD of q_c, the per-chain fraction of iterates with ROOT split on x1 | difference |
| S2m | per-chain root-variable switch count | difference of log(1+switches) |
| X1 | mean 95% interval width | reported (the width confound) |
| X2 | R-hat and IACT on per-iterate test RMSE and on sigma | reported, UNGATED |

`band_hi` = 95th percentile of the COLD arm's R over its last 500
iterates in that cell, pooled over replicates and chains - defined at
every n including 50000, where the memo's gold-derived band did not
exist. Censoring pre-registered: no qualifying t within N_SAMPLES gives
B = N_SAMPLES + 1.

S1m/S2m read the root from `extract(fit, "trees")`, whose rows are
depth-first over (chain, sample, tree, n, var, value), so row 1 of each
group is the root. `keepTrees = TRUE` is affordable only at m = 1, which
is also the only place the statistic has power.

**R-hat is not gated, by measurement.** The probe's per-replicate SD of
(R-hat_warm - R-hat_cold) is 0.198 at n=2000 and 0.482 at n=20000, at 8
chains, against absolute R-hats of 1.28 and 1.53 - neither arm converges
on that functional in 1000 draws. A 0.05 margin needs ~1500 replicates,
and more chains does not rescue it. C4 measures the same concept by
pooling p columns instead of collapsing to one scalar (probe SD
0.0011-0.0015 against a cold level of 0.0075). Memo G5/K5 are deleted.

## Power, thresholds, and the multiplicity design

Each margin clears two bars: a SCIENTIFIC ceiling (largest difference we
would call benign) and a STATISTICAL floor of 4 x SE at the planned R,
buying ~0.99 per-test power at D = 0. Floor above ceiling means the
metric leaves the gate and is reported descriptively - the rule that
removed R-hat.

Probe SD/rep is quoted small n / large n where the two differ; R is 24 at
n <= 5000 and 20 at n = 50000.

| contrast | probe SD/rep | margin delta | 4 x SE (small / large) |
|---|---|---|---|
| C1 | 0.021 / 0.024 | 0.030 rel | 0.017 / 0.021 |
| C2 | 0.052 / 0.046 | 0.060 rel | 0.042 / 0.041 |
| C3 | 0.014 / 0.0069 | 0.020 abs | 0.011 / 0.006 |
| C4 | 0.0011 / 0.0015 | 0.30 rel (0.0023) | 0.0009 / 0.0013 |
| S1m, S2m | Stage 0 | from the two nulls | Stage 0 |

Every 4 x SE sits strictly inside its margin, which is the registration
condition. S1m/S2m margins are set in Stage 0 from computable nulls:
perfect mixing gives between-chain SD ~ sqrt(0.25/n_eff), perfect
collapse ~0.5.

- **Stratify on n**, pre-registered: SMALL (n <= 5000), LARGE (n = 50000),
  independent verdicts. The probe found C1 = +0.4% at n=2000 and +4.2% at
  n=20000; pooling would average away the only place the cost appears,
  and stratifying makes a scoped recommendation expressible.
- **GREEN per stratum**: (a) BENEFIT - one aggregate one-sided test on B1
  at t=100, alpha 0.05, warm ahead; (b) for each of C1-C4, an
  inverse-variance-weighted AGGREGATE non-inferiority test across that
  stratum's cells at margin delta, alpha 0.05; (c) no cell triggers KILL.
  Aggregating rather than conjoining per-cell tests is what fixes memo
  G4/G5: SE_agg = SE_cell/sqrt(#cells), so delta/SE_agg is 5.7 to 8.9 and
  power at D = 0 is >= 0.9999 per test.
- **KILL per stratum** requires BOTH, in >= 2 cells: a one-sided harm test
  (H0: contrast <= 0) rejecting under Holm across all C1-C4 x cells
  contrasts at family-wise 0.05, AND a point estimate exceeding delta. One
  cell rejecting triggers a fresh-seed re-run of it; confirmation counts
  as the second.
- **YELLOW** otherwise: tables to VD. A SMALL-GREEN / LARGE-KILL split is
  an expected outcome, not a design failure; the fallback surface is
  documenting n.grow.sweeps as recommended in a stated regime.

Designed rates: P(GREEN | contrast 0 everywhere) >= 0.99; P(KILL |
contrast 0 everywhere) <= 0.005; P(KILL | contrast >= 1.5 delta in >= 2
cells of a stratum) >= 0.95. At exactly delta the design returns YELLOW -
delta is the indifference point.

## M0: the premise gate

Runs first, ~8 min, and can invalidate everything after it. Seed
independence is not output dispersion: `grow-from-root.md` calls the
construction "greedy-stochastic" and at large n the per-node likelihood
ratio is exp(O(n)). Branches on rho0 (M0a) and M0c:

1. rho0 >= 1 - at least as dispersed as the target; Gelman-Rubin
   satisfied; Stage 1 runs as written.
2. 0.3 <= rho0 < 1 - under-dispersed, not degenerate. Stage 1 runs and a
   GREEN verdict carries a MANDATORY caveat: warm starting spends most of
   the cold init's diagnostic headroom, so the flip ships with the
   rho0-vs-k curve in `grow-from-root.md` and the Rd may not claim
   improved convergence diagnosis.
3. rho0 < 0.3 OR M0c > 0.9 - effectively a single point; dbarts IS Tan et
   al.'s single-optimum case and the cost family would measure a
   different object than the flip would ship. STOP and re-plan around a
   dispersion-injecting variant (per-chain k jitter, subsampled grows,
   warm-starting a subset of chains), recording the finding: with
   point-mass starts, defaulting the warm start silently converts four
   chains from a convergence diagnostic into four correlated replicates.

Probe expectation is branch 2 (t=0 between-chain SD at n=10000: 1.92
cold, 0.81 at k=2, 0.53 at k=8, against a posterior SD near 0.35).

## Invariance diagnostic (replaces memo K0)

One long warm and one long cold chain at 40000 draws on S1 n=5000 must
agree on sigma (|z| <= 3 on the posterior-mean difference over its MC SE)
and on the held-out fit (mean |fhat_w - fhat_c| over the posterior SD
<= 0.1). Disagreement DOES NOT stop the study - posterior invariance is
what a sticky sampler cannot demonstrate in finite time, so a
disagreement is indistinguishable from the study's most interesting
finding. Discriminator: re-run at 4x length. Shrinking like 1/sqrt(length)
is slow mixing, REPORTED as a primary result; stable in length is a bug,
and only then does the harness stop.

## Scenarios and run plan

200 held-out points from the same DGP everywhere; `n.trees = 75` unless
stated; 8 chains; 2000 draws.

| id | DGP | cells | role |
|---|---|---|---|
| S1 | Friedman, p=10, sigma=1 | n 500/5000/50000 | gated (5000 SMALL, 50000 LARGE); 500 in Part A only |
| S2 | Friedman, p=10, sigma=5 | n 5000 | gated SMALL |
| S3 | y = 3(x1+x2+x3)+eps, p=10 | n 5000, 50000 | gated. Linear is one of the two joint-largest BART->WS coverage gains in He-Hahn Table 4 (+0.22); the memo's 0.50 cell was XBART, not BART |
| S4 | Friedman, p=50 (45 noise) | n 5000 | gated SMALL; the over-commit-to-noise regime |
| S5 | y = 4*XOR(x1>.5,x2>.5)+eps | n 5000 at m=75 and m=1 | m=75 gated SMALL; m=1 carries S1m/S2m |
| S6 | x1 == x2 exactly, y = 5*x1+eps | n 5000 at m=1 | NULL CONTROL for S1m/S2m; ungated |
| S7 | 5 factors x 8 levels carry the signal, 5 ordinal noise | n 5000 | DESCRIPTIVE (D1); prediction is a null |

S5, not S7, carries the structure family because it is the credible harm
case: at XOR every depth-1 cut has near-identical integrated likelihood,
so grow-from-root commits to a root the MH sweeps must undo, and the two
modes (root on x1 / x2) are exactly exchangeable with the analytic answer
1/2. A probe at m=1 confirmed power - three chains switched roots 8-12
times in 300 draws while a fourth locked onto x2 for 272 of 300,
between-chain SD 0.198 against a mixing null near 0.05. S6's exact
duplicates are NOT a stickiness probe (a change move between identical
columns has likelihood ratio exactly 1, so the modes are trivially
connected), hence the demotion to null control. The m=75 arms cannot
detect structural mode collapse at all - the ensemble self-averages the
label at rate 1/sqrt(#splits) with no mixing required - so the study says
so rather than gating a powerless statistic (memo G3/K3 deleted), and
m=75 cells are gated on function space plus C4.

Arms. Part A (k selection): k in {0,1,2,4,5,8,16}, 600 draws, 4 chains,
R=8, on S1 x 3n / S2 / S3@5000 / S4. k=5 is in the grid as the only
published ecosystem default (stochtree `num_gfr = 5`; that paper's "10"
is wrong about its own package). Part B: cold vs warm(k*) only, k* = the
smallest k whose B2 is within 10% of the grid minimum in every Part A
cell. The whole B-vs-k curve is reported regardless: the R5 method's 2L
is unmeasured and the flip picks a new number rather than ratifying it.
Gold standard (C4's reference, M0a's denominator): cold, 8 chains, 10000
burn + 10000 kept, 3 reps, at n = 5000 cells only; every gate consuming
it names only those cells.

## Steps

1. `benchmarks/R/grow-init.R`: header stating the LOAD-INDEPENDENT claim,
   tunables, S1-S7 constructors, `fit_arm` with `combineChains = FALSE` /
   `n.burn = 0` / pinned `n.threads`, metric extractors. `QUICK = TRUE`
   runs Stage 0 only.
2. Stage 0 (~25-30 min): M0 over the k x n grid at R=20; the m=1 S5/S6
   cells at R=8; every C1-C4 contrast at R=8 on the gated cells; the
   invariance diagnostic.
3. FREEZE. Write Stage 0's measured SDs into the power table, confirm or
   revise each delta against the floor-above-ceiling rule, set S1m/S2m
   margins from their nulls, resolve the M0 branch. Any metric whose
   floor exceeds its ceiling leaves the gate here, in writing.
4. Stage 1 Part A (~10 min): the k grid, B1/B2, select k*.
5. Stage 1 Part B (~100 min): cold vs warm(k*) at full R, plus S7
   descriptive; gold standard (~15 min).
6. Report: per-cell tables printing each contrast beside its frozen
   margin and Holm-adjusted p, the two stratum verdicts, the B-vs-k and
   rho0-vs-k curves, X1/X2 as context.

## Verification

- `QUICK = TRUE` completes in ~25-30 min printing every Stage 0 SE; full
  run ~2.5-2.75 h, one invocation, any load.
- Harness shape assertions: under `combineChains = FALSE`,
  `dim(fit$yhat.test)` is chains x draws x 200 and `dim(fit$varcount)` is
  chains x draws x p. Stop if not.
- S6 m=1 null control returns pooled p1 within MC error of 1/2 in BOTH
  arms with non-zero switch counts in every chain; otherwise the S1m/S2m
  estimators are wrong and the structure family is void.
- Two invocations at the same BASE_SEED and different `n.threads` produce
  identical tables.
- The script exits nonzero only on the invariance diagnostic's bug branch
  (disagreement stable in chain length), never on a study finding.

## Results (run 2026-08-08; harness, results.md, and all 98 checkpoints
## are untracked, not retained)

VERDICT: KILL in BOTH strata, each on a mandated fresh-seed re-run.
The grow-from-root warm start does NOT default; n.grow.sweeps stays
opt-in. Reopening requires a changed design (for example warm-starting
a subset of chains, or post-init tempering) plus a re-run of this
battery on an extended k grid.

- M0 landed branch 1, not the anticipated branch 2: rho0 >= 1 in every
  grid cell (minimum 1.116 at n=500/k=16; 1.65-4.06 at n=50000), modal
  agreement max 0.442. The warm start stays overdispersed relative to
  the target everywhere measured; Stage 1 ran with no caveat.
- SMALL: aggregate benefit passed (B1 at t=100 +0.199, z=25.5) and all
  four aggregate non-inferiority tests passed, but S4-5000 (Friedman
  p=50, noise-heavy) carried C2 = +11.10% against the 6% margin (Holm
  p=2e-17); re-run +11.58% (z=10.84). KILL.
- LARGE: aggregate benefit passed (+0.152, z=27.9), all aggregates
  non-inferior, but S1-50000 carried C1 = +4.47% (margin 3%) and C2 =
  +10.66% (margin 6%), both Holm-significant; re-run C1 +3.45%
  (z=6.33), C2 +9.87% (z=7.87). KILL.
- The substantive finding: the plateau cost concentrates in specific
  regimes - noise-heavy designs and large n - and pooling averages it
  away. The per-cell KILL clause, not any aggregate, is what caught it.
- S7 matched its pre-registered null (all contrasts inside the margins);
  the categorical renormalization stays closed as documented behavior.
- The S6 duplicate-column null control FAILED its non-zero-switch
  clause: 147/192 warm chains had zero root switches, because dbarts'
  change move cannot rotate the root of a deep tree. Per the
  Verification clause the S1m/S2m structure estimators are VOID (they
  had also independently left the gate on floor-above-ceiling). The
  root-rotation limitation is itself a stickiness datum, recorded here
  for any future structure-metric design.
- k* = 16 is a grid artefact: no k satisfied the within-10%-everywhere
  selection rule, and 16 is the largest grid value with B2 still
  falling. Any future flip proposal re-runs Part A on an extended grid.
- The invariance diagnostic landed on the slow-mixing branch and the
  sigma difference flips sign between 40k and 160k draws (+0.00138 ->
  -0.00188): mixing noise, not init-dependent bias.
- 17 deviations logged in results.md, plus a process note: an
  intermediate analysis pass printed GREEN before the mandated re-runs
  had run; the verdict logic was corrected to a PENDING state and the
  GREEN was never a valid result. Total compute 1.93 h.
