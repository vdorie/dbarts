# approximate-draws: may a Polya-Gamma latent whose exact cost grows with its shape be drawn another way?

Status: DESIGN, revised 2026-11-03 after its blind critique (verdict: revise first). No study fit has run. Part A,
an exact bounded-cost route checked in place, is the gate; the approximation study (Parts C and D) runs only
where Part A leaves a question. Admitting any approximate draw stays the maintainer's call (dec-B567).

agent: opus for Part A's engine scratch build, its checks and the harnesses' statistics; one opus review of
Part A's code and checks before its verdict; if C or D run, one opus review of their harnesses and of the
reading, and two independent referees on the report if a default would change ([Studies](README.md#studies)).
rng: NEUTRAL for the study: scratch builds and harnesses that never land. What it recommends lands as its own
slice, shifting for the exact route (draws change, the posterior does not), posterior-changing for an
approximation, both above a threshold that keeps every recorded draw bitwise.
window: after 1.0-0 unless the maintainer pulls it in.
budget: Part A about 1,000 lines (scratch engine C and C++ about 750 with tests, R checks about 250) and about
11 core-hours; the whole arc, if C and D run, about 2,500 lines and 35 core-hours, stop at 50
([Compute](#compute)).

## Summary for the maintainer

The question: when a row's hidden Polya-Gamma value is expensive to draw exactly because the row carries many
trials, a large frequency weight or a large count, can it be drawn in time that does not grow with that number;
and where the only fast way is approximate, is its effect on what users read small enough to state and accept?

Why it matters. The exact draw is a sum of one small draw per trial, about 160 nanoseconds each. A logistic fit
of 200 rows spends about 2 microseconds a row per sweep on everything else, so at 100 trials a row the sweep is
9 times slower and at 10,000 trials 800 times. Your binomial ruling makes large trial counts ordinary: county
disease rates over thousands of births, school results over hundreds of pupils. Negative binomial fits with
counts in the thousands pay the same way. Real-valued weights and a real-valued negative binomial shape wait on
a draw for a shape that is not a whole number.

What the design found:

- An exact route whose cost does not grow with the trials. Propose the hidden value from a fast approximation
  that is already very close (the first few terms of its exact series drawn exactly, the rest replaced by one
  gamma with the rest's exact mean and variance), and keep the proposal with a probability computed from the
  exact density, otherwise keep last sweep's value, as tree proposals do. At 16 to 1,000 trials at least 9,999
  proposals in 10,000 are kept. The fit's long-run distribution is exactly the posterior. Evaluating the
  exact density takes about 1.2 to 1.5 microseconds whatever the trials, measured on a heavily loaded machine.
- That route is not exact everywhere without care. In a negative binomial fit with an estimated shape, the shape
  is drawn by a step that does not look at the hidden values; keeping a previous value is valid only on a
  sweep where the shape did not move, and a toy check confirms both halves of that. The same holds when you, or
  code running BART inside a larger sampler, change weights, offsets or responses between sweeps. In those
  cases the affected rows draw exactly, as today.
- The same close approximation, used on its own without the accept step, is off by a few parts in a million at
  64 trials and less above. Used directly it would be faster than the exact route by an estimated 1.6 to 2.3
  times per sweep at 75 trees, less with more trees. The other approximations in the literature are worse or
  dearer: a plain normal is off by 3 percent at 64 trials, and the saddle-point sampler costs several times
  more for a similar accuracy.

The plan: Part A builds the exact route in a scratch copy of the engine, times whole sweeps, and checks it is
right (against the exact sum, a joint-distribution check, an exact single-tree posterior and simulation-based
calibration). If it is within your speed tolerance of the fastest approximation at every shape, and its
acceptance holds up, the approximation study does not run, and what goes to you is whether to land the exact
route. About 11 core-hours. If not, Parts C and D test the approximations on exact posteriors and on real
covariates with real trial counts at equal running time, about 24 core-hours more, and a short research report
follows if a default would change.

Three calls are yours: whether the accept-or-keep step counts as an exact draw under your ruling that every
draw stays exact (recommended: yes, within the limits above); the speed tolerance that decides whether the
approximation study runs at all (recommended: 2 times); and, only if Part A leaves it open, whether an
approximation that passes is a default or an opt-in.

## Goal

A measured answer, per kind of shape (whole shapes above a threshold; real shapes of 2 and above; real shapes
below 2), to which draw dbarts should use: today's exact sum, the exact bounded-cost route, or an approximation
above a stated threshold; with each reading fixed before the fits it reads.

## Context

Rulings, quoted, not reopened:

- dec-B14, quoted whole under
  [Open call 1: does the accept-or-keep step count as an exact draw?](#open-call-1-does-the-accept-or-keep-step-count-as-an-exact-draw).
- dec-B143: the large-shape cost is documented for 1.0-0; the maintainer, 2026-09-29: "I don't know
  necessarily that we'll want to build it, we may instead want to permit approximations under some reasonable
  criteria for when to use them."
- dec-B567: the binomial response joins logistic, one draw of shape n per row; the maintainer, 2026-11-03:
  "there's a much higher likelihood that users will want to use it with large numbers of trials." "Whether any
  approximate draw is admitted remains the maintainer's."
- dec-B563: frequency weights reach the Polya-Gamma families as the shape times w, "a sweep's cost grows with
  it"; dec-B568: real weights where the draw is exact, the Polya-Gamma families whole numbers "until a
  real-shape draw exists"; dec-B566 and dec-B570 fix what a weighted or zero-weight row draws and means, which
  no route here changes.
- dec-B15: integer dispersion, "fully exact"; dec-B202: trials and fixed shapes capped at one million ("if
  something just isn't going to run due to wall time constraints it's hard to argue that should support it").
  A bounded draw removes that reason; revisiting the cap is a consequence, not in scope.

Engine: the shape seam is [`simulatePolyaGammaShape`](../../src/bartcore/model.hpp), a sum of round(b) draws
of [`ext_rng_simulatePolyaGamma`](../../src/external/random.c) (Devroye's method), used by
[`NBResponse::drawOmega`](../../src/bartcore/model.hpp) and the bridge's
[`drawAugmentation`](../../src/R_interface_bartcore.cpp). Logistic and hazard.logistic weights
([`LogisticResponse::refreshLatents`](../../src/bartcore/model.hpp)) and multinomial trials
([`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp)) sum in their own loops. An
inactive row's draw is skipped; a reactivated row and a row whose count is swapped are redrawn exactly; a
cold start is deterministic, at the latent's mean ([`NBResponse::coldStart`](../../src/bartcore/model.hpp)). With
an estimated shape, [`NBResponse::refreshLatents`](../../src/bartcore/model.hpp) draws r from its grid
conditional collapsed over omega, then omega at the new r. The binomial row of dec-B567 is a logistic row of
shape n times w (docs/design/frequency-weights.md on branch frequency-weights).

Background: [nbinom-large-shape-pg.md](nbinom-large-shape-pg.md) designs an exact bounded-cost construction
(the saddle-point envelope accepted against the certified density, resting on WPS Conjecture 15), records that
the saddle-point sampler is approximate and that BayesLogit's fails at large shape. Its option (c), a
Metropolis step, was set aside because "the chain stays exact but the draw is not, and the bridge's replay has
no previous omega".

### Where large or real shape arises, and today's cost

Measured 2026-11-03 on the laptop (arm64, load 5 to 9, per-round minima of 2 or 3 rounds), shipped build of
main; standalone C timing of the engine's primitive linked against that build.

| source | shape per row | example users meet | cost today |
|---|---|---|---|
| binomial trials (dec-B567) | n w | NC SIDS: births per county 248 to 21,600; schools: 101 to 4,120 pupils | 160 ns per trial |
| logistic, hazard.logistic frequency weights (dec-B563) | w | identical rows collapsed: Titanic's 24 cells to 670; a census collapsed over discrete covariates | the same |
| nbinom | (y + r) w, r estimated at most 50 or fixed | Migration (carData): counts to 99,400, 830,000 in all | the same |
| multinomial trials | n_i per category | votes over several parties | 480 ns per trial over 3 categories |
| real shape | real w (power priors, fractional counts); real r (five of the nbinom study's twelve datasets have a shape below the grid's floor of 1) | | refused today |

- Primitive: PG(1, psi) 146 to 170 ns (psi 0 to 8).
- Logistic fit, 200 rows, 75 trees, one chain: 2 us a row per sweep at w = 1; 18 at 100; 166 at 1,000;
  1,630 at 10,000. nbinom fixed r: 4 at 10, 160 at 1,000, 1,590 at 10,000; estimated r on counts averaging 224,
  36. Multinomial, 3 categories: 3.3 at 1 trial, 50 at 100, 480 at 1,000.
- Beyond about 10 trials a row the latent is the sweep. A default fit (4 chains, 1,000 sweeps each) of 1,000
  rows of 10,000 trials is about 1.8 core-hours; at 100,000 trials, 18.

## Candidates

"TV" is total variation distance from PG(b, psi), computed against the exact density by Fourier inversion of
the Laplace transform (probe laterr.R, accurate to about 1e-13), the largest over psi in {0, 2, 8}. Costs per
draw on the laptop; a gamma draw is about 41 ns.

| candidate | source's account of its error | TV at b = 64 / 1e3 / 1e4 | cost per draw | shapes |
|---|---|---|---|---|
| today: sum of b Devroye draws | exact (Polson, Scott and Windle 2013) | 0 | 160 ns times b | whole b |
| truncated sum plus exact-moment remainder, K terms: the first K of the series PG = sum_k g_k / d_k drawn exactly, the rest one gamma with the rest's exact mean and variance (the PG moments less the first K terms') | constructed here; the series is PSW 2013 eq. 1 | K = 3: 3.1e-5 / 7.6e-6 / 2.4e-6; K = 5: 3.2e-6 / 8.0e-7 / 2.5e-7; K = 10: 1.2e-7 / 2.9e-8 / 9.1e-9; K = 200: 1e-13 or less (rounding) | K + 1 gamma draws: about 250 ns at K = 5 | any real b; at b = 0.5, K = 10 gives 4.2e-5, K = 20 1.4e-6 |
| saddle-point sampler (Windle, Polson and Scott 2014, sec 5) | "generates approximate J*(n, z) random variates"; no bound stated | 2.2e-5 / 3.5e-7 / 1.1e-8 | about 9 Devroye draws by WPS Table 1, about 1.5 us (not run here) | b >= 13 in BayesLogit's hybrid |
| moment-matched gamma | constructed here | 7.7e-3 / 1.9e-3 / 6.1e-4 | about 40 ns | any b |
| moment-matched normal (WPS sec 6; BayesLogit above 170; Python polyagamma above 50) | "The central limit theorem suggests that this is a reasonable approximation"; polyagamma: "equally accurate" as the saddle point | 3.1e-2 / 7.8e-3 / 2.5e-3 | 35 ns | large b; can go negative |
| truncated gamma sum, K = 200, no remainder (PSW 2013; BayesLogit) | "quite accurate for such small values of the first parameter ... may not be the case ... for arbitrary b" | 2.6e-2 / 0.10 / 0.32 (grows with b) | 8 to 13 us | small real b only |
| alternate sampler (WPS sec 4) | exact for h >= 1 "if Conjecture 8 holds", verified on [1, 4]; acceptance falls with h | 0 | linear in b | a real remainder in [1, 2) |
| auxiliary mixture (Fruhwirth-Schnatter et al. 2009, cited by PSW 2013) | "a fast approximation, at the cost of introducing a more complex latent-variable structure" | - | - | replaces the augmentation, not this seam |

The exact routes, each exact given its stated conditions:

| route | what it needs | cost per row per sweep | status |
|---|---|---|---|
| A: independence Metropolis on omega, proposal the K = 3 truncated sum plus remainder | the exact density f of PG(b, psi) at the proposal, and the proposal's density q there (both by inversion, sharing the K finite factors); f(omega) and q(omega) at the current value cached per row; psi enters f and q as an analytic tilt | measured: one f evaluation 1.2 to 1.5 us at b 16 to 1e5 (33-node trapezoid, load 19 to 25, the critique's pgdens.c); a q evaluation in the same pass, cost to measure; 4 gamma draws | Part A builds it |
| the exact construction (nbinom-large-shape-pg.md) | the saddle-point envelope, the certified density, WPS Conjecture 15 closed | estimated 2 to 4 us | not built; kept as the fallback for independent exact draws |
| Metropolis on the remainder only (first K terms drawn exactly, the rest proposed from its moment-matched gamma) | one density, of the remainder | one inversion | rejected: its weight f/q reaches 11 to 330 at b 16 to 64 and is unbounded at b = 1 (probe remmh.R) |
| Metropolis with the moment-matched gamma as proposal | one density | one inversion | rejected: the gamma's tail rate (6 at psi = 0) exceeds PG's (pi^2/2 + psi^2/2 = 4.93), so f/q is unbounded (220 on the grid at b = 16, psi = 8) and the step is not geometrically ergodic (Mengersen and Tweedie 1996) |

Route A's proposal has a bounded weight: its first term carries the exact tail rate, and on the grid sup f/q is
at most 1.004 with acceptance at least 0.9999 at b 16 to 1,000, psi 0 to 8 (the critique's probe-mh.R), so
the independence step is uniformly ergodic. "Certified" means the accept decision is made from a density whose
error is bounded down to the floating-point rounding floor, 1e-11 to 1e-9 at b >= 64 with 33 nodes; below b
of about 64 the 33-node rule is good only to 1e-7 (b = 16) and needs more nodes, which Part A measures. That
is the standard Devroye's sampler meets: its series is also evaluated in floating point.

Findings that shape the design:

- The approximation that matters is the truncated sum with an exact-moment remainder: at K = 5 it is 7 times
  closer than the saddle point at b = 64, for a sixth of its cost, and its error falls with b. The normal's
  claimed parity with the saddle point does not hold (1,400 times further at b = 64); the moment-matched gamma
  dominates the normal at equal cost; the truncated sum without its remainder is no large-shape candidate.
- Against route A, the best approximation saves the density evaluations and little else, so the question
  becomes the speed ratio at the sweep level: about (2 + 1.5 to 3 + 0.17) / (2 + 0.25), 1.6 to 2.3 at 75
  trees and 200 rows, lower with more trees. Part A measures it in place.
- Perturbation bounds do not settle an approximation's effect. Mitrophanov's bound (Alquier, Friel, Everitt
  and Boland 2016, Theorem 2.1) limits the chain's error by the kernel's TV times a mixing constant, and a
  sweep's kernel redraws n latents: at 1,000 rows and 1e-3 a row the bound is vacuous.

Toy pilot (probe toy.R): one logit parameter shared by n rows of b trials at p = 0.2, prior N(0, 9), 2e5
Gibbs iterations against the posterior by quadrature; Monte Carlo error about 0.003 posterior sd.

| n, b | normal: mean shift, sd ratio | gamma: mean shift, sd ratio |
|---|---|---|
| 1, 10 | -0.048, 1.105 | -0.010, 1.018 |
| 1, 64 | +0.001, 1.007 | -0.001, 1.002 |
| 1, 1000 | -0.000, 1.002 | +0.003, 1.000 |
| 20, 10 / 64 / 1000 | all within 0.005 and 0.2 percent | all within 0.003 and 0.3 percent |

### Where route A is exact, and where it is not

A Metropolis step on omega leaves the target invariant only if the current omega is a draw from omega's
conditional given everything the step conditions on. Inside a sweep that holds: the trees are drawn given
omega, then omega given the trees. It fails after a step that changed omega's conditional without looking at
omega, which is the situation van Dyk and Jiao (2015) analyze: in a partially collapsed sampler, "reducing the
conditioning in an MH within Gibbs sampler can change the stationary distribution".

- Estimated-r nbinom: r is drawn collapsed over omega, then omega at the new r. When r moved, the old omega was
  drawn at the old r and a rejection would keep it: route A is invalid there, and the row is drawn exactly.
  When r did not move, the current omega is a draw from its conditional at that r, because the r step never
  reads omega, so the event that r stayed is independent of omega given the trees and r; route A is valid.
  This plan departs from its critique here, which held the failure "even when r does not move". Evidence: the
  argument above, and a toy partially collapsed sampler (probe pcg.R: r in {1, 4}, omega given r normal, a
  parameter given omega normal, 1e6 sweeps, three seeds): the always-Metropolis arm biases the parameter's
  variance by about 2 percent (4.33 to 4.34 against 4.25), the moved-then-exact arm matches it (4.247 to
  4.255) as the exact arm does. Part A's joint-distribution check runs both arms on the engine. The exact draw
  when r moves costs linear in y + r; how often r moves is measured in Part A, and if it moves on most sweeps,
  estimated-r nbinom keeps the exact sum.
- Changes from outside a sweep: setWeights, setOffset, setResponse and setPredictor change a row's conditional
  without reading omega, by the same argument. Under route A the first refresh after such a change draws the
  affected rows exactly; a sampler whose weights change every sweep pays the linear cost each sweep, as today.
  A reactivated row is redrawn exactly, as today. A cold start is a deterministic value, not a draw; route A
  starts from it with its density cached, which is ordinary burn-in.
- The augmentation door the code running BART inside a larger sampler calls has no previous value and keeps
  the exact sum.
- Real b: the density at b below about 64 needs more nodes, and below 2 more again (its tails decay like
  exp(-b sqrt|y|)); Part A measures cost and accuracy at b in {0.3, 0.7, 1.5, 2, 4, 16}. Where the density can
  be certified at bounded cost, route A gives an exact draw at real b, below 1 included, for which no exact
  sampler is published.

## Part A, the gate

On a scratch build carrying route A at the seam above a compile-time threshold of 64 (so every recorded draw
stays bitwise), routing the logistic and multinomial loops through the seam:

1. Density and proposal: the trapezoid inversion of f and q in one pass, with its error bound, at b in
   {0.3, 0.7, 1.5, 2, 4, 16, 64, 1e3, 1e5, 1e6}; node counts per range set so the bound reaches the rounding
   floor. Acceptance over b x psi in {0, 2, 8, 20}.
2. Sweep-level speed, in place: logistic with binomial rows and nbinom fixed and estimated r, 200 and 2,000
   rows, 75 and 200 trees, b in {64, 1e3, 1e4, 1e5}: per-sweep time of route A, of today's sum and of the K = 5
   approximation built as a third arm of the same scratch build. Same machine, alternating arms, per-round
   minima, load stamped.
3. Correctness, each with the exact sum as its control passing first:
   - the chain's omega at fixed (b, psi), thinned, against the exact sum by two-sample KS (b <= 1e3, 1e5
     draws), and moments to the third cumulant against the exact formulas at 1e4 and above;
   - the joint-distribution check (Geweke 2004) as [geweke-mc.R](../../benchmarks/R/geweke-mc.R) runs it, for
     logistic with binomial rows, fixed-r nbinom, estimated-r nbinom with the moved-then-exact rule, and as
     positive controls the always-Metropolis estimated-r arm and a build whose cached density goes stale after
     psi changes, both of which must fail;
   - the exact single-tree posterior of [logistic-reference.R](../../benchmarks/R/logistic-reference.R) part
     1 with rows of b trials, within that gate's own tolerances;
   - simulation-based calibration (Talts et al. 2018) through [sbc.R](../../benchmarks/R/sbc.R)'s logistic
     machinery with binomial rows at b = 64 and 1e3, 300 replications, its exact arm passing first so BART's
     own mixing is not read as route A's fault.
4. The share of sweeps on which estimated r moves, on the nbinom datasets of the earlier study and Migration.

Part A's verdict: route A passes if every check in step 3 passes, its acceptance is at least 0.999 at every
measured (b, psi) where it is used, and its sweep is within the maintainer's speed tolerance (Open call 2) of
the K = 5 approximation's at every measured setting. If it passes, the approximate branch ends for the shapes
it covers, and what goes to the maintainer is the landing of route A (a shifting slice). Parts C and D run only
for shapes or families where it fails, with route A wherever it passed as their exact arm.

## Parts C and D, if Part A leaves a question

### How the effect is measured

On the posterior of what users read, against a known truth: exact posteriors (the one-parameter toy by
quadrature; single-tree BART by enumeration with rows of b trials), then semi-synthetic data on real
covariates and real trials at equal running time. Every probability is transformed per draw and averaged last.

### The reading that decides (fixed before any C or D fit)

The bar is set on the posterior scale. Under a calibrated posterior, a shift of the posterior mean by d sd
moves 95 percent coverage by about 0.11 points at d = 0.1 and 7e-5 at d = 0.025, while a width error of
2 percent moves it by 0.44 points at 95 percent and 0.67 at 90. Coverage is the scale the maintainer judges
sampler changes on (calibration, over-coverage counting as error); a mean shift also enters a point estimate's
error, and a shift of 0.1 sd adds 0.5 percent to the RMSE of a calibrated posterior mean (sqrt(1.01)).

Reading 1, admissibility of candidate X at shape b, in every C cell at b:
- the exact posterior's mass inside X's 90 and 95 percent intervals within 0.5 points of nominal;
- X's posterior mean within 0.1 exact-posterior sd of the exact mean;
- each cell's equivalence interval at level 1 - 0.05/m, m the number of cells, so the chance of any false
  failure across the family is at most 5 percent.

These bars are this plan's, argued as above: half a coverage point is what a coverage study of about 2,000
independent intervals resolves at one standard error (sqrt(0.95 x 0.05 / 2000) = 0.005), so an effect under it
cannot show in coverage evidence the package or a user would produce; a 0.1-sd shift moves coverage by a fifth
of that. They are absolute, so they bind equally at any run length; a fixed bias weighs more against Monte Carlo
error as runs grow, not less. X's threshold b* is the smallest scanned b from which every larger scanned b
passes.

Reading 2, at the user level, in every D cell: the paired difference from the exact arm in coverage (90 and
95, per unit, per group, per contrast) has its 3-standard-error band inside +-0.01, and its mean interval
width within 2 percent. RMSE and bias against the truth are reported for both arms, not gated: model error,
which both arms share, dominates them. A candidate passing Reading 1 and failing here is not admitted at that
shape, and the report says what Reading 1 missed.

Recommendation rule: recommend admitting X above b* where Readings 1 and 2 pass, with its measured sweep gain
over route A beside it; the maintainer rules (Open call 3).

### Designs and data

C1, the toy: candidates K = 5 truncated sum plus remainder, saddle point if route A failed on cost, the gamma
and the normal as comparators; b in {16, 32, 64, 128, 256, 1e3, 1e4, 1e5}; success probability 0.5, 0.05,
0.001 (at b = 1e4 and 0.001 a row holds about 10 events, psi near 7); n = 1 and 20 rows; real b in {0.3, 0.7,
1.5} with K = 20. The exact sum runs as a control at b <= 1e3 only (at b = 1e5, n = 20 it would cost about
18 core-hours a cell); above, route A is the control where Part A passed it.

C2, single tree: logistic-reference.R's one-predictor, three-cut design with binomial rows, b in {16, 64, 256,
1e3}, 2 rows per cell (near-isolated) and 25 (pooled); a fixed-r nbinom twin with counts of b. Sized at ESS 1e5
per cell, so an interval endpoint's Monte Carlo error (2.67/sqrt(ESS) sd for a 2.5 percent quantile) moves
coverage by about 0.05 points and the family-wise interval stays inside the 0.5-point bar.

D, real covariates and real trials. The datasets were chosen by the three occasions where large shapes arise,
aggregated binomial counts, rows collapsed into weighted cells, and large counts, taking for each occasion the
datasets available offline (R, its recommended packages, installed CRAN packages, adult through the k study's
loader) that carry it with real covariates; this is not a search over datasets for a result.

| design | rows | shape | outcome the truth is fitted to | estimands |
|---|---|---|---|---|
| North Carolina SIDS, 1974 (sf's nc) | 100 counties | births 248 to 21,600, rate about 0.002 | SIDS deaths of births | county rate; rate by region; rate ratio, top against bottom quartile of nonwhite births |
| California schools (survey's apipop), 1,000 sampled | 1,000 | enrollment 101 to 4,120 | pupils eligible for subsidized meals, round(meals x enroll / 100), of enroll; covariates without meals | school proportion; district means; contrast by school type |
| adult census, collapsed over its discrete covariates | thousands of cells | cell counts to the thousands | income above 50,000 | cell probability; probability by education and by sex; sex contrast within occupation |
| Titanic (datasets) | 24 cells | weights to 670 | survival | survival by cell; class contrasts |
| Migration (carData), nbinom | 90 | counts to 99,400 | migrants | flow mean; destination totals; distance ratio |

Truths: on each dataset's real covariates and real trials, two surfaces fitted to the real outcome, a BART fit
and a smooth additive model, so the answer does not rest on BART's own prior; outcomes drawn afresh per
replication. Arms: the exact arm (route A where it passed, else the sum), each candidate passing Reading 1, at
default length and at equal time (the exact arm's wall time, capped at 8 times the default length).

Real shape at the user level (the five nbinom datasets with a shape below 1) needs the real-shape r step of
TODO negbin-real-dispersion and is not in this study.

### Sizing

| component | question | size | cost (core-h) | precision reached | smallest decision-changing difference | stop |
|---|---|---|---|---|---|---|
| A1 | density cost, accuracy; acceptance | 1e6 evaluations per (b, psi) | 0.2 | timing to 5 percent; bound to rounding | acceptance 0.999; a bound above the floor | none |
| A2 | sweep speed in place | 4 b x 2 n x 2 trees x 3 arms, 5 rounds | 0.5 | ratio to 5 percent | the tolerance of Open call 2 | none |
| A3 | correctness | KS 1e5 draws; Geweke 7 arms; single tree; SBC 2 arms x 300 | 0.5 + 4 + 0.5 + 5 | Geweke about 0.05 sd at 3 se; SBC about 0.12 sd | a detection; the positive controls must be detected | a failing control stops Part A |
| A4 | how often r moves | the nbinom datasets, default length | 0.2 | share to 1 point | whether moved sweeps dominate the cost | none |
| C1 | posterior effect, toy | 2e5 iterations per cell, about 250 cells | 0.5, + 1 for the sum at b <= 1e3 | mean 0.012 sd, coverage 0.2 points, family-wise | 0.1 sd; 0.5 points | rerun a straddling cell at 4x, once |
| C2 | posterior effect, single tree | ESS 1e5 per cell, 12 cells | 10 | coverage about 0.3 points family-wise | 0.5 points | the same |
| D | user level, equal time | 10 replications per truth (pilot 3), to 20 if a band straddles | 12, 24 doubled | paired coverage band about +-0.005 estimated | 0.01 coverage | sequential at 10 |

The Geweke and SBC arms check route A's code (a stale cached density is the likely defect); they cannot see
an approximation's small error and are not part of Reading 1. D's paired band is an estimate the pilot
measures; the reviewer audits the table before any main run.

## Compute

Part A: about 11 core-hours (A3's SBC 5 and Geweke 4 dominate), the laptop at 2 cores except SBC and Geweke on
an x86-64 Linux host; about 1,000 lines. If C and D run in full: C about 11.5, D about 12 (24 doubled), so the
whole arc about 35 core-hours, 47 on the doubled path, stop at 50, and about 2,500 lines. D per replication:
exact arms about 1,200 core-seconds with today's sum (Migration 530, schools 400, SIDS 210, adult 70), less
with route A; candidates about 450 each.

## Report

If a default would change, a PDF of at most 5 pages of text, readable alone, refereed twice before it reaches
the maintainer: the problem and today's cost; each candidate with its source's account of its error and its
cost, and the exact routes; the evaluation (exact posteriors, semi-synthetic data on real covariates and
trials, the bars and their argument); the results (latent TV against b; coverage and mean effects against b;
coverage, width, RMSE and bias at equal time by dataset; speed); the recommendation and its limits. Otherwise
the finding goes in docs/design/approximate-draws.md in the plan-and-finding form.

## Decision

### Open call 1: does the accept-or-keep step count as an exact draw?

Background. Today each hidden value is a fresh draw from its exact distribution. The new route proposes a value
from a close, fast approximation and keeps it with a probability worked out from the exact density, otherwise
keeping last sweep's value; more than 9,999 proposals in 10,000 are kept. The fit's long-run distribution is
exactly the posterior, as it is today; what changes is that a hidden value is occasionally carried over one
sweep, which tree proposals already do whenever they are refused. In 0.9-34 none of these families existed;
its probit hidden values were fresh exact draws and its trees were accept-or-keep steps, as now.

Your ruling, whole: "Every draw in 1.0 stays exact: an approximate Polya-Gamma draw, the data augmentation that
would make a logistic likelihood conjugate under arbitrary real weights, is declined. The alternative was
accepting that approximation to unlock the models behind it. A real-valued negative-binomial dispersion and
arbitrary real weights on a binary fit therefore cannot be fit in 1.0-0, and both wait on that one unlock. The
same ruling kept the hurdle model's sampler-only mode refused, and struck, as reopening a settled point, a
proposal to remove the second, redundant sampler creation that ordinal and negative-binomial fits perform,
about 3 milliseconds a fit." It declines approximate draws for 1.0, which comes before this work; it does not
speak to accept-or-keep steps.

The route is exact only where nothing has changed a row's distribution behind its back: in negative binomial
fits with an estimated shape, on sweeps where the shape moved, and after weights, offsets or responses are
changed from outside, those rows are drawn fresh exactly, at today's cost.

- (a) It counts as exact, within those limits: logistic, binomial, multinomial and negative binomial with a
  fixed shape always; with an estimated shape on sweeps where the shape stayed, if Part A's check confirms it,
  otherwise always fresh. About 1,000 lines for Part A, a landing slice after it; the hidden-value cost about
  1.5 to 3 microseconds a row whatever the trials.
- (b) Only fresh exact draws count: the exact construction designed earlier, 1,000 to 2,000 lines, resting on
  a conjecture still to be proved, or an admitted approximation.

Recommendation: (a). What your ruling protects is a posterior computed exactly, and that holds; every BART draw
already depends on the one before it.

### Open call 2: how much faster must an approximation be before it is studied at all?

Background. With the exact route in place, the best approximation (the same close proposal used without the
accept step) still skips the density evaluation, so a sweep would be faster by a ratio Part A measures,
estimated at 1.6 to 2.3 times at 75 trees and 200 rows and less with more trees, since the trees' share of a
sweep grows. This is a value call: how much speed justifies giving up exactness.

- (a) 2 times: if the exact route's sweep is within 2 times the approximation's at every shape measured, no
  approximation is offered and Parts C and D do not run (about 24 core-hours, 1,500 lines of harness and a
  report saved).
- (b) A tighter tolerance, such as 1.25 times: the approximation study probably runs, and an approximation may
  be offered for a gain between 1.25 and 2 times.
- (c) Exact only, whatever the ratio: no approximation study for whole shapes.

Recommendation: (a). Every fit stays exact under one rule across families and the help needs no caveat, for at
most half a sweep's time at the estimated ratio.

### Open call 3, deferred until after Part A: an approximation that passes, default or opt-in?

This goes to the maintainer only if Part A leaves an approximation live. Background for then: BayesLogit 2.4's
default sampler switches silently to a normal above a shape of 170, to the saddle point above 13 and to a sum
of 1,000 gammas from 3 to 13, and its manual does not say so (read in its source by this plan's critique, not
re-read here); the Python polyagamma package states its switch to a normal above 50 and calls it "equally
accurate" as the saddle point, which the table above contradicts by a factor of 1,400 at 64; a threshold default
changes draws whenever a later version moves the threshold or the form, which bears on reproducibility across
versions; logistic weights are not in 0.9-34, so no compatibility is owed. Options then: (a) a default above
b*, stated in the help, with a setting forcing exact draws; (b) opt-in only; (c) neither. Provisional
recommendation: (a), as dec-B143's "criteria for when to use them" describes.

## Constraints

- No package code lands from this study; scratch builds only. A slice the maintainer rules on follows the
  process for its class ([Process by risk](README.md#process-by-risk)).
- Any route lands above a threshold at or above the largest shape in a recorded snapshot, baseline or pinned
  test, so those stay bitwise.
- Draws come only from the chain's RNG; the per-row density cache is scratch, recomputed from omega, never
  state.
- Nothing from BayesLogit's source (GPL 3; dbarts is GPL 2 or later).
- Out of scope: the probit binomial (dec-B567), the real-shape r step of negbin-real-dispersion, the cap of
  dec-B202, whether an outer sampler's own steps are valid.

## Steps

1. Part A steps 1 to 4 above; review of the scratch build and checks; Part A's verdict to the orchestrator.
2. Open calls 1 and 2 to the maintainer before Part A's verdict is read against the tolerance.
3. If Part A passes everywhere: the finding and route A's landing plan; the approximate branch closes.
4. Otherwise C1 and C2 for the shapes left, Reading 1; then D's pilot, sizing audit and main run, Reading 2.
5. Report or finding, referees if a default would change, Open call 3 if live.

## Verification

- The design's probes (scratch/apx-design/probes, scratch/apx-crit) reproduce the tables above to their Monte
  Carlo error; Part B of the earlier design is folded into Part A step 1 and the candidates table.
- Positive controls: Part A's always-Metropolis estimated-r arm and its stale-cache build fail the
  joint-distribution check; in C, the truncated sum without remainder fails Reading 1 at large b and the normal
  at b = 10 with n = 1. A harness passing a control has no power and stops the study.

## Where this plan departs from its critique

- Taken: the remainder's exact moments and the re-ranked candidates (rerun, confirmed); the truncated-sum
  proposal and the gamma proposal's unbounded weight; the partially collapsed failure; exact draws after
  swaps and reactivation; the measured density cost and the meaning of certified; the bar on the posterior
  scale with the ESS-400 reading dropped and the run-length direction corrected; C2's size; family-wise
  control; the exact sum dropped as a C1 control above 1e3; Geweke and SBC pointed at route A's code with
  SBC's exact arm first; Part D's choice of datasets stated, the schools outcome defined, the bias-ratio gate
  dropped; Part A as the gate; the tolerance as the maintainer's value call; Open call 3 deferred.
- Not taken: that route A fails on estimated-r nbinom "even when r does not move" (argument and toy above;
  Part A's joint-distribution check decides on the engine). Added: Metropolis on the remainder alone, tried as
  a cheaper route and rejected on its weight; outside changes other than weights (offsets, responses,
  predictors) treated as swaps are.
- Not rerun: the density timing (the machine was at load 35); the critique's figure, taken at load 19 to 25,
  probably overstates it.

## Not verified

- Route A's sweep cost: one f evaluation is measured under heavy load; q's in the same pass and the cache are
  not.
- BayesLogit's switching points and its manual's silence, as read by the critique.
- How often estimated r moves per sweep, which sets estimated-r nbinom's cost under route A.
- The density's node counts and accuracy at real b below 16.
- SBC's and Geweke's detection limits: computed from their designs, checked in their pilots.
- van Dyk and Jiao (2015) read in its abstract only; the argument for the moved-then-exact rule is this plan's,
  checked on a toy, not taken from the paper.
- No published exact sampler for real shape below 1 or published error bound for the normal or saddle point
  was found (WPS, PSW, polyagamma's and pgdraw's documentation, D'Angelo and Canale 2022).
