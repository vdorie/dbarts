# approximate-draws: may a Polya-Gamma latent whose exact cost grows with its shape be drawn another way?

Status: DESIGN, 2026-11-03, before its blind critique. No study fit has run; the probes below are latent-level
measurements and one toy posterior. The arc was reopened by dec-B567; admitting any approximate draw stays the
maintainer's call.

agent: opus for the prototype, the harnesses' statistics and the engine scratch builds; one opus review of the
harnesses before Part C runs and of the reading before it reaches the maintainer; two independent referees on
the report if a default would change ([Studies](README.md#studies)).
rng: NEUTRAL for the study itself: benchmark harnesses and scratch builds that never land. Whatever the study
recommends lands as its own slice: shifting for an exact route, posterior-changing for an approximation, both
above a threshold that keeps every recorded draw bitwise.
window: after 1.0-0 unless the maintainer pulls it in; dec-B567 opened it on 2026-11-03.
budget: prototype about 600 lines of C (scratch), harnesses about 1,200 lines of R and 300 of C; compute about
25 core-hours expected, 40 on the doubled path, stop at 50 ([Compute](#compute)).

## Summary for the maintainer

The question: when a row's hidden Polya-Gamma value is expensive to draw exactly because the row carries many
trials, a large frequency weight or a large count, can it be drawn in time that does not grow with that
number, and if the only fast way is approximate, is its effect on what users read small enough to state and
accept?

Why it matters now. The exact draw is a sum of one small draw per trial, about 160 nanoseconds each. A logistic
fit of 200 rows spends 2 microseconds a row per sweep on everything else; at 100 trials a row the sweep is
9 times slower, at 10,000 trials 800 times. A binomial response (your ruling of today) makes large trial counts
the ordinary case: county disease rates over thousands of births, school results over hundreds of pupils, vote
counts in the hundreds of thousands. Negative binomial fits with counts in the thousands pay the same way, and
the negative binomial study already left one real dataset out because a replication took 1.5 hours. Real-valued
frequency weights and a real-valued negative binomial shape wait on a draw for a shape that is not a whole
number.

What was found while designing (measured here, not yet the study):

- Approximations differ by orders of magnitude. Against the exact distribution, computed by numerical
  inversion, a normal with the right mean and variance is off by 3 percent of probability mass at 64 trials and
  0.25 percent at 10,000; a gamma with the right mean and variance, at the same cost, is about 5 times closer;
  the saddle-point sampler of Windle, Polson and Scott is off by 0.002 percent at 64 and far less above; the
  truncated sum of gammas that the record had named gets worse as the trials grow (30 percent at 10,000), so it
  is no candidate for large shapes, only for shapes below 2.
- On a one-parameter toy where the exact posterior is known, the normal made a single row's interval 10 percent
  too wide at 10 trials and under 1 percent at 64 or more; the gamma 2 percent at 10. Pooled over 20 rows
  neither was distinguishable from exact.
- There is a third route the record set aside: keep the chain exact by drawing from the cheap gamma and
  accepting or rejecting it against the true density, a Metropolis step like the ones the trees already take.
  It accepts about 99 percent of proposals at 64 trials and 99.7 percent at 1,000, its cost does not grow with the
  trials, and it needs neither the saddle-point envelope nor the unproven conjecture the exact construction
  rests on. Its one hard part, evaluating the true density to certified accuracy, is the part the exact
  construction already designed. Its draw is not independent of the last one, so whether it meets "every draw
  exact" is your call (Open call 1).

The study, in four parts, each run only if the last leaves a question:

- A. Build and time the exact routes in a standalone prototype: does a bounded exact draw cost a few
  microseconds, as estimated? If so, an approximation could make a sweep at most about 2.5 times faster than
  that, which bounds what approximating could buy.
- B. The latent-level table above, completed over shapes from 16 up, to place each candidate's threshold.
- C. Checks against exact posteriors: the one-parameter toy over trials, success rates down to 1 in 1,000 and
  pooling; single-tree BART fits whose posterior is computed exactly by enumeration; a joint-distribution check
  and simulation-based calibration for errors those miss. A candidate is admissible from the smallest number
  of trials at which its shift of every posterior mean and quantile stays under 2.5 percent of a posterior
  standard deviation and its interval widths within 2 percent: half the Monte Carlo error of a run at the
  effective sample size Vehtari and colleagues recommend, so it adds at most about a tenth to the error a user's
  run already has.
- D. On real covariates with made-up outcomes and the real trial counts (North Carolina SIDS by county, California
  schools, the adult census collapsed into weighted cells, Titanic, and migration counts for negative binomial),
  at equal running time: coverage and width of 90 and 95 percent intervals for each unit's probability or rate,
  group rates and a contrast an analyst would report, with RMSE and bias, against the exact sampler.

If a candidate passes C and D, the report states its threshold, its measured effect and its speed gain over the
best exact route, and you decide whether to admit it. If a default would change, the result comes as a short
research report under your standard for such studies, with two referees first. About 25 core-hours.

## Goal

A measured answer, per kind of shape (large whole shapes; real shapes at or above 2; real shapes below 2), to
which draw dbarts should use: today's exact sum, a bounded-cost exact route, or an approximation above a stated
threshold; with the criterion fixed before any fit and the evidence the maintainer needs to rule.

## Context

Rulings, quoted, not reopened:

- dec-B14: "Every draw in 1.0 stays exact: an approximate Polya-Gamma draw ... is declined." Real dispersion
  and real binary weights wait on that unlock.
- dec-B143: the large-shape cost is documented for 1.0-0; the maintainer, 2026-09-29: "I don't know
  necessarily that we'll want to build it, we may instead want to permit approximations under some reasonable
  criteria for when to use them."
- dec-B567: the binomial response joins logistic, one draw of shape n per row; the maintainer, 2026-11-03:
  "there's a much higher likelihood that users will want to use it with large numbers of trials." The arc is
  opened on that ground; "whether any approximate draw is admitted remains the maintainer's."
- dec-B563: frequency weights reach the Polya-Gamma families as the shape times w, "a sweep's cost grows with
  it"; dec-B568: real weights where the draw is exact, the Polya-Gamma families whole numbers "until a
  real-shape draw exists"; dec-B566 and dec-B570 fix what a weighted or zero-weight row draws and means, which
  no route here changes.
- dec-B15: integer dispersion, "fully exact"; dec-B202: trial counts and fixed shapes capped at one million
  because a run that cannot finish is not supported ("if something just isn't going to run due to wall time
  constraints it's hard to argue that should support it"). A bounded draw removes that reason; revisiting the
  cap is a consequence, not in scope.

Engine: the shape seam is [`simulatePolyaGammaShape`](../../src/bartcore/model.hpp), a sum of round(b) draws
of [`ext_rng_simulatePolyaGamma`](../../src/external/random.c) (Devroye's method), used by
[`NBResponse::drawOmega`](../../src/bartcore/model.hpp) and the bridge's
[`drawAugmentation`](../../src/R_interface_bartcore.cpp). Logistic and hazard.logistic weights
([`LogisticResponse::refreshLatents`](../../src/bartcore/model.hpp)) and multinomial trials
([`MultinomialForestCombiner::drawForestGlue`](../../src/bartcore/combiner.hpp)) sum in their own loops. An
inactive row's draw is skipped, not discarded, and the draws poll for an interrupt. The binomial row of
dec-B567 is a logistic row of shape n times w (docs/design/frequency-weights.md on branch frequency-weights,
its binomial section).

Background: [nbinom-large-shape-pg.md](nbinom-large-shape-pg.md) designs the exact bounded-cost
construction, records that the saddle-point sampler is approximate, measures BayesLogit's sampler failing at
large shape, and derives the bound f_b / sp_b <= C_b used below. Its option (c), a Metropolis step, was set
aside because "the chain stays exact but the draw is not, and the bridge's replay has no previous omega".

### Where large or real shape arises, and today's cost

Measured 2026-11-03 on the laptop (arm64, load 5 to 9, so per-round minima of 2 or 3 rounds), shipped build of
main; standalone C timing of the engine's primitive linked against that build.

| source | shape per row | example users meet | cost today |
|---|---|---|---|
| binomial trials (dec-B567) | n w | NC SIDS: births per county 248 to 21,600; schools: 101 to 4,120 pupils; Florida 2000 votes per county 2,400 to 625,000 | 160 ns per trial |
| logistic, hazard.logistic frequency weights (dec-B563) | w | identical rows collapsed: Titanic's 24 cells to 670, a census collapsed over discrete covariates | the same |
| nbinom | (y + r) w, r estimated at most 50 or fixed | Migration (carData): counts to 99,400, 830,000 in all; any fixed r in the thousands | the same |
| multinomial trials | n_i per category | votes over several parties | 480 ns per trial over 3 categories |
| real shape | real w (power priors, fractional counts); real r (five of the nbinom study's twelve datasets have a shape below the grid's floor of 1) | | no exact draw below 1; refused today |

- Primitive: PG(1, psi) 146 to 170 ns (psi 0 to 8).
- Logistic fit, 200 rows, 75 trees, one chain: 2 us a row per sweep at w = 1; 18 at 100; 166 at 1,000;
  1,630 at 10,000. nbinom fixed r: 4 at 10, 160 at 1,000, 1,590 at 10,000; estimated r on counts averaging 224,
  36. Multinomial, 3 categories: 3.3 at 1 trial, 50 at 100, 480 at 1,000.
- So beyond about 10 trials a row the latent is the sweep. A default fit (4 chains, 1,000 sweeps each) of
  1,000 rows of 10,000 trials is about 1.8 core-hours; at 100,000 trials, 18.

## Candidates

Errors as the sources state them, then as measured here. "TV" is total variation distance from PG(b, psi)
computed against the exact density by Fourier inversion of the Laplace transform (probe laterr.R, grid to
about 1e-8), maximum over psi in {0, 2, 8}. Costs per draw on the laptop.

| candidate | source's account of its error | TV at b = 64 / 1e3 / 1e4 | cost per draw | shapes |
|---|---|---|---|---|
| today: sum of b Devroye draws | exact (Polson, Scott and Windle 2013) | 0 | 160 ns times b | whole b |
| exact construction: saddle-point envelope, accepted against the certified density | exact given the bound and WPS Conjecture 15 (nbinom-large-shape-pg.md) | 0 | estimated 2 to 4 us, not measured | b > 64; real b >= 2 |
| independence Metropolis: moment-matched gamma proposal, accepted against the certified density | exact stationary distribution (Metropolis-Hastings); constructed here | 0 (chain); acceptance 0.989 / 0.997 / -; 0.979 at 16 | gamma draw plus one density, estimated 1 to 4 us, not measured | b where the density is certified; real b |
| saddle-point sampler (Windle, Polson and Scott 2014, sec 5) | "generates approximate J*(n, z) random variates"; no bound stated; Daniels' approximation "improves as n increases" | 2.2e-5 / 3.5e-7 / 1.1e-8 (falls like b^-1.5) | about 9 Devroye draws by WPS Table 1 (10 to 12 times faster than the sum at n = 100), about 1.5 us | b >= 13 in BayesLogit's hybrid |
| moment-matched normal (WPS sec 6, BayesLogit above 170, Python polyagamma above 50) | "The central limit theorem suggests that this is a reasonable approximation"; polyagamma: "equally accurate" as the saddle point | 0.031 / 0.0078 / 0.0025 (falls like b^-0.5) | 35 ns | large b; can go negative |
| moment-matched gamma | not found in the literature; constructed here | 0.0077 / 0.0019 / 0.0006 | about 40 ns | any b; positive |
| truncated gamma sum, K = 200 (Polson, Scott and Windle 2013 eq. 1; BayesLogit) | "quite accurate for such small values of the first parameter ... may not be the case ... for arbitrary b"; TODO recorded "bias bounded ~1e-3" | 0.026 / 0.10 / 0.32 (grows with b); 1.3e-3 to 3.9e-3 at b = 0.5 | 8 to 13 us | real b, small |
| truncated gamma sum plus its remainder as one moment-matched gamma | constructed here | 2.6e-5 / 1.0e-4 / 3.2e-4; 1e-6 to 4e-6 at b = 0.5 | 8 to 13 us | real b below 2 |
| alternate sampler (WPS sec 4) | exact for h >= 1 "if Conjecture 8 holds", verified numerically on [1, 4]; acceptance falls with h | 0 | linear in b | real b >= 1 as a remainder |
| auxiliary mixture (Fruhwirth-Schnatter et al. 2009, cited by PSW 2013 for large counts) | "a fast approximation, at the cost of introducing a more complex latent-variable structure" | - | - | not this seam: replaces the augmentation |

Findings that shape the design:

- The normal's claimed parity with the saddle point does not hold: at b = 64 it is 1,400 times further in TV.
  At equal cost the moment-matched gamma dominates the normal everywhere measured, so the normal is carried
  only as the published comparator.
- The truncated gamma sum is not a large-shape candidate; its mean bias is a fixed fraction (about 2 / (pi^2 K)
  at psi = 0), which grows in standard-deviation units like sqrt(b).
- The saddle-point sampler is within about 2x of the estimated cost of both exact routes, so it buys nothing
  over them if Part A confirms the estimate; it stays in only if Part A finds the exact routes slow.
- Perturbation bounds do not settle the question. Mitrophanov's bound (Alquier, Friel, Everitt and Boland
  2016, Theorem 2.1) limits the chain's error by the kernel's TV times a mixing constant, but a sweep redraws
  n latents, so the kernel's TV is up to n times a row's: at 1,000 rows and 1e-3 a row the bound is vacuous.
  The effect must be measured on the posterior.
- Per draw, error falls with b for every large-shape candidate, so admissibility is hardest at the smallest
  shape it would be used for; that is where the study spends its fits, and it is also where the exact sum is
  cheap enough to be the comparator.

Toy pilot (probe toy.R): one logit parameter shared by n rows of b trials at p = 0.2, prior N(0, 9), 2e5
Gibbs iterations against the posterior by quadrature; Monte Carlo error about 0.003 posterior sd.

| n, b | normal: mean shift, sd ratio | gamma: mean shift, sd ratio |
|---|---|---|
| 1, 10 | -0.048, 1.105 | -0.010, 1.018 |
| 1, 64 | +0.001, 1.007 | -0.001, 1.002 |
| 1, 1000 | -0.000, 1.002 | +0.003, 1.000 |
| 20, 10 / 64 / 1000 | all within 0.005 and 0.2 percent | all within 0.003 and 0.3 percent |

A shared parameter pools its rows' latents, so the worst case is a quantity informed by one row; BART can
isolate a row when its trials dominate the data, so the study includes isolated rows.

### What a user would see

- A threshold default: above b* a row draws the admitted form, the help states it with its measured effect,
  and a control setting forces the exact path. This is how the default samplers of BayesLogit and the Python polyagamma
  package behave, and what "criteria for when to use them" (dec-B143) describes.
- An opt-in: exact everywhere unless the user asks; the cost stays linear for anyone who does not read the help.
- An exact route: nothing to state; the help's cost sentence (dec-B143) and possibly the cap (dec-B202) change.

Help wording if a form is admitted, for example: "Above a shape of 64 (trials times weight, or count plus
shape), the Polya-Gamma draw is approximated by a gamma with the exact mean and variance; checks against exact
posteriors found it moves posterior means and interval ends by less than 2.5 percent of a posterior standard
deviation. control(exact.latents = TRUE) draws exactly at a cost linear in the shape." The control's name is
decided with the slice.

## How the effect is measured

On the posterior of what users read, against a known truth, never only on the latent:

- Exact posteriors (score against the truth itself): the one-parameter toy by quadrature, and single-tree BART
  by enumeration of the tree space with 1-D quadrature of the leaf marginals, as
  [logistic-reference.R](../../benchmarks/R/logistic-reference.R) part 1 does for 0/1 rows; here its rows carry
  b trials. Scored: posterior mean, sd and 2.5, 50, 97.5 percent quantiles of the logit and of the probability
  (the probability transformed per draw, averaged last).
- Joint-distribution check (Geweke 2004), as [geweke-mc.R](../../benchmarks/R/geweke-mc.R) runs it for the
  gaussian sampler: the successive-conditional chain is exact from step 0 only if every conditional is, so an
  approximate latent shows as disagreement between the arms. Structural, cheap, no convergence claim.
- Simulation-based calibration (Talts, Betancourt, Simpson, Vehtari and Gelman 2018) through
  [sbc.R](../../benchmarks/R/sbc.R)'s logistic machinery with binomial rows. Its resolution is coarse: at 300
  replications the mean normalized rank has standard error 1/sqrt(12 x 300) = 0.017, so it resolves shifts of
  about 0.12 posterior sd; it guards against gross errors (rare events, small b), and does not decide.
- User level at equal running time (VD's rule for sampler changes): coverage and width of 90 and 95 percent
  intervals, pointwise and simultaneous, RMSE and bias, against a semi-synthetic truth; paired with the exact
  arm on the same data. Equal time means each arm is given the exact arm's wall time, capped at 8 times the
  default length, beyond which the faster arm's gain is Monte Carlo error a user would not see.

## The reading that decides (fixed before any fit)

Reading 1, admissibility of candidate A at shape b. In every Part C exact-posterior cell at b, for each scored
quantity, the deviation from the exact posterior in units of the exact posterior sd has its 99 percent interval
inside +-0.025 (means and quantiles), and the sd ratio's inside [0.98, 1.02]. A cell whose interval straddles a
bar is run again at up to 4 times its length, and still straddling counts as failing. A detection by the Geweke
or SBC arm at b (p below 0.01 after Holm's correction over that arm's functionals) fails A at b. A's threshold b*
is the smallest scanned b from which every larger scanned b passes.

Source and argument: Vehtari, Gelman, Simpson, Carpenter and Burkner (2021) recommend requiring a rank-normalized
effective sample size above 400, "typically sufficient to get a stable estimate of the Monte Carlo standard
error". At 400 the Monte Carlo error of a posterior mean is 0.05 posterior sd and of a posterior sd about 3.5
percent. An approximation error of half that adds at most sqrt(1 + 0.25) - 1, about 12 percent, to the total
error of a summary at the recommended precision, and less for the longer runs dbarts's default gives most
quantities. This is the approximate-MCMC literature's framing, error weighed against computation (Johndrow,
Mattingly, Mukherjee and Dunson 2015; Alquier et al. 2016), with the yardstick a user's own run supplies.

Reading 2, at the user level. In every Part D cell, the candidate's paired difference from the exact arm in
coverage (90 and 95, per unit, per group, per contrast) has its 3-standard-error band inside +-0.02, and its
RMSE and absolute bias ratios to the exact arm's are at most 1.05 at equal time. A candidate passing Reading 1
and failing here is not admitted at that shape, and the report says what Reading 1 missed.

Reading 3, speed: reported, not gated. In-situ per-sweep time against today's sum and against the bounded
exact route, same machine, alternating arms, per-round minima, at the shapes Part D's data carry.

Recommendation rule: recommend admitting A above b* where Readings 1 and 2 pass, with Reading 3's gain over the
best exact route stated beside it; recommend nothing approximate where an exact route is within 2x of A at the
sweep level, since then the approximation buys less than a second chain would. The maintainer rules.

## Designs and data

Part C toy cells (C1): candidates normal, gamma, gamma-sum-plus-remainder, saddle point (if Part A keeps it),
and the Metropolis route as a positive control; b in {16, 32, 64, 128, 256, 1e3, 1e4, 1e5}; success
probability 0.5, 0.05, 0.001 (the rare-event cells: at b = 1e4 and 0.001 a row holds about 10 events, and |psi|
is near 7, where the gamma's TV is largest); pooling n = 1 and n = 20 rows; real b in {0.3, 0.7, 1.5} for the
real-shape candidates. Chosen to span where each candidate's error and the posterior's shape differ: the toy is
the worst case of pooling and the only design with a truth exact to quadrature at every shape.

Part C single-tree cells (C2): logistic-reference.R's one-predictor, three-cut design with binomial rows,
b in {16, 64, 256, 1e3} per row, 2 rows per cell (near-isolated) and 25 (pooled); a nbinom twin with counts of
b. On scratch builds whose seam calls the candidate above a compile-time threshold.

Part C structural (C3, C4): Geweke on logistic with binomial rows at b in {16, 64}, n = 50; SBC at each candidate's
provisional b* from C1 and C2, n = 50, from the prior.

Part D datasets, chosen as the places large shapes occur in practice, each available offline (R, its
recommended packages, or installed CRAN packages; adult through the k study's loader):

| design | rows | shape | why | estimands |
|---|---|---|---|---|
| North Carolina SIDS, 1974 (sf's nc) | 100 counties | births 248 to 21,600, rate about 0.002 | disease mapping, the classic aggregated binomial; rare events | county rate; rate by region; rate ratio for the top against the bottom quartile of nonwhite births |
| California schools (survey's apipop), 1,000 sampled | 1,000 | enrollment 101 to 4,120 | administrative counts with rich covariates | school proportion; district means; contrast by school type |
| adult census, collapsed over its discrete covariates | thousands of cells | cell counts to the thousands | frequency weights as collapsed identical rows | cell probability; probability by education and by sex; sex contrast within occupation |
| Titanic (datasets) | 24 cells | weights to 670 | every row isolated, weights large | survival by cell; class contrasts |
| Migration (carData), nbinom | 90 | counts to 99,400 | large counts; left out of the nbinom study for cost | flow mean; destination totals; distance ratio |

Truths: on each dataset's real covariates and real trials, two surfaces fitted to the real outcome, a BART fit
and a smooth additive model, so the answer does not rest on BART's own prior (as binary-interval-study.md
did); outcomes drawn from them afresh per replication. Real outcomes are a diagnosis outside the rule. Arms:
the exact sum (and the bounded exact route if Part A builds it), each candidate passing Reading 1 at default
length and at equal time.

Real shape at the user level (the five nbinom datasets with a shape below 1) needs the real-shape r move of
TODO negbin-real-dispersion and is not in this study; C1's real-b cells carry the latent-level evidence it
would build on.

## Sizing

| component | question | size | cost | precision reached | smallest decision-changing difference | stop |
|---|---|---|---|---|---|---|
| A | exact routes' cost and acceptance | 1e6 draws per (b, psi), b 16 to 1e6 | 0.2 core-h | timing to 5 percent | a bounded route within 2x of the gamma at the sweep level | none |
| B | latent TV, thresholds | 8 b x 3 psi x 6 forms | 0.1 | TV to 1e-8 | b where TV crosses 1e-3 | none |
| C1 | posterior deviation, toy | 2e5 iterations per cell, about 240 cells | 0.5 | 99% half-width 0.008 sd (mean), 0.006 (sd ratio) | 0.025 sd, 2% | per cell: rerun at 4x only if straddling |
| C2 | posterior deviation, single tree | chains to ESS 2e4 per cell, set by a pilot of the autocorrelation | 2 | 99% half-width about 0.018 sd | the same bars | the same |
| C3 | joint distribution | 1e5 steps per arm, 4 runs | 4 | resolves about 0.05 sd at 3 se, pilot-checked | a detection | stop a run at its planned length |
| C4 | calibration | 300 replications, 100 draws | 6 | about 0.12 sd | a detection | none |
| D | user level, equal time | 10 replications per truth (pilot 3), to 20 if a band straddles | 12, 24 doubled | paired coverage band about +-0.01 estimated | 0.02 coverage | sequential at 10: stop when every band clears |

C1's size: ESS about 1e5 at an autocorrelation time near 2 (measured in the pilot), so the standard errors are
0.0032 sd and 0.0022 in the sd ratio. D's paired band is an estimate; the 3-replication pilot measures it and
the reviewer audits the table before the main run ([Compute](#compute)).

## Compute

About 25 core-hours expected, 40 on the doubled path, stop at 50. Parts A to C2 run on the laptop at 2 cores;
C3, C4 and D on an x86-64 Linux host. D per replication: exact arms about 1,200 core-seconds (Migration 530,
schools 400, SIDS 210, adult 70), candidates at default and equal time about 450 each. Parts run in order and
each stops the arc if it settles the question: A finding a bounded exact route within 2x of the gamma ends the
approximate branch for large whole shapes (Reading 3's rule), leaving C1's real-shape cells.

## Report

If a default would change, a PDF of at most 5 pages of text, readable alone, refereed twice before it reaches
the maintainer:

1. The problem: where large and real shapes arise, today's cost (the table above, re-measured on the host).
2. Methods: each candidate with its source's account of its error and its cost; the exact routes.
3. Evaluation: exact posteriors, the joint-distribution check, SBC, semi-synthetic data on real covariates and
   trials; the criterion and its source (Vehtari et al. 2021; Johndrow et al. 2015).
4. Results: latent TV against b; posterior deviation against b with the bars; coverage, width, RMSE and bias at
   equal time by dataset; speed.
5. Recommendation and limits: the threshold, the help's wording, what was not tested (real shape at the user
   level, the probit binomial, multinomial at the user level).

If no default changes, the finding goes in docs/design/approximate-draws.md in the plan-and-finding form.

## Decision

### Open call 1: does a Metropolis-corrected latent meet "every draw exact"?

Background: today each hidden value is an independent draw from its exact distribution. The Metropolis route
proposes a value from the cheap gamma and keeps it with a probability computed from the true density, or keeps
the previous value. The chain then targets exactly the same posterior; what changes is that about 1 row in 100
at 64 trials (1 in 300 at 1,000) keeps its last value for a sweep, as tree proposals already do when rejected.
dec-B14 names approximate draws, which this is not; its words "every draw exact" can be read either way.

- (a) It counts as exact: the arc builds this route first (about 600 lines of C with its tests, no conjecture
  to close, cost estimated 1 to 4 us a row); the approximate study then runs only for what it does not reach.
  The augmentation door the code running BART inside a larger sampler calls, which has no previous value, keeps
  the exact sum unless that caller passes the last latent; decided with the slice.
- (b) Only independent exact draws count: the exact construction (1,000 to 2,000 lines, resting on closing
  WPS Conjecture 15) or an admitted approximation; the study runs in full.

Recommendation: (a). The posterior is exact, which is what dec-B14 protects, and the package's tree moves are
Metropolis steps already.

### Open call 2: if a form passes, is it the default above its threshold or an opt-in?

Background: the default samplers of BayesLogit and the Python polyagamma package switch to an approximation
above a shape; an opt-in leaves anyone who does not read the help with the cost linear in the trials. Under a default
the report standard applies (a PDF, two referees) and the help states the effect; under an opt-in the finding
form suffices.

- (a) Default above b*, stated in the help, a control setting forcing the exact path. Report about 1 extra
  agent-day.
- (b) Opt-in only, exact by default.
- (c) Neither: exact routes only, whatever the evidence.

Recommendation: (a), as the criteria "for when to use them" of dec-B143 describe; under Open call 1 (a) the
recommendation rule above already withholds any approximation an exact route nearly matches.

## Constraints

- No package code lands from this study; scratch builds only, prefix apx-*. A slice the maintainer rules on
  follows the process for its class ([Process by risk](README.md#process-by-risk)).
- Any route lands above a threshold at or above the largest shape in a recorded snapshot, baseline or pinned
  test, so those stay bitwise; logistic and multinomial reach it only through the seam.
- Nothing from BayesLogit's source (GPL 3; dbarts is GPL 2 or later); reimplement from the papers.
- The exact gates keep exercising the exact path; an admitted form gets its own gate arm at its threshold.
- Out of scope: the probit binomial (truncated normal, dec-B567), the r move of negbin-real-dispersion, the
  cap of dec-B202.

## Steps

1. Part A: standalone C prototype in scratch: the certified density of PG(b, 0) by trapezoid inversion with its
   error bound (nbinom-large-shape-pg.md, its Design), the gamma proposal and the Metropolis step caching the
   current value's density per row; the exact construction's stage 1 only if call 1 is (b). Time it; check the
   chain's draws against the exact sum (two-sample KS, b <= 1e3) and its moments above.
2. Part B: extend probe laterr.R to the scan and record the table.
3. Part C1 in C; C2 on scratch builds; review of both harnesses; then C3 and C4.
4. Reading 1 per candidate; if none passes anywhere it is used, stop and report.
5. Part D: 3-replication pilot, sizing audit, main run with the sequential stop; Readings 2 and 3.
6. Report (or finding), referees, the maintainer.

## Verification

- The design's probes (latent TV by inversion, the toy, the fit costs) sit in the arc's scratch until Part B
  moves them to benchmarks/R; rerun there, they reproduce the tables above to their Monte Carlo error.
- Positive controls: the Metropolis route and the exact sum pass Reading 1 in every C cell; the truncated gamma
  sum (K = 200) fails it at large b, and the normal at b = 10 with n = 1 (the toy pilot's 10 percent width
  error): a harness that passes either has no power and stops the study.

## Not verified

- The exact routes' cost: both estimates (2 to 4 us, 1 to 4 us) rest on counting operations, not on a timing.
- The saddle-point sampler's cost here: taken from WPS Table 1 (a 2014 laptop), not run.
- SBC's and Geweke's detection limits: computed from their designs, to be checked in their pilots.
- Single-tree autocorrelation with binomial rows, hence C2's run lengths.
- Whether any later published work gives an exact sampler for real shape below 1 or a bound for the normal or
  saddle-point error: searched (WPS, PSW, polyagamma's and pgdraw's documentation, D'Angelo and Canale 2022),
  none found.
- The adult census collapsing and the schools' outcome surface: the cell structure and the real outcome used
  for each truth are fixed in the pilot.
