# k-move-pricing: pricing the collapsed scale move under logistic

Status: planned, in critique (2026-11-06); amended after its blind critique ([Amendments](#amendments)).
Nothing has run. The sequence is dec-B530's; every rule below is written before any fit of the step it decides.

agent: one opus implementer for step 0 (the move rebuilt on main), one opus review of it that runs the move's
mutants; one opus implementer for the harness and steps 1 to 3, run on the bench box; one review of the harness
and its smoke (the gate-only tier of [Process by risk](README.md#process-by-risk)); one opus review of the step 2
reading before it goes to the maintainer; two independent referees on the report when one is owed
([The report](#the-report)).
rng: the study changes no package code. Step 0 is NEUTRAL with the engine's toggle off and NEUTRAL with it on for
every fit outside the move's reach; where the move reaches (a logistic fit drawing k under a finite prior scale)
it is POSTERIOR-CHANGING by the README's rule once it lands on, the stationary distribution unchanged by
construction.
window: after the k study, before 1.0-0 (dec-B530). Engine slices stay serial.
budget: step 0, the move's own diff (about 300 engine lines, 950 of tests/cpp, 290 of tinytest, 255 of its
exact gate) replayed onto main less the probit half, plus the conflict work below and about 250 lines of
tests/cpp rewritten for the move's new placement; the pricing harness about 1,000 lines of R, most of it reused
from the k study's. Compute about 33 core-hours expected through step 2, about 39 with step 3, 52 at most with
every contingent arm, stop at 55 ([Sizing and stops](#sizing-and-stops)).

## Goal

After this plan runs, the logistic collapsed scale move has one of four outcomes, each reached by a rule
written here first: on for every logistic fit it can reach, with no trigger; on where a check on the data
fires; an opt-in; or parked. Each outcome rests on the move as it would ship on main's engine, and on
calibration of what users estimate at equal running time.

## Context

- The finding that opened pricing: [binary-k-mixing-need.md](../design/binary-k-mixing-need.md). Under logistic,
  where the signal is strong, the shipped sampler's k does not reach the small values the data want in a default
  run and intervals come out too narrow; the move removes that harm in every cell where it is shown (seven
  readings, all on probabilities and bands on them, none on an effect). Under probit no harm from k is shown.
  Its plan: [binary-k-mixing-need.md](binary-k-mixing-need.md), whose regimes, selection bands and per-fit
  times this plan reuses ([Regimes](binary-k-mixing-need.md#regimes),
  [Pre-registration](binary-k-mixing-need.md#pre-registration)).
- The sequence and its framing rulings: dec-B530 (the four steps, each rule written first); dec-B455 (a gain is
  less distance from nominal on either side, corrected for noise); the maintainer's rule of 2026-10-17 (a change
  to a hyperparameter's mixing is judged by what users estimate, calibration the target, over-coverage an error,
  at equal running time, on interesting real problems); dec-B555 ("k's" means getting k right closes the gap
  with today's tree moves, not that the trees are clear).
- What the move is: leaves times a factor r and k divided by it, r drawn by one slice step from its conditional
  with the latents integrated out (Liu and Sabatti 2000), taken on a coin of one in four, then the latents drawn
  at the scaled fit. Survey: [k-tail-samplers.md](../design/k-tail-samplers.md); logistic mixing:
  [k-mixing-pg-families.md](../design/k-mixing-pg-families.md). It is built on the branch collapsed-move-main,
  with its plan there as docs/plans/collapsed-scale-move.md.
- The probit step it sits beside: [`Chain::drawForestRescale`](../../src/bartcore/chain.hpp),
  [`Chain::forestRescaleApplies`](../../src/bartcore/chain.hpp), the test-only toggle
  [`probitRescaleForest`](../../src/bartcore/chain.hpp); the sweep's opening augmentation draw
  [`drawAugmentation`](../../src/bartcore/chain.hpp) (dec-B586).
- The pilot that measured the move at equal time and fixed how equal time is spent:
  [binary-interval-pilot.md](../design/binary-interval-pilot.md). The September prior study:
  [binary-hyperprior.md](binary-hyperprior.md).
- What 0.9-34 did: it had no logistic family. Its probit fits drew k under the improper chi(1.25, Inf) with no
  scale move of any kind; 1.0-0 draws k under chi(1.5, 2) for both links and adds the probit step.

## The sequence

| step | question | decided by | pass | fail |
|---|---|---|---|---|
| 0 | the move on main's engine | its correctness gates, and no draw moved with it off or out of reach | pricing starts | stop; fix and regate |
| 1 | is the move's small-k tail right? | SBC ranks of k and three functionals, under the default prior and one tilted to small k; the small-k share; near-0-or-1 coverage | step 2 | pricing stops |
| 2 | does it win at equal time? | gain on the k study's problem cells; every class anywhere shown harmless or improving; cost | on for every logistic fit | step 3 or 4 |
| 3 | only if the cost is too high: a check on the data | a trigger fixed on step 2's simulated fits, tested on real datasets not yet seen | on where it fires | step 4 |
| 4 | an opt-in | the maintainer, on steps 1 to 3 | | |

## Step 0: the move on main's engine

### What moves onto main

collapsed-move-main carries the move's sixteen commits on main as of 83b6e2b99 (2026-10-29), 596 commits behind
main today. Its diff against that base is 2,628 lines added in 16 files: the plan (993), tests/cpp (948), the
engine in chain.hpp and model.hpp (299), tinytest (about 290 over six files), the exact gate script (255),
benchmarks/README.md and man/dbarts.Rd. It adds one virtual to the response model in model.hpp and nothing to
the bridge, the C API in dbarts.h or the stored state, so it is no ABI event; installs take `--preclean`.

A three-way merge of main and the branch, run read-only, conflicts in five files: chain.hpp (two hunks), the
exact gate script (three), benchmarks/README.md, test-logistic-weight-swap.R and the plans index (one each). The
textual conflicts are small; three changes on main since the base are not, and the rebuild answers each:

1. **The probit step now scales the forests that share one prior total** (forests-shared-total, its calls
   recorded as dec-A251): [`Chain::forestRescaleApplies`](../../src/bartcore/chain.hpp) admits a combining
   chain's sharers, and [`Chain::drawForestRescale`](../../src/bartcore/chain.hpp) scales their summed fit
   against the offset plus every other forest's contribution. The branch's shared predicate and leaf-scaling
   helper were written for one forest. The helper is rewritten over the sharers, and the shared conditions
   split from the family test as the branch did, so the two moves cannot disagree on them. Whether the
   logistic move also reaches sharers is call 2.
2. **A sweep draws its augmentation variables first** (dec-B586), after the per-sweep hook and before the level
   step, the probit step and the trees. The move reads the fit, the response, the offset, the counts and the
   mask, never a latent. So it goes after the hook and before [`drawAugmentation`](../../src/bartcore/chain.hpp),
   and drops its own latent redraw: the sweep's own draw, at the scaled fit, is that redraw, and every recorded
   channel is still written after it. That saves one Polya-Gamma pass on every sweep that takes the move. Left
   as built, the move would follow the sweep's draw and redraw the latents again, paying the pass twice. Exactness
   is the same argument as on the branch (a kernel on the marginal of trees, leaves and k, then the latents from
   their conditional before anything reads them), so the exact gate and start-from-truth decide it, not a new
   derivation. A canceled run can then stop between the move and the draw. The next sweep draws the latents
   before anything reads them, so no stale latent is read; review checks that no reader outside a sweep takes
   them in that window.
3. **Leaf scales per tree** (blocks() groups, dec-B506; the trees' sd on every forest, dec-B550). The move's
   density reads no leaf scale. The group acts on (mu, k) and leaves each standardized leaf invariant whatever
   its scale, as the probit step's comment states for sharers. Review confirms that every leaf scale main
   carries is a multiple of 1/k.

The probit half (the branch's commit adding the probit evaluator and its tests) is not replayed (call 1). The
branch's three test and help commits that came with it (count-weighted logistic expectations, hurdle seeds, the
weights help) are replayed only where main does not already carry them. Two tests/cpp checks of the branch,
`checkCollapsedMapping` and `checkLogisticRedraw`, test the move's own latent redraw, which the rebuild drops; they
become sweep-level tests (the move, then [`drawAugmentation`](../../src/bartcore/chain.hpp)), about 250 lines
rewritten. The plan file is brought up to date with main: its status, its recheck against main, and its open calls
answered by this plan's calls 1, 2, 5 and 9.

The rebuild is a fresh branch off main replaying the commits, not a merge of main into the old branch. It stays a
branch until pricing decides (call 6). Pricing builds two private libraries from one commit of it: "move" as
committed, and "shipped" with the engine toggle `collapsedRescaleForest` flipped off by a one-line patch kept
beside the harness, as the k study did. Step 2's burn-in-only arm is a third library, from a patch that takes the
move only on burn-in sweeps.

### Gates before pricing starts

Pricing spends compute only on a move that is correct and touches nothing else. On the rebuilt tip, against its
own libraries, independently of the implementer:

- **Nothing moves with it off.** The shipped library draws bit for bit what main draws: the current equivalence
  baseline (equivalence-a6e83f38, 55 scenarios) 55 of 55 identical, bcf-equivalence and multinomial-equivalence
  identical, the four reference snapshot files passing, and the k study's 22 paths (both links; numeric, factor
  and missing predictors; test rows; 4 and 10 chains; an offset; fixed k; the sampler loop) bitwise against main.
  The probit step's loop over sharers moves into the shared helper, and none of the 55 scenarios fits sharers, so
  forests-blocks-equivalence's scenarios also run on fixed seeds, the shipped library against main, every fit and k
  bitwise, with the tests/cpp sharer cases of the probit step.
- **Nothing outside its reach moves with it on.** The move library also reproduces all 55 scenarios bitwise: the
  two logistic scenarios draw k under chi(1.5, Inf), where the move declines before the coin, and every other
  scenario is another family or a fixed k. bcf-equivalence (its latent_logistic scenarios are BCF, which the move
  does not reach) and multinomial-equivalence (its two-category case is the multinomial family) are bitwise too.
  With call 2 taken, the forests-blocks equivalence is bitwise too
  unless one of its scenarios is a logistic fit with sharers and a drawn k at a finite scale, which is then
  named as an expected mover.
- **It moves where it should.** It differs from the shipped library on every logistic path with a drawn k at a
  finite scale, and a fixed-k fit is bitwise equal on both libraries. A longer run's first draws equal a shorter
  run's on both.
- **It is right.** The move's tests/cpp block plain and under the address and undefined-behavior sanitizers; the
  exact gate's five logistic arms (pure, mixed, offset, mask, counts) quick and full, at |z| 4.5, with the
  logistic mixing bound set on both machines (the move's plan, early measurement 3); start-from-truth under
  logistic at |z| 3, R = 2,000; the logistic census bar (tau of k at most 45 median and 175 at the 99th
  percentile, the fit's worst at most 90); with call 2 taken, an exact gate arm of two forests sharing one
  total under logistic; [outer-step-exact.R](../../benchmarks/R/outer-step-exact.R) and
  [mask-redraw-exact.R](../../benchmarks/R/mask-redraw-exact.R) in full on both libraries, the gates of where a
  sweep draws its latents against an outer step that reads them, whose logistic arms (one tree, k drawn under
  chi(1.5, 2)) the move reaches: they are expected to move on the move library and must pass. A failure stops
  pricing until fixed and regated.
- The cheap gates of the engine tier: lintr, air, doc-freshness, the mutation battery's anchors.

Gates run on the old base do not carry over; every gate above reruns on the rebuilt tip.

### Gates at landing

After pricing decides, the move lands in the form the decision takes, through the engine tier's full battery for
a posterior-changing slice ([RNG classes and their gates](README.md#rng-classes-and-their-gates)): the review
with the move plan's 21 mutants (less the probit ones; its 9 and 15 replaced by two for the new placement: the
move after `drawAugmentation` with the latents left as drawn at the pre-move fit, caught by the exact gate and
outer-step-exact; the move before the per-sweep hook, caught by outer-step-exact's logistic-conditional arm),
the full tinytest suite, R-loaded ASan over the new and touched test files, the leak sanitizer on the x86 box,
`R CMD check --as-cran`, every exact gate in quick mode and the move's in full, the logistic-k SBC arm joining
the matrix with its confirmation run, the `logistick` equivalence scenario and the re-recorded baseline with its
oracle, the consumers' suites. These are the landing's, not pricing's. Before the merge, a few fits of each
priced regime run on the rebased tip and must be bitwise equal to pricing's fits on the same seeds; any that
differ rerun step 2's first look on that regime.

## Common to steps 1 to 3

- **Logistic only.** Every arm fits `family = "logistic"` under the package defaults (75 trees, chi(1.5, 2), 500
  burn-in, 500 draws a chain; 4 chains, or 10 for the effect models as bartCause runs them), `n.threads = 1`, on
  the bench box.
- **Regimes**: the k study's, on its real covariates, with its selection bands fixed in its pre-registration
  (moderate separation score 0.574 to 0.807, strong above 0.807 and below 0.853, near-separated 0.853 to 0.915,
  0.995 or more rejected as deterministic). Truths are drawn from the model's own prior (SBC), so an exact
  sampler covers at nominal on average over replications, and selection reads the data alone
  ([`sbcReplication`](../../benchmarks/R/sbc.R)'s route, as the k study built it).

  | id | regime | covariates, n | selection | step 1 | step 2 set |
  |---|---|---|---|---|---|
  | S1 | near separation | tictactoe, banknote, wdbc, ionosphere; 500 or all rows | near-separated band | yes | A |
  | S2 | rare outcome, small n | pima, german, creditapproval; 500 | rate in [0.01, 0.05] or [0.95, 0.99], score at most moderate | yes | A |
  | S3 | rare outcome, large n | adult, 10,000 subsampled | as S2 | no | A |
  | S4 | small n, strong signal | kyphosis, birthwt, haberman, cleveland, infert; 81 to 200 | strong band | yes | A |
  | S5 | many predictors | ACIC 2016's 58 covariates and 192 permuted copies; 1,000 | none | no | B |
  | S6 | moderate control | pima, german; 500 to 768 | moderate band, rate in [0.2, 0.8] | no | B |
  | S7 | effects, moderate | IHDP, 747, its real treatment | none | no | B |
  | S8 | effects, rare outcome | ACIC 2016, 1,000, setting 31's treatment model | as S2 | no | B |
  | E3 | the code that runs the BART sampler inside a larger sampler: 30 percent of outcomes imputed every sweep | pima, german; 500 | moderate band | no | B |

  A is the regimes with strong signal or rare outcomes, where the k study's seven problems are; B is fits with
  no strong signal. The probit-only regimes (S9, D2, E1, E2), the multiple-treatments designs (D1) and the real
  outcomes of B are not priced ([Left out](#left-out)).
- **Seen and held out.** All 22 real prediction datasets of the k study have been seen, with their real outcomes
  (tictactoe, banknote, mushroom, spambase, ionosphere, sonar, wdbc, climate, biopsy, magic, bank, adult,
  creditapproval, cleveland, pbc, kyphosis, pima, german, transfusion, infert, birthwt, haberman), and so have
  lindner, ACTG 175 and lalonde, and IHDP, ACIC 2016 and 2019 and Twins as designs. Steps 1 and 2 use some of
  them for covariates only, under drawn truths. Real outcomes enter only through two sets of datasets not yet
  looked at, R2 (step 2's held-out score guard, 6 datasets) and R3 (step 3's test, 12), drawn from these
  candidates, UCI binary classification sets with no time-to-event status: AI4I 2020 predictive maintenance,
  breast cancer Coimbra, breast cancer (Ljubljana), cervical cancer risk factors, chronic kidney disease,
  congressional voting records, default of credit card clients, diabetic retinopathy Debrecen, early-stage
  diabetes risk, hepatitis, HTRU2, ILPD, mammographic mass, myocardial infarction complications (lethal outcome),
  online shoppers purchasing intention, phishing websites, Polish companies bankruptcy (third year), QSAR
  biodegradation, raisin, red wine quality (quality 7 or more), rice (Cammeo and Osmancik), SPECTF heart,
  Taiwanese bankruptcy prediction, thoracic surgery, vertebral column (two classes), wilt. They are chosen to
  span outcome rates from about 0.03 to 0.6, sizes from about 100 to 30,000 rows and the fields people fit
  binary outcomes in: health, credit and commerce, and the physical sciences. Struck before any check, from the
  UCI descriptions: Statlog heart (largely Cleveland's patients, and cleveland is seen); EEG eye state, occupancy
  detection and ozone (time series); musk (conformations of the same molecules); parkinsons (repeated recordings
  of the same people); seismic-bumps (sequential shifts). Before any fit each candidate is checked: rows
  independent (no time series, no repeated subjects or instances), not drawn from the same source as any seen
  dataset, at least 100 rows, a recorded sha256 for the file; training rows are capped at 4,000 as the September
  study did. These facts are from the UCI listings as known and are confirmed at this check, before the split;
  if fewer than 24 pass, the list is extended by the same rule until 24 do. Then each separation score is
  computed, from the data alone, and kept sealed in the harness: it is used for the split and not written out
  until step 3 starts, since it is the trigger's own feature. Sorted by score ascending, in consecutive blocks of
  three, one dataset of each block goes to R2 and two to R3, which position goes to R2 set by a seed recorded
  beside the list, until R2 holds 6 and R3 12; the rest are reserves, taken from the same block if a dataset
  fails. The assignment is committed before step 2 runs.
- **Equal running time.** Per regime, the move's time ratio r (move over shipped) on the box, idle: a full
  default fit as the user runs it, three replications a regime, three alternating repetitions per library,
  per-round minima, the load stamped; r is the median of the nine pairs. The shipped arm at equal time takes
  ceil(1000 r) - 500 draws a chain after the same 500 burn-in, so its total sweeps equal the move's total time at
  its per-sweep cost, the pilot's frozen choice of spending the extra time on draws
  ([What was compared](../design/binary-interval-pilot.md#what-was-compared)). Its default reading is the first 500
  draws of the same run, by the prefix property. The burn-in-only arm's equal-time comparator is the prefix of
  that run at ceil(1000 r_b) - 500 draws, r_b its own ratio. The sampler loop (E3) is timed by total fit time.
  These timings are also the cost benchmark of step 2, together with time-only cells at 10,000 and 50,000 rows,
  75 and 200 trees, both signals, from the move plan's large-cell generator, timed as the move plan's early
  measurement 1 did: both libraries restored from one stored state before every timed segment of 100 sweeps,
  ten interleaved repetitions in alternating order. Person-period cells of hazard.logistic and a count-weighted
  cell are timed the same way (call 2).
- **Measures.** For each class, the 90 percent equal-tailed interval's coverage of the truth, computed on the
  latent scale. Coverage is pooled over replications before the distance from 0.90 is taken, the distance on
  either side corrected for noise as dec-B455 defines it, and over-coverage counts as error. Every coverage
  reading is made two ways: pooled over rows, and per fit (each fit's own coverage, averaged). Any weighting by
  the data keeps an exact sampler on nominal, so both are valid targets; a gain counts only if it holds under
  both, and a harm counts if it shows under either. Standard errors are bootstrapped over replications, 4,000
  resamples, every cell of a replication resampled together. Beside coverage, for the report: root-mean-square
  error and bias of the posterior mean against the truth, the interval score (Gneiting and Raftery 2007), and
  width.
- **Classes**: interior probability and probability near 0 or 1 (labels from the shipped arm's equal-time reading
  pooled with the move's, a function of the data and the chains and never of the truth, the k study's rule); the
  band on probabilities; the log risk ratio (S2, S3, S7, S8); the average effect overall and on the treated,
  individual effects, subgroups and the log odds ratio (S7, S8). A class a regime does not have is not read there.
- **No told fit.** No rule here compares with a fit told the true k: pricing does not attribute a gap to k, it
  asks whether the move closes it, against nominal. So none runs.

## Step 1: the move's small-k tail is right

**Question.** Does the move sample k's posterior correctly where k is small, the region the shipped sampler
rarely reached and the move now reaches (on the k study's long runs a mean of 7.7 percent of the move's draws
fell under 0.05, and most of them on some real data)? A sampler can only be shown wrong here by a truth that lands
in that region often. Under the default prior only 5 of 1,040 of the k study's true k fell under 0.05, so a test
under that prior alone could hardly fail.

**Design.** Simulation-based calibration (Talts et al. 2018) of the move library, one chain a replication, 99
draws at a thin T after a burn-in B, on the covariates and sizes of S1, S2 and S4, under two priors used both
for the truth and for the fit:
- the default chi(1.5, 2), with the k study's selection for each regime;
- chi(1.5, 0.1), which puts 21.7 percent of true k under 0.05. Here the separation bands are dropped, since a
  correctness test must reach the near-deterministic data where small k lives; the deterministic rejection
  goes too; S2 keeps its outcome-rate test. Any rule that reads the data alone keeps SBC exact.
- contingent, chi(1.5, 0.5) (16 percent of true k under 0.2, 53 percent under 0.5): run only if the pilot finds
  over half of the tilted arm's replications at a separation score of 0.995 or more, where the likelihood is
  nearly flat along the scale and the move mostly samples the prior.

T and B come from a pilot of 20 replications a regime and prior, run long, which also reports each arm's
distribution of separation scores. T is the smallest multiple of 50 at which every pilot replication keeps at
least 50 effective draws of the 99 on each of the four ranked functionals below, a stuck average fit being as
likely as a stuck k near separation (dec-B463). B is at least nine times the slowest tau of those functionals in
the pilot. If T would pass 500, that arm is reported as a mixing finding and is not read as a failure. R = 600 a
prior, 200 a regime.

**Power, shown before the main run.** The README's detection floor was measured on the probit-k arm (150 rows,
3 predictors, 50 trees, thin 700), and its error is in a Gibbs draw of k, not in this move. So the tilted arm
runs first at R = 150 on two mutant libraries, the move plan's mutants 1 (the density in r against dr instead
of dr / r) and 3 (the k term dropped or underflowed); both must flag under rule 1. If either does not, step 1 is
restated as unable to see that class of error before its main run. About 1 core-hour.

**Rule, fixed now.** Step 1 passes when all of these hold:
1. Ranks. For each prior, four functionals are ranked: k, the mean latent fit over training rows, and the
   latent fit at two held-out rows fixed by index in advance. None of the eight (two priors times four, twelve
   with the contingent arm) may leave its ecdf band (Sailynoja et al. 2022) or show a chi-square p (20 bins over
   the 100 ranks) below alpha = 0.05 / 16 over the 16 rank tests (0.05 / 24 with the contingent arm), the SBC
   harness's two rules for k. At R = 600 the harness flags an error of one half in the shape of k's conditional
   about 80 percent of the time on its own arm
   ([R/sbc.R - simulation-based calibration](../../benchmarks/README.md#rsbcr---simulation-based-calibration));
   the mutant runs above say what this arm sees.
2. The small-k share. The statistic is the mean over replications of P(k < c | y) - 1{true k < c}, within 3
   bootstrap standard errors of 0, at cutoffs 0.2 and 0.5 under the default prior (2.0 and 8.0 percent of true
   k under them), 0.02 and 0.05 under chi(1.5, 0.1) (5.7 and 21.7 percent), and 0.2 and 0.5 under the contingent
   arm.
3. Coverage near 0 or 1, at the thinned length. For an exact sampler the interval between the 5th and 95th of
   99 draws covers in expectation 0.898 when the draws are independent, about 0.8945 at an effective-draw ratio
   of 0.5 and 0.885 at 0.25 (an AR(1) computation); the target is that figure at the arm's measured ratio. Rows
   near 0 or 1, labeled from the fit's own draws, cover within 3.29 standard errors of that target under each
   prior, pooled over regimes. Labels taken from the same draws share their Monte Carlo noise with the interval,
   a small bias the reading states.
4. Numerics: every draw finite in every fit; log risk and odds ratios computed on the log scale and finite; the
   largest latent magnitude reported.

A failure of rules 1 to 3 is confirmed by rerunning the failing prior on fresh seeds named before the rerun, and
fails only if the same functional, or the same rule, fails again. A non-finite draw fails at once.

Default-length coverage near 0 or 1 is reported here and decided in step 2.

**Stop.** One early look at R = 150 a prior: a p below alpha / 10 on any functional goes straight to the
confirming rerun; a non-finite draw ends step 1 as failed. Otherwise it runs to R = 600.

**After a pass**: step 2. **After a fail**: pricing stops (dec-B530) and the move is not turned on in any form.
A failure of an exact move is a defect in the code or its numerics, so the finding goes to the maintainer with
its diagnosis. Whether a fix reopens pricing is the maintainer's call, with the recommendation that a fix be
tested by rerunning step 1 on fresh seeds, since the step tests the code and not the idea.

## Step 2: it wins at equal time

**Question.** At the same running time, are intervals better calibrated with the move than without it, where
the k study showed harm, and no worse anywhere, on fits with strong signal or without? This also measures what
the move costs.

**Arms**, every replication: the shipped library at equal time (its default reading being the prefix); the move
library at the default length; the burn-in-only library at the default length. In E3 the burn-in-only arm is read
as the loop runs it: where the loop takes no burn-in through the engine, that arm is the shipped sampler, which is
the point call 4 turns on.

**Rule, fixed now.** With d(arm) the noise-corrected distance from 0.90 and the bounds one-sided at 0.0025 a look
(two looks; z 2.81):
1. **Gain.** The move's d is smaller than the shipped arm's at equal time by at least 0.02, with its lower bound
   above 0, under both weightings, on at least one of the k study's five non-band problem cells: S4 interior and
   near 0 or 1, S1 interior and near 0 or 1, S2 near 0 or 1. A fix must remove at equal time at least one harm
   the evidence showed. The bar is the house's noticeable coverage difference
   ([Bars and where they come from](binary-k-mixing-need.md#bars-and-where-they-come-from)).
2. **No harm.** For every regime of A and B and every class, under each weighting, the worsening d(move) minus
   d(shipped at equal time) is *harmless* when its upper bound at one-sided 0.025 a look is below 0.02 for the
   row classes, individual effects and subgroups (the house's noticeable difference, the k study's one bar), and
   below 0.035 for bands. *Harm* is shown when the worsening is at least that bar with its lower bound above 0 at
   z 2.81. One-value estimands (average effects, ratios) are read for shown harm only, as dec-B533 reads them on
   calibration. The ratios of root-mean-square error and interval score, move over shipped, are harmless when
   their upper bound is below 1.02. R2's held-out log score (20 splits a dataset, 80 to 20) is harmless when the
   upper bound of its worsening is below 0.01 nats a row, September's noticeable log score
   ([The default for now](binary-hyperprior.md#the-default-for-now)). The move is on only when every class is
   harmless or improves. Any class unresolved at full counts goes to the maintainer with the counts. This is the
   burden the binary interval study put on its guards: a class passes on a bound shown, never on a harm not seen
   ([Decision rule](binary-interval-study.md#decision-rule)).
3. **Cost.** The cost is too high when no-harm fails only in B while the gain holds, since at equal time the
   shipped sampler's extra draws then buy more than the move does on fits that do not need it. Or when any timed
   logistic cell, the large cells included, runs above 1.25 times the shipped library: the move plan's stop 4,
   set before any measurement (call 3).

**After the reading.**
- Gain, every class harmless or improving, cost not too high: **the move is on for every logistic fit it
  reaches**, with no trigger. The report follows ([The report](#the-report)), then the landing.
- Gain, harm shown only in B, or the cost bar: the cost is too high. If the burn-in-only arm meets rules 1 and 2 on
  its own, it goes to the maintainer as an alternative before step 3 (call 4); otherwise, step 3.
- Gain with harm shown in A, or no gain: neither a default nor a trigger is supported. Step 3 is skipped, since a
  trigger can only fire where the move gains, and step 4 goes to the maintainer with the recommendation to park
  the move.
- Gain, no harm shown, but a class unresolved at full counts: the maintainer decides, with the counts each
  unresolved class would need ([Sizing and stops](#sizing-and-stops)). Bands in B are expected here. Step 3 is
  fitted meanwhile if the cost bar fails.

**Read beside the rule, deciding nothing**: the burn-in-only arm on every class (its gain against the move's
says how much of the harm is the start, which the k study left open); effective draws per second of the row
estimands under both kernels, on 10 replications each of S1, S2, S4 and S6, both libraries, 4 chains of 5,000
draws after 1,000 burn-in; bias and width by class; the shipped default against the move, for comparison with the
k study.

**Expectation, written before any fit.** The k study's stored fits, read after its results and on the old
engine, have the shipped sampler at equal time covering per fit 0.779 against the move's 0.904 near separation,
0.839 against 0.914 for rare outcomes at 500 rows and 0.854 against 0.887 in small samples. So the gain half is
expected to pass, and the step will be decided by no-harm and cost.

**Stop.** Looks at the counts in [Sizing and stops](#sizing-and-stops). A regime stops at its first look once
every class has resolved: harmless, improving, or a harm shown. The gain rule reads at each look and stops
counting once met. A harm shown ends step 2 as failing no-harm, but the regimes reached so far keep their counts,
since step 3 is fitted on them. If a first look's standard error on any class exceeds 1.5 times its planned value,
power is restated, with the counts that would resolve it, before the regime goes on. No replication is added past
the full counts.

## Step 3: a check on the data, only if the cost is too high

**Question.** Can a rule that reads only the data find the fits that need the move, so that only they pay for it?
A rule on the data alone keeps the kernel exact. A fit that crosses it through `setData` changes kernel
mid-chain, and both kernels are exact.

**Form, fixed on the simulated fits of step 2 before any R3 dataset is read.** The candidates are the separation
score at or above t (t from 0.70 to 0.95 by 0.01), the smaller class's count at or below e (e in 10, 20, 40,
80), since the score misses rare outcomes (the k study's
[Could a rule find the fits that need it?](../design/binary-k-mixing-need.md#could-a-rule-find-the-fits-that-need-it)),
and the union of one of each. A replication's gain is d(shipped at equal time) minus d(move) on the row classes.
For each candidate the gain is pooled over the replications it fires on and, apart, over those it spares, so
noise in single replications cancels rather than being truncated into gain. The trigger is the candidate that
fires on the fewest of B's fits while the upper bound on pooled gain over its spared set is below 0.02 in every A
cell; the single-feature form takes a tie. The k study's 0.8 was chosen after its data and is not carried over.

**Test on R3.** For each of the 12 datasets: a long run of 8 chains of each library, 5,000 sweeps of burn-in and
20,000 recorded at thin 20, as the k study's references ran; then 20 default-length fits of each library.

A dataset needs the move when its long runs pass the k study's acceptance (both libraries' halves agree, k's
quantiles included, dec-B462) and then, scored against the pooled reference, 20 shipped fits at equal time hold
a reference mass on the row classes at least 0.02 from 0.90 (its t bound excluding 0), while 20 move fits at the
default length are within 0.02. It does not need the move when the shipped fits' upper bound is below 0.02. This
reads the default-length harm users meet, not long-run disagreement, which the trees' mixing drives as much as
k's (dec-B555). Datasets whose references fail acceptance, or that read neither way, are unresolved and counted
apart. Expected, from the k study: most high-score datasets fail acceptance, as 10 of its 12 did. The trigger is
computed on each dataset's training rows.

**Rule, fixed now.** The trigger passes when it fires on every R3 dataset that needs the move and spares at least
half of those that do not. Below half, a second kernel would buy back too little of the cost to be worth
carrying. With fewer than 3 resolved datasets on either side, sensitivity or specificity cannot be read and
step 3 reads unresolved, which goes to the maintainer with the counts. With 12 datasets this test is coarse: it
can show a trigger that misses a dataset needing the move, and little more, and the finding says so.

**After a pass**: the move is on where the trigger fires; the report follows. **After a fail**: step 4.

## Step 4: an opt-in

Only if steps 2 and 3 both fail, and step 1 passed. The opt-in's form (a logical `dbartsControl` argument fixed
at creation, default off, as the probit step's withdrawn switch was, about 120 lines; or a family setting) and
whether to offer it at all go to the maintainer with steps 1 to 3's evidence. The recommendation is written then,
against what they show.

## Sizing and stops

Per-fit seconds are the k study's for a shipped logistic default fit on the box (its cost per accepted
replication). A step 2 replication costs about 3.3 shipped fits (equal time, the move, burn-in-only at r near
1.1), plus 20 percent for scoring and start-up and 15 percent for shared memory bandwidth. Spreads are the k
study's per replication, of the paired difference in distance from 0.90; where a regime was not in the k study
the spread is assumed and marked. "Resolves" means harmless at a true difference of 0 with probability 0.9: the
upper bound at one-sided 0.025 below the bar, which needs a standard error of at most the bar over 3.24.

| component | question it serves | size (looks) | cost, core-hours | precision at full size | smallest difference that changes the call | minimum size that resolves it | sequential stop |
|---|---|---|---|---|---|---|---|
| step 0, correctness gates | is the rebuilt move right and inert off | as listed | about 5 (laptop and box) | gates | any gate failing | as listed | none: a failure stops |
| step 1 pilot | T, B and the score distribution | 20 a regime and prior, run long | about 1 | sets T | | | |
| step 1 power check | can the tilted arm see the move's errors | 150, two mutant libraries | about 1 | flags or not | either mutant unflagged | 150 | none |
| step 1 SBC | is the small-k tail right | 600 a prior, 99 draws, thin T (150, 600) | about 5 at T 200; 10 at most | the README's detection floor, checked by the mutants | half the shape of k's conditional | 600 | fail at 150 on p below alpha / 10, confirmed by a rerun |
| step 1, third arm (contingent) | the tail away from near-determinism | 600, as above | 2.5 at T 200; 5 at most | as above | as above | 600 | as above |
| step 2, S1 | gain, no harm | 200 (100, 200) | 0.65 | interior se 0.007 at a spread of 0.10; near 0 or 1 0.014 at 0.20 | 0.02 gain; 0.02 harm | 263 interior, 1,050 near 0 or 1, at true 0; 100 for the expected gain | first look, once resolved |
| step 2, S2 | gain, no harm | 480 (240, 480) | 1.75 | interior se 0.006 at a spread of 0.135 | as S1 | 480 interior | as S1 |
| step 2, S4 | gain, no harm | 200 (100, 200) | 0.3 | interior se 0.004 at a spread of 0.054 | as S1 | 77 interior | as S1 |
| step 2, S3 | no harm at large n | 100 (50, 100) | 4.8 | interior se 0.006 at an assumed spread of 0.06 | 0.02 harm | 95 interior; 213 to 3,600 near 0 or 1 | as S1 |
| step 2, S5 | no harm | 220 (110, 220) | 1.25 | interior se 0.003 at an assumed 0.05; near 0 or 1 0.006 at 0.09 | 0.02 harm; 0.035 on bands | 213 near 0 or 1; 770 to 1,050 bands | as S1 |
| step 2, S6 | no harm | 220 (110, 220) | 0.6 | interior se 0.003 at 0.043; near 0 or 1 0.006 at 0.09 | as S5 | as S5 | as S1 |
| step 2, S7 | no harm on effects | 200 (100, 200) | 3.5 | individual effects se 0.004 at an assumed 0.05; near 0 or 1 0.006 at 0.09 | as S5 | 66 effects; 213 near 0 or 1; bands as S5 | as S1 |
| step 2, S8 | no harm on effects | 200 (100, 200) | 4.0 | individual effects se 0.004 at 0.053 | as S5 | 77 effects; otherwise as S7 | as S1 |
| step 2, E3 | no harm in a loop | 220 (110, 220) | 1.7 | as S5, assumed | as S5 | as S5 | as S1 |
| step 2, R2 held-out score | no harm on real outcomes | 6 datasets, 20 splits | about 1 | log score se about 0.002 a dataset | 0.01 nats | 20 splits | none |
| step 2, timing | r, and the cost bar | 9 regimes, 8 large cells, hazard and counts cells | about 1 as 100-sweep segments (3.5 as full fits), box idle | ratio to about 0.01 | 1.25 on a cell | 9 pairs; 10 segments | none |
| step 2, effective draws a second | reported | 10 a regime, S1, S2, S4, S6, both libraries | 0.3 | | decides nothing | | none |
| step 3, R3 long runs | which datasets need the move | 12 datasets, 8 chains a library | about 5 | per dataset | | 3 resolved on each side | none |
| step 3, R3 default fits | the default-length harm on each | 20 shipped at equal time and 20 move fits a dataset | about 1 | t bound a dataset | 0.02 | 20 | none |
| step 4 | | | 0 | | | | |

**Audit, regime by regime, at a true difference of 0.** The bars are now 0.02 for every class but bands and
0.035 for bands, and a class passes only on a bound shown, so each count is set to resolve where it can.
- S1: the interior resolves 81 percent of the time at 200, near 0 or 1 29 percent. Both are cells where the
  shipped sampler's harm is known (3 to 19 points), and an improvement resolves them at the first look; at a true
  0 they would need 263 and about 1,050. Kept at 200.
- S2: raised from 200 to 480 (looks 240 and 480), which resolves the interior 90 percent of the time at its spread
  of 0.135; near 0 or 1 is a gain cell, read as S1's.
- S4: resolves at 200 (interior over 99 percent); near 0 or 1 is a gain cell.
- S3: the interior resolves 92 percent of the time at 100. Near 0 or 1, at spreads from 0.09 to 0.37, it would
  need 213 to 3,600, which is 5 to 170 more core-hours; at 100 it resolves only through improvement, the
  expected outcome for rare outcomes, and otherwise goes to the maintainer with its count.
- S5, S6, E3: raised from 100 (S6 from 200) to 220 (looks 110 and 220), which resolves the interior over 99
  percent of the time and near 0 or 1 91 percent at a spread of 0.09.
- S7, S8: raised from 100 to 200 (looks 100 and 200); individual effects and subgroups resolve over 99 percent of
  the time, near 0 or 1 88 percent at 0.09.
- Bands in B, at spreads of 0.3 to 0.35, resolve at 220 only 32 to 41 percent of the time; resolving them at
  90 percent needs 770 to 1,050 a regime, about 30 to 45 core-hours more over B. They are not sized to resolve:
  an unresolved band goes to the maintainer with its count, as step 2's reading says.
- Showing a true harm of exactly 0.02 at 90 percent power needs about 420 at a spread of 0.10 and 1,700 at 0.20;
  under this rule a true harm that is not shown still leaves the class unresolved, so the move does not turn on
  by it.
- The spreads for S3, S5, S7 and E3 are assumed. If a first look's standard error exceeds 1.5 times the planned
  value, power is restated before going on (step 2's stop).

**Totals.** About 33 core-hours expected through step 2: step 0 about 5, step 1 about 7, step 2 about 21, about 12
of it if the regimes stop at their first look. About 39 with step 3. 52 at most, with step 1 at its longest thin,
its third arm, and the timing run as full fits. That is within dec-B530's rough guess of tens of core-hours, mostly
step 2. Against the plan before its critique (25 expected, 30 with step 3, 38 at most), step 2's regimes grew from
12 to 19 core-hours under the shown-bound rule (S2 to 480, B's regimes to 200 or 220), and step 1 gained the
power check (1) and the contingent arm (up to 5), step 3 its default-length fits (1). The stop is 55, about 1.4
times the expected total with step 3. Under it, cuts go in this order: the contingent arm's thin to its pilot's
T, S3 to its first look, the other step 2 regimes to their first looks, R3 to its first 8 datasets in score
order. A smaller design, B at 100 a regime and S2 at 200 (about 7 core-hours less), would resolve B's near-0-or-1
classes at a true 0 only about 60 percent of the time and S2's interior 55 percent, so more classes would reach
the maintainer unresolved instead of deciding.

Where it runs: every timing and every fit of steps 1 to 3 on the bench box, at most four jobs at once, the box
idle for timing ([The x86 leg](../../benchmarks/README.md#the-x86-leg)); step 0's gates on the laptop, with its
x86 leg on the box.

## The report

If step 2 or step 3 turns the move on by default, the maintainer gets a research report before the ruling, to
the standard for changing a default ([Studies](README.md#studies)). The maintainer set that standard for default
priors and said mixing might need it too. A change of kernel for every logistic fit is a default users meet, so
the plan holds it to that standard. The report is a short PDF, at most 5 pages of text with graphs as needed,
self contained. It argues as the literature argues for a change of sampler. Exactness comes first: the
generalized Gibbs step of Liu and Sabatti (2000), the exact gate, SBC (Talts et al. 2018). Then accuracy and
calibration of the estimands users report at equal computation, and effective draws per second (Vehtari et al.
2021). Its estimands are natural to each regime: probabilities for prediction data, risk ratios at rare rates,
average, individual and subgroup effects where there is a treatment. It scores root-mean-square error and bias
beside coverage and width, and heterogeneous effects through the individual and subgroup classes. It says how
the regimes and datasets were chosen, and what 0.9-34 did. Two independent referees read it on a neutral brief
before it reaches the maintainer. If the verdict is park or opt-in, the finding is a design note in the usual
form.

## The proposals in TODO

Each proposal recorded under TODO k-move-pricing, with what this plan does with it.

- **Step 1, a prior tilted toward small k.** Adopted, amended. chi(1.5, 0.1) in truth and fit (21.7 percent of
  true k under 0.05), with the separation bands and the deterministic rejection dropped, since the test must
  reach where small k lives. The default prior is kept as a second arm, because dec-B530 names it, and a third,
  chi(1.5, 0.5), runs if the tilted arm turns out mostly near-deterministic. The thin is set by a pilot, because
  one core-hour is not enough at a thin that yields independent draws, and the arm's power is shown on two mutant
  libraries before it runs; about 7 core-hours. Rank uniformity, near-0-or-1 coverage and numerics are rules 1, 3
  and 4 of step 1.
- **The length of any told fit read at a small true k.** Adopted in a stronger form: no rule reads a told fit, so
  none runs. If a diagnosis later needs one under probit or under the tilted prior, its length is set first by
  the 10x and 40x check the k study ran.
- **A burn-in-only arm.** Adopted, as an arm of step 2 read by the same rule. It costs one default-length fit a
  replication and answers the k study's open question about the start. Its place in the sequence is call 4,
  because a schedule tied to burn-in never runs where the sampler is driven one sweep at a time with no burn-in,
  as the code that runs the BART sampler inside a larger sampler drives it.
- **Held-out scores as a gate, with the datasets named.** Adopted, amended. It is a no-harm guard in step 2 at
  0.01 nats a row, not a gain criterion: under logistic the held-out probabilities converge where k does not
  ([binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md)), so a gain there is not expected. It runs
  on R2, six datasets none of the probes has seen, named by rule above; about 1 core-hour.
- **What "wins" means.** Adopted: noise-corrected distance from nominal on either side, by class, pooled over
  rows and per fit, a gain counting only under both weightings and a harm under either.
- **The no-harm margin against nominal, not the told fit.** Adopted: every bar compares the two arms'
  distances from 0.90, and no told fit enters.
- **How equal time is measured.** Adopted, under Common to steps 1 to 3.
- **A number for "cost too high".** Adopted, as rule 3 of step 2, with the number put to the maintainer (call 3).
- **Whether the probit half stays out.** Answered: out (call 1).
- **The note from dec-B555.** Taken as framing. Step 2 asks for a gain only on the cells where the k study showed
  harm, all on probabilities, and only no harm on effects, since a gain on effects is not what that evidence
  predicts.
- **Whether the move's code size is a gate.** Declined as a gate, reported. The code is written (about 1,600
  lines with its tests, less the probit half), and the shape, consistency and what users want decide, not the
  work.

## Calls for the maintainer

1. **The probit half.** dec-B530 speaks of logistic fits and is silent on probit. Options: (a) out of pricing and
   of 1.0-0, probit keeping the step it ships (recommended); (b) priced the same way, adding probit arms to
   steps 1 and 2, about 15 to 20 core-hours, the cost of the move reaching every default bart() binary fit and every
   binary stan4bart fit, and the move's plan's exposure of a bitwise compare to the runner's erfc. The shape
   argument for (a): each link then carries one scale move for k, probit its exact rescaling step and logistic
   the collapsed one. The evidence agrees: no harm from k is shown under probit, and the move's largest probit
   gain, 3.4 points near separation (se 1.4), is not past the k study's rule. 0.9-34 had no scale move for probit
   and no logistic family.
2. **The move's reach under logistic.** Options: (a) the forests that share one prior total, as the probit step
   now reaches them, besides a single forest (recommended), about 40 to 60 engine lines with tests/cpp cases and
   an exact gate arm; (b) a single forest only, as built. (a) makes "every logistic fit" mean every fit with a
   drawn k at a finite scale, as the probit step's reach does; the density is unchanged, the offset taking the
   other forests' contributions, and the forests-shared-total plan's D5 already says the logistic move covers the
   group the same way. Pricing measures single forests; the sharers rest on the same density and an exact-gate
   arm. Second, hazard.logistic and count-weighted logistic fits, which the move reaches through the same
   likelihood: (c) in reach and timed, their scoring unpriced (recommended), since person-period data are large
   rare-outcome fits, S3's shape, which step 2 scores, and the timing reads them at the cost bar's worst rows;
   (d) scored as a regime of their own, person-period data under a drawn truth, about 4 core-hours; (e) out of
   reach, the move declining on them, at the cost of a kernel that differs between a logistic fit and the same
   fit expanded to person-period rows. 0.9-34 had neither logistic fits nor forests sharing a total.
3. **What "cost too high" is.** Options: (a) no harm failing at equal time only on fits with no strong signal, or
   any timed logistic cell above 1.25 times the shipped sampler, the move plan's own stop, set before any
   measurement (recommended); (b) equal time alone, with no wall-time bar, since equal time already charges the
   move its cost; (c) a ratio of 1.05, bench-sampler's regression gate, which the k study's timings (mostly 1.04
   to 1.14, on the old placement) make likely to send pricing to step 3. The 1.25 bar matters only for the large
   cells, which are timed and not scored, since equal time already charges the scored ones. 0.9-34 had no move
   to pay for.
4. **Burn-in-only, if the full move costs too much.** Options: (a) put it to you as the alternative before step
   3 when it meets step 2's rule on its own (recommended); (b) step 3 first, burn-in-only reported. Its cost:
   a fit run one sweep at a time from outside would never take the move, so a user's kernel would depend on how
   the run is split. The trigger has no such dependence. 0.9-34 had no move to schedule.
5. **A switch, if the move is on by default.** Options: (a) none, as dec-B397 ruled for the probit step
   (recommended unless step 2's cost bar was close), one kernel to test and document; (b) a logical
   `dbartsControl` argument, default on, about 120 lines. 0.9-34 had no switch for any move.
6. **Where the move waits.** Options: (a) on its rebuilt branch until pricing decides, landing once in the form
   decided (recommended, dec-B494's rule that a feature merges complete); (b) on main now with the engine toggle
   off and no R surface, so pricing reads main itself, at the cost of an engine change on main whose fate is
   undecided. Either way the bitwise check at landing applies. 0.9-34 shipped no engine toggle of this kind.
7. **The held-out datasets.** The candidate list and the rule that splits it into R2 and R3 are above; you may
   add, strike or reorder before step 2 runs, and nothing is fitted on them before then. No 0.9-34 study used
   held-out real datasets for a sampler change.
8. **The September prior study's revisit under logistic.** TODO binary-hyperprior-revisit waits for k's scale
   move to land (dec-B390), and its first condition is that k mixes. Options: (a) read that condition free from
   step 2's move fits at the default length, k's split-Rhat below 1.05 and its effective draws in the hundreds
   (recommended); if it holds and the move lands for logistic, a logistic rerun of the study goes to you as its
   own item after landing; (b) leave the revisit as it stands, a probit item, with nothing read here. The study
   was probit, and with the probit half out (call 1) its probit condition is not met by this plan either way.
   0.9-34 had no logistic family.
9. **Negative binomial's item.** The move's plan leaves negative binomial to its own item, and there is none.
   Options: (a) name TODO k-mixing-pg-families as that item, which already holds the Polya-Gamma families' k
   mixing (recommended); (b) open a new item. 0.9-34 had no negative binomial family.

## Left out

- **Probit**, every probit-only regime (S9, D2, E1, E2) and the probit half: call 1.
- **Negative binomial.** Its item is TODO k-mixing-pg-families (call 9). The negative binomial default study
  measured k under a scratch build of that half, and its default is still open.
- **The multiple-treatments designs.** Their miss is shared by both samplers and points to the trees and the model,
  and every design with a real effect failed its first long run.
- **The colleagues' case.** It has not been shared; it has a slot in the k study's plan (dec-B529).
- **Real outcomes scored by long-run agreement.** They carry no truth, and their 22 datasets are spent; real
  outcomes enter only as R2's held-out score and R3's trigger test.
- **Count-weighted logistic and hazard.logistic** are reached by the move, gated by its exact gate's counts arm and
  timed, but not scored as regimes: the likelihood is the same, and the counts were not a k study regime. Call 2
  puts their reach to the maintainer.
- **Multi-forest fits with amplitudes, BCF and the latent BCF families**: out of the move's reach, as of the probit
  step's.
- **The September prior study's rerun** (TODO binary-hyperprior-revisit, dec-B390): it was a probit study, and with
  the probit half out pricing does not supply the probit mixing fix its revisit condition asks for. It stays its own
  item; whether step 2 reads its first condition under logistic is call 8.
- **The trees.** The gaps both samplers share go to TODO tree-mixing-proposals whatever pricing finds.
- **Rare outcomes past 10,000 rows**: timed at 50,000, not scored.

## Verification

- Step 0: the gate list above, each with its count (55 of 55 identical twice, the 22 paths, the exact gate's |z|
  and bounds, the census figures, start-from-truth's |z|), on the rebuilt tip's own libraries, independently of the
  implementer.
- Steps 1 to 3: the harness's `check` action on its job lists (every replication of every arm present, no reading
  written to another's path, no label reading the truth), a smoke of every regime at one replication and a tenth
  of every length reviewed before step 1, and each step's reading made by the harness's `read` action, the rule
  applied as written here.
- `Rscript tools/check-doc-freshness.R .` on this file and the indexes.

## Not verified

- That the rebuild's three semantic changes are as small as read: the two conflicting hunks in chain.hpp and
  the placement before the augmentation draw were read on main, not built.
- The move's cost on main with the redraw dropped; the k study's 1.04 to 1.14 is the old placement's.
- That a thin of 500 or less gives 50 effective draws on each of step 1's four functionals under chi(1.5, 0.1),
  and the selection costs of the tilted arm's S2 rate test.
- The candidate datasets' sizes, rates and files: from the UCI listings as known, checked at the data stage.
- Per-replication spreads in step 2 between the move and the shipped arm at equal time; the sizes use the k study's
  spreads against its told fit.

## Amendments

The blind critique of the first draft made 22 points; each is listed by its number with what changed.

1. Taken. The shipped arm at equal time takes ceil(1000 r) - 500 draws after the 500 burn-in, so the move's whole
   cost is charged; the burn-in-only comparator likewise at r_b.
2. Taken. No harm is now a bound shown: harmless when the upper bound of the worsening is below 0.02 (0.035 on
   bands), harm when the worsening reaches the bar with its lower bound above 0, one-value estimands read for
   shown harm only, unresolved classes to the maintainer. The sizes were re-audited at a true difference of 0
   (S2 to 480, S5, S6 and E3 to 220, S7 and S8 to 200); bands in B and S3's near-0-or-1 class are not sized to
   resolve and say what would.
3. Taken. outer-step-exact and mask-redraw-exact run in full on both libraries; two placement mutants replace the
   move plan's 9 and 15.
4. Taken. forests-blocks-equivalence on fixed seeds and the tests/cpp sharer cases, bitwise with the move off.
5. Taken. Statlog heart and the six candidates without independent rows are struck, the list extended
   to 26, each candidate's facts checked at the data stage before the split, the list extended until 24 pass.
6. Taken. A dataset needs the move only when its long runs pass acceptance and the default-length harm shows
   against the reference; failed references are unresolved; fewer than 3 resolved on a side reads unresolved.
7. Taken. alpha is 0.05 / 16, and a failure of rules 1 to 3 is confirmed on fresh seeds before it stands.
8. Taken. The tilted arm's power is shown on mutants 1 and 3 before the main run, about 1 core-hour.
9. Taken. The thin keeps 50 effective draws on all four functionals; a thin past 500 is a mixing finding.
10. Taken. Cutoffs 0.2 and 0.5 under the default prior, 0.02 and 0.05 under the tilted one, read as a mean
    difference within 3 bootstrap standard errors.
11. Taken. The count of declines is dropped, since no build counts them; every draw finite stays.
12. Taken. The trigger is fitted on pooled gain over the fired and spared sets, not truncated per-fit gains.
13. Taken. The split is by score ascending in blocks of three, the position in each block by a recorded seed.
14. Taken. The table has the minimum size that resolves each component; the large cells are timed as
    restored-state 100-sweep segments (about 1 core-hour, 3.5 as full fits); a first look whose standard error
    passes 1.5 times planned restates power.
15. Taken. hazard.logistic and count-weighted reach is folded into call 2 as options (c) to (e); the September
    revisit's first condition is call 8; negative binomial's item is call 9.
16. Taken. `checkCollapsedMapping` and `checkLogisticRedraw` become sweep-level tests, about 250 lines budgeted.
17. Taken. The citation for the trees' sd on every forest is dec-B550.
18. Taken. The pilot reports the score distribution, and chi(1.5, 0.5) runs as a third arm, alpha to 0.05 / 24,
    when over half the tilted arm scores 0.995 or more: a correctness test should reach the tail where the
    likelihood still informs k, and the pilot costs nothing more to read.
19. Taken in part. The target is stated at the arm's measured effective-draw ratio. The label bias is noted, not
    removed: labeling from odd draws and scoring even ones would halve the draws behind each interval and move
    its target further than the bias does.
20. Taken. The separation scores stay sealed in the harness until step 3. Stratifying by rate and size instead
    was not taken, since the split by score is what spreads R3 over the range the trigger is tested on.
21. Taken. Call 2 cites D5 and says pricing measures single forests, the sharers resting on the same density and
    an exact-gate arm.
22. Taken. bcf-equivalence and multinomial-equivalence are listed as bitwise with the move on.
