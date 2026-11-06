# k-move-pricing: pricing the collapsed scale move under logistic

Status: planned, in critique (2026-11-06). Nothing has run. The sequence is dec-B530's; every rule below is
written before any fit of the step it decides.

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
exact gate) replayed onto main less the probit half, plus the conflict work below; the pricing harness about 900
lines of R, most of it reused from the k study's. Compute about 25 core-hours expected through step 2, about 30
with step 3, 38 at most, stop at 45 ([Sizing and stops](#sizing-and-stops)).

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
| 2 | does it win at equal time? | gain on the k study's problem cells; no harm on any class anywhere; cost | on for every logistic fit | step 3 or 4 |
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
3. **Leaf scales per tree** (blocks() groups, dec-B506; the trees' sd on every forest, dec-B562). The move's
   density reads no leaf scale. The group acts on (mu, k) and leaves each standardized leaf invariant whatever
   its scale, as the probit step's comment states for sharers. Review confirms that every leaf scale main
   carries is a multiple of 1/k.

The probit half (the branch's commit adding the probit evaluator and its tests) is not replayed (call 1). The
branch's three test and help commits that came with it (count-weighted logistic expectations, hurdle seeds, the
weights help) are replayed only where main does not already carry them. The plan file is brought up to date with
main: its status, its recheck against main, and its open calls answered by this plan's calls 1, 2 and 5.

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
- **Nothing outside its reach moves with it on.** The move library also reproduces all 55 scenarios bitwise: the
  two logistic scenarios draw k under chi(1.5, Inf), where the move declines before the coin, and every other
  scenario is another family or a fixed k. With call 2 taken, the forests-blocks equivalence is bitwise too
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
  total under logistic. A failure stops pricing until fixed and regated.
- The cheap gates of the engine tier: lintr, air, doc-freshness, the mutation battery's anchors.

Gates run on the old base do not carry over; every gate above reruns on the rebuilt tip.

### Gates at landing

After pricing decides, the move lands in the form the decision takes, through the engine tier's full battery for
a posterior-changing slice ([RNG classes and their gates](README.md#rng-classes-and-their-gates)): the review
with the move plan's 21 mutants (less the probit ones), the full tinytest suite, R-loaded ASan over the new and
touched test files, the leak sanitizer on the x86 box, `R CMD check --as-cran`, every exact gate in quick mode
and the move's in full, the logistic-k SBC arm joining the matrix with its confirmation run, the `logistick`
equivalence scenario and the re-recorded baseline with its oracle, the consumers' suites. These are the landing's,
not pricing's. Before the merge, a few fits of each priced regime run on the rebased tip and must be bitwise
equal to pricing's fits on the same seeds; any that differ rerun step 2's first look on that regime.

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
  looked at, R2 (step 2's held-out score guard) and R3 (step 3's test), drawn from these candidates, all UCI
  binary classification sets with no time-to-event status: breast cancer Coimbra, default of credit card
  clients, diabetic retinopathy Debrecen, EEG eye state, HTRU2, ILPD, mammographic mass, musk (version 1),
  occupancy detection, online shoppers purchasing intention, ozone (eight-hour), parkinsons, phishing websites,
  QSAR biodegradation, raisin, rice (Cammeo and Osmancik), seismic-bumps, SPECTF heart, Statlog heart, thoracic
  surgery, vertebral column (two classes), wilt. They are chosen to span outcome rates from about 0.05 to 0.6,
  sizes from about 100 to 30,000 rows and the fields people fit binary outcomes in: health, credit and
  commerce, and the physical sciences. Before any fit each candidate is checked: rows independent, at least 100
  rows, a recorded sha256 for the file, training rows capped at 4,000 as the September study did. Then its
  separation score is computed, from the data alone. Sorted by score, the eligible ones go alternately to R3
  and R2 until R2 holds 6 and R3 12, with the order within ties set by a seed recorded beside the list. The rest
  are reserves, taken in score order if a dataset fails. The assignment is committed before step 2 runs.
- **Equal running time.** Per regime, the move's time ratio r (move over shipped) on the box, idle: a full
  default fit as the user runs it, three replications a regime, three alternating repetitions per library,
  per-round minima, the load stamped; r is the median of the nine pairs. The shipped arm at equal time takes ceil(r
  times 500) draws a chain after the same 500 burn-in, the pilot's frozen choice of spending the extra time on
  draws ([What was compared](../design/binary-interval-pilot.md#what-was-compared)). Its default reading is the
  first 500 draws of the same run, by the prefix property. The burn-in-only arm's equal-time comparator is the
  prefix of that run at its own ratio. The sampler loop (E3) is timed by total fit time. These timings
  are also the cost benchmark of step 2, together with time-only cells at 10,000 and 50,000 rows, 75 and 200
  trees, both signals, from the move plan's large-cell generator.
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

T and B come from a pilot of 20 replications a regime and prior, run long. T is the smallest multiple of 50 at
which every pilot replication keeps at least 50 effective draws of the 99 on k, the rule the move's plan set
for its logistic-k SBC arm. B is at least nine times the slowest tau of k in the pilot. If T would pass 500,
T is 500, R falls to what 10 core-hours buy, and the power the README states for R = 600 is restated for that
R. R = 600 a prior, 200 a regime.

**Rule, fixed now.** Step 1 passes when all of these hold:
1. Ranks. For each prior, four functionals are ranked: k, the mean latent fit over training rows, and the
   latent fit at two held-out rows fixed by index in advance. None of the eight (two priors times four) may
   leave its ecdf band (Sailynoja et al. 2022) or show a chi-square p (20 bins over the 100 ranks) below
   alpha = 0.05 / 8, the SBC harness's two rules for k. At R = 600 this flags an error of one half in the shape
   of k's conditional about 80 percent of the time, and an error that scales with the leaf count at R = 50
   ([R/sbc.R - simulation-based calibration](../../benchmarks/README.md#rsbcr---simulation-based-calibration)).
2. The small-k share. Under each prior and at cutoffs 0.05 and 0.2, the mean posterior share of k under the
   cutoff equals the share of true k under it within 3 bootstrap standard errors.
3. Coverage near 0 or 1, at the thinned length: the interval between the 5th and 95th of 99 draws holds 0.90
   in expectation for an exact sampler with independent draws. Rows near 0 or 1, labeled from the fit's own
   draws, cover within 3.29 standard errors of 0.90 under each prior, pooled over regimes.
4. Numerics: every draw finite in every fit; log risk and odds ratios computed on the log scale and finite;
   the move declines on a non-finite log density only at a state of zero prior density, so the count of such
   declines outside those states is 0; the largest latent magnitude is reported.

Default-length coverage near 0 or 1 is reported here and decided in step 2.

**Stop.** One early look at R = 150 a prior: a p below alpha / 10 on any functional, or a non-finite draw, ends
step 1 as failed. Otherwise it runs to R = 600.

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
2. **No harm.** In every regime of A and B and every class, harm is shown when the lower bound of the worsening,
   d(move) minus d(shipped at equal time), passes its margin, under either weighting: 0.010 for the row
   classes, 0.015 for individual effects and subgroups, 0.040 for average effects and ratio estimands, 0.035
   for bands. The ratios of root-mean-square error and interval score, move over shipped, may not have a lower
   bound above 1.02. On R2, the move's held-out log score may not be worse than the shipped arm's at equal time
   with a lower bound past 0.01 nats a row on any dataset (20 splits a dataset, 80 to 20). The margins are the
   binary interval study's guards, fixed before its pilot ran, and September's noticeable log score
   ([The default for now](binary-hyperprior.md#the-default-for-now)).
3. **Cost.** The cost is too high when no-harm fails only in B while the gain holds, since at equal time the
   shipped sampler's extra draws then buy more than the move does on fits that do not need it. Or when any timed
   logistic cell, the large cells included, runs above 1.25 times the shipped library: the move plan's stop 4,
   set before any measurement (call 3).

**After the reading.**
- Gain, no harm, cost not too high: **the move is on for every logistic fit it reaches**, with no trigger. The
  report follows ([The report](#the-report)), then the landing.
- Gain, harm only in B, or the cost bar: the cost is too high. If the burn-in-only arm meets rules 1 and 2 on
  its own, it goes to the maintainer as an alternative before step 3 (call 4); otherwise, step 3.
- Gain with harm in A, or no gain: neither a default nor a trigger is supported. Step 3 is skipped, since a
  trigger can only fire where the move gains, and step 4 goes to the maintainer with the recommendation to park
  the move.

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
every class has resolved: a harm shown, or the upper bound of its worsening below its margin. The gain rule
reads at each look and stops counting once met. A harm shown ends step 2 as failing no-harm, but the regimes
reached so far keep their counts, since step 3 is fitted on them. No replication is added past the full counts.

## Step 3: a check on the data, only if the cost is too high

**Question.** Can a rule that reads only the data find the fits that need the move, so that only they pay for it?
A rule on the data alone keeps the kernel exact. A fit that crosses it through `setData` changes kernel
mid-chain, and both kernels are exact.

**Form, fixed on the simulated fits of step 2 before any R3 dataset is read.** The candidates are the separation
score at or above t (t from 0.70 to 0.95 by 0.01), the smaller class's count at or below e (e in 10, 20, 40,
80), since the score misses rare outcomes (the k study's
[Could a rule find the fits that need it?](../design/binary-k-mixing-need.md#could-a-rule-find-the-fits-that-need-it)),
and the union of one of each. A replication's gain is its per-fit decrease in distance from 0.90 on the row
classes, move against the shipped arm at equal time, and only positive gains count. The trigger is the
candidate that captures at least 90 percent of the summed gain over A and B while firing on the fewest of B's
fits, the single-feature form taking a tie. The k study's 0.8 was chosen after its data and is not carried over.

**Test on R3.** For each of the 12 datasets: a long run of 8 chains of each library, 5,000 sweeps of burn-in and
20,000 recorded at thin 20, as the k study's references ran. The dataset needs the move when the halves disagree
in substance: k's quantiles disagree, or 5 percent or more of row quantities do (the k study's category, which
needs no acceptance and so reads failed runs too). The trigger is computed on each dataset's training rows.

**Rule, fixed now.** The trigger passes when it fires on every R3 dataset that needs the move and spares at least
half of those that do not. Below half, a second kernel would buy back too little of the cost to be worth
carrying. If no R3 dataset needs the move, sensitivity cannot be tested and step 3 reads unresolved, which
goes to the maintainer with the counts. With 12 datasets this test is coarse: it can show a trigger that misses
a dataset needing the move, and little more, and the finding says so.

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
study's per replication.

| component | question it serves | size (looks) | cost, core-hours | precision at full size | smallest difference that changes the call | sequential stop |
|---|---|---|---|---|---|---|
| step 0, correctness gates | is the rebuilt move right and inert off | as listed | about 5 (laptop and box) | gates | any gate failing | none: a failure stops |
| step 1 pilot | T and B | 20 a regime and prior, run long | about 1 | sets T | | |
| step 1 SBC | is the small-k tail right | 600 a prior, 99 draws, thin T (150, 600) | about 5 at T 200; 10 at most | the README's detection floor | half the shape of k's conditional | fail at 150 on p below alpha / 10 |
| step 2, S1 | gain, no harm | 200 (100, 200) | 0.65 | interior se about 0.007 at a spread of 0.10 | 0.02 gain; harm past 0.010 | first look, once resolved |
| step 2, S2 | gain, no harm | 200 (100, 200) | 0.7 | as S1 | as S1 | as S1 |
| step 2, S4 | gain, no harm | 200 (100, 200) | 0.3 | as S1 | as S1 | as S1 |
| step 2, S3 | no harm at large n | 100 (50, 100) | 4.8 | interior se about 0.006 at an assumed spread of 0.06 | harm past 0.010 | as S1 |
| step 2, S6 | no harm | 200 (100, 200) | 0.55 | se about 0.003 at a spread of 0.04 | harm past 0.010 | as S1 |
| step 2, S5 | no harm | 100 (50, 100) | 0.6 | se about 0.005 at an assumed spread of 0.05 | harm past 0.010 | as S1 |
| step 2, S7 | no harm on effects | 100 (50, 100) | 1.8 | individual effects se about 0.005, assumed as S5 | harm past 0.015 | as S1 |
| step 2, S8 | no harm on effects | 100 (50, 100) | 2.0 | as S7 | as S7 | as S1 |
| step 2, E3 | no harm in a loop | 100 (50, 100) | 0.8 | as S5 | harm past 0.010 | as S1 |
| step 2, R2 held-out score | no harm on real outcomes | 6 datasets, 20 splits | about 1 | log score se about 0.002 a dataset | 0.01 nats | none |
| step 2, timing | r, and the cost bar | 9 regimes and 8 large cells, 9 pairs | about 1 (box idle) | ratio to about 0.01 | 1.25 on a cell | none |
| step 2, effective draws a second | reported | 10 a regime, S1, S2, S4, S6, both libraries | 0.3 | | decides nothing | none |
| step 3, R3 long runs and default fits | does the trigger find the fits that need it | 12 datasets | about 5 | per dataset | one dataset missed | none |
| step 4 | | | 0 | | | |

Totals: about 25 core-hours expected through step 2 (step 0 about 5, step 1 about 6, step 2 about 14, about 9 of
it if the regimes stop at their first look); about 30 with step 3; 38 at most. The stop is 45, about 1.5 times the
expected total with step 3. Under it, cuts go in this order: S3 to its first look, the other step 2 regimes to
their first looks, R3 to its first 8 datasets in score order.
That is within dec-B530's rough guess of tens of core-hours, mostly step 2.

How the sizes were audited. Rows of A are sized so that the gain, expected at 3 to 19 points per fit, resolves at
the first look, and so that the row classes' guard sees a harm of about 0.04 at 90 percent power. A harm of 0.02
would need about 1,700 replications a regime at these spreads, which no affordable size reaches, so in A the
guard catches large harm only, and the finding says so. Rows of B are where no-harm decides. At S6's spread, 200
replications show a harm of about 0.022 against the 0.010 margin; at S5, S7, S8 and E3, 100 replications show
about 0.03. The step 1 sizes follow the SBC harness's measured detection floor. More replications would refine
an estimate without changing a call, so none are planned.

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
  reach where small k lives. The default prior is kept as a second arm, because dec-B530 names it. The thin is
  set by a pilot, because one core-hour is not enough at a thin that yields independent draws; about 6 core-hours.
  Rank uniformity, near-0-or-1 coverage and numerics are rules 1, 3 and 4 of step 1.
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
- **The no-harm margin against nominal, not the told fit.** Adopted: every margin compares the two arms'
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
   other forests' contributions. The pricing regimes are single-forest fits, and sharing forests are the same
   model as blocks() of one forest by the forests-blocks equivalence, so pricing reads both. 0.9-34 had neither
   logistic fits nor forests sharing a total.
3. **What "cost too high" is.** Options: (a) no harm failing at equal time only on fits with no strong signal, or
   any timed logistic cell above 1.25 times the shipped sampler, the move plan's own stop, set before any
   measurement (recommended); (b) equal time alone, with no wall-time bar, since equal time already charges the
   move its cost; (c) a ratio of 1.05, bench-sampler's regression gate, which the k study's timings (mostly 1.04
   to 1.14, on the old placement) make likely to send pricing to step 3. 0.9-34 had no move to pay for.
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

## Left out

- **Probit**, every probit-only regime (S9, D2, E1, E2) and the probit half: call 1.
- **Negative binomial.** The move's plan leaves it to its own item. The negative binomial default study measured k
  under a scratch build of that half, and its default is still open.
- **The multiple-treatments designs.** Their miss is shared by both samplers and points to the trees and the model,
  and every design with a real effect failed its first long run.
- **The colleagues' case.** It has not been shared; it has a slot in the k study's plan (dec-B529).
- **Real outcomes scored by long-run agreement.** They carry no truth, and their 22 datasets are spent; real
  outcomes enter only as R2's held-out score and R3's trigger test.
- **Count-weighted logistic and hazard.logistic** are reached by the move and gated by its exact gate's counts arm,
  but not priced as regimes: the likelihood is the same, and the counts were not a k study regime.
- **Multi-forest fits with amplitudes, BCF and the latent BCF families**: out of the move's reach, as of the probit
  step's.
- **The September prior study's rerun** (TODO binary-hyperprior-revisit, dec-B390): it was a probit study, and with
  the probit half out pricing does not supply the probit mixing fix its revisit condition asks for. It stays its own
  item.
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
- That a thin of 500 or less gives 50 effective draws of k under chi(1.5, 0.1), and the selection costs of the
  tilted arm's S2 rate test.
- The candidate datasets' sizes, rates and files: from the UCI listings as known, checked at the data stage.
- Per-replication spreads in step 2 between the move and the shipped arm at equal time; the sizes use the k study's
  spreads against its told fit.
