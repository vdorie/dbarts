# k-move-pricing: pricing the collapsed scale move under logistic

Status: planned, in critique (2026-11-06); amended after its blind critique, and on 2026-11-07 for the
maintainer's rulings dec-B617 to dec-B626 and the convergence measure ([Amendments](#amendments)).
Nothing has run. The sequence is dec-B530's; every rule below is written before any fit of the step it reads.
The rules state what a reading shows and recommend; the maintainer decides from a research report (dec-B624).

agent: one opus implementer for step 0 (the move rebuilt on main), one opus review of it that runs the move's
mutants; one opus implementer for the harness, steps 1 to 3 and the convergence measure, run on the bench box;
one review of the harness and its smoke (the gate-only tier of [Process by risk](README.md#process-by-risk)),
which also checks the convergence measure's argument against the truth; one opus review of the step 2 reading
before its report is written; two independent referees on each report ([The report](#the-report)).
rng: the study changes no package code. Step 0 is NEUTRAL with the engine's toggle off and NEUTRAL with it on for
every fit outside the move's reach; where the move reaches (a logistic fit drawing k under a finite prior scale)
it is POSTERIOR-CHANGING by the README's rule once it lands on, the stationary distribution unchanged by
construction.
window: after the k study, before 1.0-0 (dec-B530). Engine slices stay serial.
budget: step 0, the move's own diff (about 300 engine lines, 950 of tests/cpp, 290 of tinytest, 255 of its
exact gate) replayed onto main less the probit half, plus the conflict work below and about 250 lines of
tests/cpp rewritten for the move's new placement; the pricing harness about 1,300 lines of R, most of it reused
from the k study's, the convergence measure's streams and scoring about 300 of them. Compute about 38
core-hours expected through step 2 with the convergence measure, about 44 with step 3, about 70 at most with
every contingent arm and every second look, stop at 66 ([Sizing and stops](#sizing-and-stops)).

## Goal

After this plan runs, the maintainer has chosen one of four outcomes for the logistic collapsed scale move,
from a research report on each step read (dec-B624), each outcome recommended by a rule written here first: on
for every logistic fit it can reach, with no trigger; on where a check on the data fires; an opt-in; or parked.
Each rests on the move as it would ship on main's engine, on calibration of what users estimate at equal
running time, and, read beside it, on how long each sampler takes to reach its stationary distribution from the
start a user gets (dec-B622).

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
- The maintainer's rulings on this plan, 2026-11-07: logistic only, probit keeping its own step (dec-B617); the
  move reaching the forests that share one prior total (dec-B618); hazard.logistic and count-weighted logistic
  fits in reach, timed and not scored (dec-B619); the two cost bars (dec-B620); step 2 as the move against the
  shipped sampler head to head at equal time, no arm switching samplers within a run (dec-B621); a measure of how
  long each takes to converge, taken from published comparisons (dec-B622); no switch if the move goes on by
  default, provisionally (dec-B623); bands reported with no bar, and every step's results reaching the maintainer
  as a research report, the maintainer deciding (dec-B624, dec-B625); the report's contents taken from published
  reports of its kind (dec-B626).
- How published comparisons measure the time to reach the stationary distribution:
  [sampler-convergence-measures.md](../design/sampler-convergence-measures.md).
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

| step | question | read by | what follows |
|---|---|---|---|
| 0 | the move on main's engine | its correctness gates, and no draw moved with it off or out of reach | pricing starts on a pass; on a fail, stop, fix and regate |
| 1 | is the move's small-k tail right? | SBC ranks of k and three functionals, under the default prior and one tilted to small k; the small-k share; near-0-or-1 coverage | step 2 on a pass; on a fail pricing stops (dec-B530) and the finding goes to the maintainer |
| 2 | does it win at equal time, and how long does each sampler take to converge? | gain on the k study's problem cells; every class anywhere shown harmless or improving; cost (dec-B620); beside them, the convergence measure (dec-B622) | a report to the maintainer whatever it shows, recommending the move on for every logistic fit, step 3, or parking; the maintainer decides (dec-B624) |
| 3 | a check on the data, only if the maintainer sends pricing here after step 2's report | a trigger fixed on step 2's simulated fits, tested on real datasets not yet seen | a report recommending the move on where the trigger fires, or step 4; the maintainer decides |
| 4 | an opt-in | the maintainer, on steps 1 to 3 | |

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
   split from the family test as the branch did, so the two moves cannot disagree on them. The logistic move
   reaches the sharers too (dec-B618).
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

The probit half (the branch's commit adding the probit evaluator and its tests) is not replayed (dec-B617). The
branch's three test and help commits that came with it (count-weighted logistic expectations, hurdle seeds, the
weights help) are replayed only where main does not already carry them. Two tests/cpp checks of the branch,
`checkCollapsedMapping` and `checkLogisticRedraw`, test the move's own latent redraw, which the rebuild drops; they
become sweep-level tests (the move, then [`drawAugmentation`](../../src/bartcore/chain.hpp)), about 250 lines
rewritten. The plan file is brought up to date with main: its status, its recheck against main, and its open calls
answered by dec-B617 (the probit half), dec-B618 (the sharers), dec-B623 (the switch) and this plan's call 9.

The rebuild is a fresh branch off main replaying the commits, not a merge of main into the old branch. It stays a
branch until pricing decides (call 6). Pricing builds two private libraries from one commit of it: "move" as
committed, and "shipped" with the engine toggle `collapsedRescaleForest` flipped off by a one-line patch kept
beside the harness, as the k study did. No other library is built: no arm switches samplers within a run
(dec-B621).

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
  With the sharers in reach (dec-B618), the forests-blocks equivalence is bitwise too
  unless one of its scenarios is a logistic fit with sharers and a drawn k at a finite scale, which is then
  named as an expected mover.
- **It moves where it should.** It differs from the shipped library on every logistic path with a drawn k at a
  finite scale, and a fixed-k fit is bitwise equal on both libraries. A longer run's first draws equal a shorter
  run's on both.
- **It is right.** The move's tests/cpp block plain and under the address and undefined-behavior sanitizers; the
  exact gate's five logistic arms (pure, mixed, offset, mask, counts) quick and full, at |z| 4.5, with the
  logistic mixing bound set on both machines (the move's plan, early measurement 3); start-from-truth under
  logistic at |z| 3, R = 2,000; the logistic census bar (tau of k at most 45 median and 175 at the 99th
  percentile, the fit's worst at most 90); with the sharers in reach, an exact gate arm of two forests sharing one
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
  draws of the same run, by the prefix property. The sampler loop (E3) is timed by total fit time.
  These timings are also the cost benchmark of step 2, together with time-only cells at 10,000 and 50,000 rows,
  75 and 200 trees, both signals, from the move plan's large-cell generator, timed as the move plan's early
  measurement 1 did: both libraries restored from one stored state before every timed segment of 100 sweeps,
  ten interleaved repetitions in alternating order. Person-period cells of hazard.logistic and a count-weighted
  cell are timed the same way, their timings entering the cost bar (dec-B619). The same timings give each
  library's seconds per sweep, which the convergence measure's readings in seconds use.
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

**Arms**, every replication, two mechanisms head to head (dec-B621): the shipped library at equal time (its
default reading being the prefix) and the move library at the default length. Both keep bart()'s fixed 500
burn-in, since it is the default; how long each takes to converge is read apart, in
[The convergence measure](#the-convergence-measure) (dec-B622).

**Readings, fixed now.** They state what step 2 shows and frame the recommendation; they decide nothing, the
maintainer deciding from the report (dec-B624). The bars below are the differences step 2's replications are
sized to resolve and the frame of the recommendation (dec-B625). With d(arm) the noise-corrected distance from
0.90 and the bounds one-sided at 0.0025 a look (two looks; z 2.81):
1. **Gain.** The move's d is smaller than the shipped arm's at equal time by at least 0.02, with its lower bound
   above 0, under both weightings, on at least one of the k study's five non-band problem cells: S4 interior and
   near 0 or 1, S1 interior and near 0 or 1, S2 near 0 or 1. A fix must remove at equal time at least one harm
   the evidence showed. The bar is the house's noticeable coverage difference
   ([Bars and where they come from](binary-k-mixing-need.md#bars-and-where-they-come-from)).
2. **No harm.** For every regime of A and B and every class, under each weighting, the worsening d(move) minus
   d(shipped at equal time) is *harmless* when its upper bound at one-sided 0.025 a look is below 0.02 for the
   row classes, individual effects and subgroups (the house's noticeable difference, the k study's one bar).
   *Harm* is shown when the worsening is at least that bar with its lower bound above 0 at z 2.81. One-value
   estimands (average effects, ratios) are read for shown harm only, as dec-B533 reads them on calibration.
   Bands carry no bar and no label (dec-B625): the report gives each sampler's empirical band coverage against
   nominal, with its interval and the replications behind it, and the difference between the samplers with its
   interval. The ratios of root-mean-square error and interval score, move
   over shipped, are harmless when their upper bound is below 1.02. R2's held-out log score (20 splits a dataset,
   80 to 20) is harmless when the upper bound of its worsening is below 0.01 nats a row, September's noticeable
   log score ([The default for now](binary-hyperprior.md#the-default-for-now)). The move is recommended on only when
   every class is harmless or improves. Any class unresolved at full counts is reported with the counts. This is the
   burden the binary interval study put on its guards: a class passes on a bound shown, never on a harm not seen
   ([Decision rule](binary-interval-study.md#decision-rule)).
3. **Cost.** The cost is too high when no-harm fails only in B while the gain holds, since at equal time the
   shipped sampler's extra draws then buy more than the move does on fits that do not need it. Or when any timed
   logistic cell, the large cells included, runs above 1.25 times the shipped library: the move plan's stop 4,
   set before any measurement (dec-B620).

**After the reading.** Whatever step 2 shows, its results go to the maintainer as a research report
([The report](#the-report)), which sets each sampler's empirical coverage against nominal, with widths,
root-mean-square error and bias, for every class and band, and the convergence measure beside them; the
maintainer decides (dec-B624). The readings frame the report's recommendation:
- Gain, every class harmless or improving, cost not too high: recommended, **the move on for every logistic fit
  it reaches**, with no trigger and no switch (dec-B623, provisional: the switch goes back to the maintainer
  beside step 2's reading and its times before the move lands on), then the landing.
- Gain, harm shown only in B, or the cost bar: the cost is too high, and step 3 is recommended.
- Gain with harm shown in A, or no gain: neither a default nor a trigger is supported, and parking the move is
  recommended (step 4). Step 3 would not be, since a trigger can only fire where the move gains.
- Gain, no harm shown, but a class unresolved at full counts: no recommendation past what is shown; the report
  gives the counts each unresolved class would need ([Sizing and stops](#sizing-and-stops)).

Step 3 runs only if the maintainer, having read this report, sends pricing there.

**Expected reading, written before any fit.** At the planned counts some classes resolve with probability 0.9 or
less even at a true difference of 0 (the audit in [Sizing and stops](#sizing-and-stops)). If the move improves
A's problem cells as the k study's stored fits predict, about 2 classes are expected unresolved, all in S3 and B,
and about 5 if B's near-0-or-1 spreads sit at 0.2 rather than the 0.09 assumed; the chance that every class
resolves is about 0.08 at best. If the move changes nothing anywhere, about 5 to 8 are expected unresolved. So
step 2 will most likely reach the maintainer with counts, whatever the move does, and the reading says which
classes and what each would need.

**Read beside the readings, framing nothing**: the convergence measure, which reads how much of the harm is the
start, the question the k study left open ([The convergence measure](#the-convergence-measure)); effective draws
per second of the row
estimands under both kernels, on 10 replications each of S1, S2, S4 and S6, both libraries, 4 chains of 5,000
draws after 1,000 burn-in; bias and width by class; the shipped default against the move, for comparison with the
k study.

**Expectation, written before any fit.** The k study's stored fits, read after its results and on the old
engine, have the shipped sampler at equal time covering per fit 0.779 against the move's 0.904 near separation,
0.839 against 0.914 for rare outcomes at 500 rows and 0.854 against 0.887 in small samples. So the gain half is
expected to pass, and the recommendation to turn on no-harm and cost.

**Stop.** Looks at the counts in [Sizing and stops](#sizing-and-stops). A regime stops at its first look once
every class has resolved: harmless, improving, or a harm shown. The gain rule reads at each look and stops
counting once met. A harm shown is reported and does not end step 2: the other regimes run to their counts, since
the report reads every regime (dec-B624). If a first look's standard error on any class exceeds 1.5 times its
planned value, power is restated, with the counts that would resolve it, before the regime goes on. No
replication is added past the full counts.

## The convergence measure

**Question.** From the start a user gets, how many sweeps, and how many seconds, does each sampler need before
its chains draw from the stationary distribution, the posterior? Step 2's arms keep bart()'s fixed 500 burn-in;
this measure runs each sampler with no burn-in and reads the whole stream from its start (dec-B622). It reads
how much of the harm is the start, which the k study left open, and whether 500 sweeps suit either sampler
under logistic. dec-B622 keeps the fixed burn-in as the default, so the plan proposes no change to it.

**What published comparisons do, and what is ours.** They read the distance between a chain's marginal at sweep t
and the posterior: as an upper bound from coupled chains, from many independent chains against a long-run
reference, or as agreement between chains from dispersed starts (R-hat read as a trajectory), each with a
threshold on its own distance ([sampler-convergence-measures.md](../design/sampler-convergence-measures.md)). No
coupling of a BART sampler exists, and building one would be new work with its own proof; references are missing
where the shipped sampler is slowest; and R-hat misses a slow direction that chains from one start share, as k's
drift from its starting value is shared. So the measure reads three ways. Reading 1, the truth's rank among
independent chains at each sweep, is the study's own: its argument is SBC's (Talts et al. 2018), but no published
comparison reads it sweep by sweep that the survey found. Reading 2 reads many chains against long-run
references, the published method, with a statistic of the study's own. Reading 3 is published as it stands.

**Regimes.** Step 2's problem cells and its fits with no strong signal: S1, S2 and S4, where the five non-band
problem cells are, and S5, S6, S7, S8 and E3 (set B). S3 is in neither and is left out: at 10,000 rows its streams
would cost about 6 core-hours at 2,000 sweeps and about 15 with the shipped sampler's longer streams.

**Starts.**
- *Default*, the start bart() gives: every tree a single leaf with a fit of zero, k at its starting value and
  the latents at their deterministic cold start. Every chain of a replication starts there and its stream differs
  from the others' only by its random numbers. This is the reading users meet. A probe of the shipped sampler on
  three of the k study's stored near-separated and small-sample datasets had k at a median of 2.3 to 2.5 after the
  first sweep, against posterior medians of 0.12 to 0.78.
- *Dispersed*: each chain from its own draw from the model's prior, trees, leaves and k, drawn by the route step
  1 and the k study draw truths by and independent of the replication's truth. Gelman and Rubin (1992), Brooks and
  Gelman (1998) and Vehtari et al. (2021) ask for starts overdispersed relative to the posterior, and Biswas,
  Jacob and Vanetti (2019) start their Polya-Gamma logistic comparison from the prior. It reads whether a sampler
  that reaches small k from above, as from the default start, reaches it as fast from either side. On the first
  25 replications of each regime.

**Streams.** On each replication, 19 independent chains of each library from the start, as one stream with no
burn-in, every other setting bart()'s default (75 trees, chi(1.5, 2), `n.threads = 1`); the effect regimes run
the same 19 chains, not their 10. The replications are step 2's own, the first 50 of each regime in its order,
then up to 100: the same data, truths and row labels, so no selection runs again and every label comes from step
2's fits, never from the truth. Each stream is read at a grid of sweeps, every sweep to 50 and then about 40
points to each tenfold length, scored as it runs, so only the truth's ranks, the reference masses and 20 rows'
values of each class are kept. A library's stream in a regime runs to 2,000 sweeps, and on to 8,000 only where its
reading has not settled (below) by 2,000. Past 8,000 the reading is "not reached by 8,000 sweeps", read beside the
k study's long runs, which at 5,000 sweeps of burn-in and 20,000 recorded did not reach agreement between the two
samplers in most logistic fits with a small true k.

**Quantities.** Every class step 2 scores that is a function of one state of a chain: k; the latent fit at the
scored rows, interior and near 0 or 1, as step 2 labels them; the log risk ratio in S2; individual effects,
subgroup effects, the average effects overall and on the treated, the log risk ratio and the log odds ratio in
S7 and S8. Bands are left out: a band is a statement about the joint spread of many draws, which one state of each
of 19 chains does not give.

**Reading 1, against the truth (every replication).** Under step 2's truths, drawn from the model's prior, the
truth given the data is one exact draw from the posterior, independent of every chain run on those data. If at
sweep t the chains draw from the posterior, the truth and the 19 chains' states at t are 20 independent draws from
one distribution, so the truth's rank among them is uniform on 0 to 19, and it lies strictly between the smallest
and the largest of the 19 with probability exactly 18/20 = 0.90. That holds on any data, under any selection that
reads the data alone, and whatever the autocorrelation along each chain, since the 19 are separate chains. So at
each grid sweep t and each class, c(t) is the share of truths inside the 19 chains' range, pooled over rows and per
fit as step 2 weighs them, and d(t) its noise-corrected distance from 0.90 (dec-B455). A sampler still on its way
shows d(t) away from 0, most often under nominal, its chains bunched away from the truth. The rank histogram over
its 20 bins, with its ecdf band (Sailynoja et al. 2022), is reported beside it at 500, 1,000 and 2,000 sweeps and
at the stream's end. Its limit is SBC's: it reads the average over a regime's datasets, so a sampler wrong in
opposite directions on different datasets can read near nominal. Reading 2 reads datasets one at a time.

**Reading 2, against the k study's long runs, where they serve.** The k study stored its simulated datasets with
their truths and, for each long run that passed acceptance (the two samplers' long runs agreeing, dec-B462), the
pooled draws as 1,000 quantiles a row and k's draws in full. Among the regimes read here these are 2 of S1's 12
long runs, 4 of S2's 12 and 9 of S4's 24, 15 in A, and 3 of 6 in S5, 6 of 6 in S6, 5 of 6 in S7, 4 of 6 in S8
and 8 of 12 in E3, 26 in B. On each, 19 chains of each library from the default start, as in Reading 1, and at
each grid sweep the reference's mass inside the 19 chains' range, per row and for k; at stationarity its
expectation is 18/20 = 0.90 again, now on one dataset, the noise of the truth's single draw replaced by the
reference's mass. d(t) as above, per dataset and averaged in each regime. It is low in noise but covers only
datasets on which the long runs agreed, which in A are the easier ones: in the k study a logistic long run with a
true k below 0.5 passed 1 time in 16. Long runs that failed are not used, since a reference from the move's half
alone is the sampler under test (dec-B456). The references were made on the k study's engine; the posterior does
not depend on the engine as long as main fits the stored data under the same model, which the harness checks on
each dataset before use (the cut points and the prior main builds against those the k study's harness builds at
its recorded commit), dropping a dataset that differs.

**Reading 3, agreement between chains.** Rank-normalized split-R-hat (Vehtari et al. 2021) over a replication's 19
chains, on the second half of each prefix as Brooks and Gelman (1998) iterate it, for k and the 20 kept rows of
each class, and over both starts' 38 chains where both ran. It shows what a user's diagnostic says at each length
beside the readings that know the answer: R-hat near 1 while d(t) stands away from 0 is a shared slow direction
that R-hat does not see. It is plotted, 1.01 marked, and has no threshold of its own here.

**When a reading has settled, in sweeps and seconds.** For each library, regime, start, class and weighting,
t_sweeps is the first grid sweep from which d(t), smoothed over the grid points within a factor of 1.25 of t,
stays at or below 0.02 to the stream's end; a curve that never does reads "not reached" at the stream's length.
0.02 is pricing's noticeable coverage difference, the difference step 2 is sized to resolve; published thresholds
read other distances (0.25 of total variation for a mixing time, 0.01 for an honest burn-in, 1.01 for R-hat), so
the bar on this scale is the study's own, and the curves are reported whole so that any other can be read.
t_seconds is t_sweeps times the library's seconds per sweep in that regime from step 2's timing on the idle box,
the move's whole cost charged, and each stream's own clock is recorded beside it. The curves are drawn against
seconds as well as sweeps, both libraries on one axis, the form of Biswas, Jacob and Vanetti's (2019) figure 4.
d(500) and d(1,000) are reported too: how far each sampler stands from its stationary distribution where a default
fit starts and stops drawing.

**Replications and noise.** Intervals bootstrap the replications, 4,000 resamples, every reading of a replication
resampled together; t_sweeps is recomputed in each resample, which gives its interval. Planning spreads are step
2's per replication (Sizing and stops): a 19-chain range at one sweep is not a 4-chain interval of 500 draws, and
its spread for the row classes is assumed equal and measured at the first look. At 50 replications the per-fit
standard error of c(t) for the interior is then about 0.011 to 0.031 at spreads of 0.076 to 0.22, about half that
pooled; at 100, 0.008 to 0.022. k is one value a replication, its range indicator spread like a 0.9 coverage
indicator, 0.30, so its c(t) is known to about 0.04 at 50 and 0.03 at 100: k's reading leans on Reading 2 and the
rank histogram, and it is a diagnostic beside the row classes, which are what users estimate. What would change a
reading is the factor between the two libraries' t_sweeps and where each stands against 500 and 1,000 sweeps, so
t_sweeps is needed within a factor of 2, not to the last sweep.

**Stop.** A first look at 50 replications a regime, both libraries, default start. A library's reading in a
regime stops there when every row class's t_sweeps interval spans at most a factor of 2, or reads "not reached"
with its lower end at the stream's length; otherwise the regime runs to 100 at the lengths the first look needed.
No replication is added past 100. The dispersed start stays at 25, its streams lengthened by the same rule.

**How it is read**: beside step 2's readings, framing nothing; the report gives each sampler's curves and its
t_sweeps and t_seconds by regime, start and class (call 11).

## Step 3: a check on the data, if the maintainer sends pricing there

Step 3 runs only if the maintainer, having read step 2's report, sends pricing here (dec-B624); step 2's
readings recommend it when the cost is too high.

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
cell, an empty spared set meeting the condition; the single-feature form takes a tie. The k study's 0.8 was
chosen after its data and is not carried over.

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

**Reading, fixed now.** The trigger passes when it fires on every R3 dataset that needs the move and spares at
least half of those that do not. Below half, a second kernel would buy back too little of the cost to be worth
carrying. With fewer than 3 resolved datasets on either side, sensitivity or specificity cannot be read and
step 3 reads unresolved, reported with the counts. With 12 datasets this test is coarse: it can show a trigger
that misses a dataset needing the move, and little more, and the report says so.

**After the reading.** Step 3's results go to the maintainer as a research report ([The report](#the-report)),
the trigger's firing and sparing shown dataset by dataset, and the maintainer decides (dec-B624). A pass
recommends the move on where the trigger fires; a fail or an unresolved reading recommends step 4.

## Step 4: an opt-in

Only if the maintainer, having read the reports, sends pricing here, and step 1 passed. The opt-in's form (a
logical `dbartsControl` argument fixed
at creation, default off, as the probit step's withdrawn switch was, about 120 lines; or a family setting) and
whether to offer it at all go to the maintainer with steps 1 to 3's evidence. The recommendation is written then,
against what they show.

## Sizing and stops

Per-fit seconds are the k study's for a shipped logistic default fit on the box (its cost per accepted
replication). A step 2 replication costs about 2.2 shipped fits (equal time and the move, at r near 1.1), plus
20 percent for scoring and start-up and 15 percent for shared memory bandwidth. Spreads are the k
study's per replication, of the paired difference in distance from 0.90; where a regime was not in the k study
the spread is assumed and marked. Harmless must hold under each weighting, so each row is sized on the worse of
the two. In the k study's problem cells the per-fit standard error of the interior ran 1.4 to 2 times the pooled
one (near separation 0.024 against 0.012, small samples 0.007 against 0.005, the sampler loop 0.008 against
0.005), and near 0 or 1 at most equal to it. So the interior, individual effects and subgroups are sized per fit,
at twice the pooled spread where the ratio was not measured, and near 0 or 1 pooled. "Resolves" means harmless
at a true difference of 0 with probability 0.9: the upper bound at one-sided 0.025 below the bar, which needs a
standard error of at most the bar over 3.24.

| component | question it serves | size (looks) | cost, core-hours | precision at full size | smallest difference that changes the call | minimum size that resolves it | sequential stop |
|---|---|---|---|---|---|---|---|
| step 0, correctness gates | is the rebuilt move right and inert off | as listed | about 5 (laptop and box) | gates | any gate failing | as listed | none: a failure stops |
| step 1 pilot | T, B and the score distribution | 20 a regime and prior, run long | about 1 | sets T | | | |
| step 1 power check | can the tilted arm see the move's errors | 150, two mutant libraries | about 1 | flags or not | either mutant unflagged | 150 | none |
| step 1 SBC | is the small-k tail right | 600 a prior, 99 draws, thin T (150, 600) | about 5 at T 200; 10 at most | the README's detection floor, checked by the mutants | half the shape of k's conditional | 600 | fail at 150 on p below alpha / 10, confirmed by a rerun |
| step 1, third arm (contingent) | the tail away from near-determinism | 600, as above | 2.5 at T 200; 5 at most | as above | as above | 600 | as above |
| step 2, S1 | gain, no harm | 200 (100, 200) | 0.45 | per fit, interior se 0.014 at a spread of 0.20; near 0 or 1, pooled, 0.014 at 0.20 | 0.02 gain; 0.02 harm | 1,050 for each at true 0; 100 for the expected gain | first look, once resolved |
| step 2, S2 | gain, no harm | 480 (240, 480) | 1.15 | per fit, interior se 0.010 at 0.22 (measured; pooled 0.135) | as S1 | 1,280 interior at true 0; 480 at the expected improvement | as S1 |
| step 2, S4 | gain, no harm | 200 (100, 200) | 0.2 | per fit, interior se 0.005 at 0.076 (1.4 times pooled) | as S1 | 152 interior | as S1 |
| step 2, S3 | no harm at large n | 100 (50, 100) | 3.2 | per fit, interior se 0.012 at an assumed 0.12 | 0.02 harm | 380 interior; 213 to 3,600 near 0 or 1 | as S1 |
| step 2, S5 | no harm | 270 (135, 270) | 1.05 | per fit, interior se 0.006 at an assumed 0.10; near 0 or 1, pooled, 0.005 at 0.09 | 0.02 harm; bands reported with no bar | 263 interior; 213 near 0 or 1 | as S1 |
| step 2, S6 | no harm | 220 (110, 220) | 0.4 | per fit, interior se 0.006 at 0.086; near 0 or 1 0.006 at 0.09 | as S5 | 195 interior; 213 near 0 or 1 | as S1 |
| step 2, S7 | no harm on effects | 300 (150, 300) | 3.55 | per fit, individual effects, subgroups and interior se 0.006 at an assumed 0.10; near 0 or 1 0.005 at 0.09 | as S5 | 263; 213 near 0 or 1 | as S1 |
| step 2, S8 | no harm on effects | 300 (150, 300) | 4.05 | per fit, individual effects se 0.006 at 0.106 | as S5 | 296; otherwise as S7 | as S1 |
| step 2, E3 | no harm in a loop | 220 (110, 220) | 1.15 | per fit, interior se 0.005 at 0.08 (1.6 times pooled); near 0 or 1 as S5 | as S5 | 169 interior; 213 near 0 or 1 | as S1 |
| step 2, R2 held-out score | no harm on real outcomes | 6 datasets, 20 splits | about 1 | log score se about 0.002 a dataset | 0.01 nats | 20 splits | none |
| step 2, timing | r, and the cost bar | 9 regimes, 8 large cells, hazard and counts cells | about 1 as 100-sweep segments (3.5 as full fits), box idle | ratio to about 0.01 | 1.25 on a cell | 9 pairs; 10 segments | none |
| step 2, effective draws a second | reported | 10 a regime, S1, S2, S4, S6, both libraries | 0.3 | | decides nothing | | none |
| convergence, default start | sweeps and seconds to the stationary distribution, both libraries, 8 regimes | 50 a regime (50, 100), 19 chains, 2,000 sweeps, 8,000 where unsettled | 3.6 at 2,000 sweeps; 4.9 with the shipped sampler to 8,000 in A; 17 with every second look and every shipped stream to 8,000 | per fit, interior se 0.011 to 0.031 at 50, 0.008 to 0.022 at 100; k 0.04 and 0.03 | t_sweeps within a factor of 2; d(500) and d(1,000) | 50 where the curves stand apart, 100 where they are close | first look, once each t_sweeps interval spans at most a factor of 2 |
| convergence, dispersed start | does a sampler arrive as fast from a prior draw | 25 a regime, 19 chains, as above | 1.8; 2.5 with A's shipped streams to 8,000; 4.3 at most | about 1.4 times the default start's at 50 | read for its direction | | none |
| convergence, Reading 2 | the curves dataset by dataset on accepted long runs | 41 stored datasets, 19 chains a library | 0.4; 0.8 at most | per dataset | | | none |
| convergence, Reading 3 | what R-hat says at each length | the same streams | 0 | | frames nothing | | none |
| step 3, R3 long runs | which datasets need the move | 12 datasets, 8 chains a library | about 5 | per dataset | | 3 resolved on each side | none |
| step 3, R3 default fits | the default-length harm on each | 20 shipped at equal time and 20 move fits a dataset | about 1 | t bound a dataset | 0.02 | 20 | none |
| step 4 | | | 0 | | | | |

**Audit, regime by regime, at a true difference of 0, on the worse weighting.** The bar is 0.02 for every class
read for harmlessness, and a class passes only on a bound shown, so each count is set to resolve where it can.
- S1: the interior and near 0 or 1 each resolve 29 percent of the time at 200. Both are cells where the shipped
  sampler's harm is known (3 to 19 points), and an improvement resolves them at the first look; at a true 0 each
  would need about 1,050. Kept at 200.
- S2: at 480 the interior resolves 90 percent of the time pooled but 51 percent per fit, at a measured per-fit
  spread of 0.22. Per fit the k study's stored fits have the shipped sampler 4.6 points under and the move 3.0
  over, an improvement of about 0.016 in distance, at which 480 resolves 95 percent of the time; at a true 0 it
  would need about 1,280, 3 core-hours more. Kept at 480; near 0 or 1 is a gain cell, read as S1's.
- S4: resolves at 200 (interior 96 percent per fit); near 0 or 1 is a gain cell.
- S3: the interior resolves 38 percent of the time per fit at 100 (92 percent pooled); it would need about 380,
  13 core-hours more. Near 0 or 1, at spreads from 0.09 to 0.37, it would need 213 to 3,600, which is 5 to 170 more
  core-hours. Kept at 100: S3 resolves only through improvement, and otherwise goes to the maintainer with its
  counts.
- S5: raised from 220 to 270 (looks 135 and 270), which resolves the interior 91 percent of the time per fit and
  near 0 or 1 95 percent at 0.09.
- S6 and E3: 220 resolves the interior 93 and 96 percent per fit, near 0 or 1 91 percent at 0.09. Unchanged.
- S7, S8: raised from 200 to 300 (looks 150 and 300); individual effects, subgroups and the interior resolve 93
  and 90 percent per fit, near 0 or 1 97 percent at 0.09.
- B's near-0-or-1 spreads are taken at 0.09, the low end of the k study's 0.09 to 0.37. At 0.2 they would resolve
  about a third of the time at these counts.
- Bands carry no bar (dec-B625). At step 2's 135 to 270 replications a regime, one sampler's band coverage is
  known to about 0.04 to 0.06 either side at 95 percent, and the report gives it so.
- Showing a true harm of exactly 0.02 at 90 percent power needs about 420 at a spread of 0.10 and 1,700 at 0.20;
  under this reading a true harm that is not shown still leaves the class unresolved, so the move is not
  recommended on by it.
- The spreads for S3, S5, S7 and E3, and every per-fit factor of 2, are assumed. If a first look's standard error
  exceeds 1.5 times the planned value, power is restated before going on (step 2's stop).

**Audit of the convergence measure.** Its streams cost 19 chains times the stream's length at the k study's
per-sweep times on the box: 0.12 to 0.47 milliseconds a chain-sweep for the shipped sampler across the eight
regimes, the move 1.10 to 1.14 times that on its old placement; early sweeps, on smaller trees, cost less, which
the figures do not count.
- Length. The shipped sampler's streams in A are expected to run past 2,000 sweeps. A probe of the shipped sampler
  on main's engine, 20 chains of 3,000 sweeps from the default start, on three of the k study's stored datasets
  (two near-separated, one small-sample), had k's distribution within a Kolmogorov distance of about 0.15 of the
  reference only after about 1,000 to 2,000 sweeps on the two whose long runs passed acceptance, and on the
  third, with a true k of 0.02 and a reference from the move's half, still at a median k of 0.32 against the
  reference's 0.12 at 3,000. The move's time from the start is not measured; at stationarity its k's
  autocorrelation time was about 115 sweeps at the 99th percentile of its survey's census, against about 23,000
  for the shipped sampler, so its streams are expected to settle by 2,000 in most regimes. The cap is 8,000;
  20,000, the k study's length, would add about 3 core-hours expected and up to 20, and on the hardest datasets the
  k study's long runs already read past it.
- Replications. 50 place t_sweeps within a factor of 2 wherever the two samplers' curves stand well apart, as the
  probe and the k study's equal-time readings (shipped 0.78 to 0.85 per fit in A against the move's 0.89 to 0.91)
  lead one to expect in A; 100 where they are close, as they may be in B. The dispersed start at 25 is read for
  its direction, its intervals about 1.4 times the default start's at 50.
- Chains. 19 a replication make the range of the chains a 90 percent interval exactly at stationarity. Fewer
  chains with more replications would buy more precision for the same compute, but 9 chains give an 80 percent
  range, off the level every other reading of pricing uses.
- S3 is left out (about 6 core-hours at 2,000 sweeps, 15 with the shipped sampler's longer streams).

**Totals.** About 38 core-hours expected through step 2 with the convergence measure: step 0 about 5, step 1 about
7, step 2 about 17.5 (about 10 if its regimes stop at their first look), the convergence measure about 8. About
44 with step 3. About 70 at most, with step 1 at its longest thin and its third arm, the timing run as full fits,
and the convergence measure at every second look with every shipped stream at 8,000; about 14 more if the move
also needs its long streams everywhere. That is within dec-B530's rough guess of tens of core-hours. Against the
previous amendment (37 expected, 43 with step 3, 56 at most), step 2's fitting fell by a third, 25 to about 17.5,
with the burn-in-only arm dropped (dec-B621), and the convergence measure adds about 8 expected and 22 at most
(dec-B622). The stop is 66, 1.5 times the expected total with step 3. Under it, cuts go in this order: the
convergence measure's second looks, the dispersed start's long streams, the contingent arm's thin to its pilot's
T, S3 to its first look, the other step 2 regimes to their first looks, R3 to its first 8 datasets in score order.
A smaller design, B at 100 a regime and S2 at 200 (about 7 core-hours less), would resolve B's classes at a true 0
only about 50 to 60 percent of the time and S2's interior per fit 25 percent, so more classes would reach the
maintainer unresolved.

Where it runs: every timing, every fit of steps 1 to 3 and every stream of the convergence measure on the bench
box, at most four jobs at once, the two libraries' streams interleaved, the box idle for timing
([The x86 leg](../../benchmarks/README.md#the-x86-leg)); step 0's gates on the laptop, with its x86 leg on the
box.

## The report

Whatever a step shows, its results go to the maintainer as a research report, and the maintainer decides
(dec-B624). Step 2's report also carries step 1's results and the convergence measure; a step 1 that fails
reaches the maintainer on its own, and step 3, if it runs, has its own report. Each is to the standard for changing a default ([Studies](README.md#studies)) (dec-B559): the
maintainer set that standard for default priors and said mixing might need it too, and a change of kernel for
every logistic fit is a default users meet. Each report is a short PDF, at most 5 pages of text with graphs as
needed, self contained. Every report gives each sampler's empirical coverage against nominal, with its interval
and replications, widths, root-mean-square error and bias, for every class and band, and the difference between
the samplers with its interval (dec-B625); step 2's adds the convergence curves, t_sweeps and t_seconds. Past
that, its contents are taken from published comparisons of MCMC samplers, read when the report is written:
their readings, tables, figures and diagnostics. The report follows them and cites them, and names as the
study's own anything it shows with no published model, such as Reading 1 of the convergence measure
(dec-B626). The survey notes which published comparisons can serve
([Models for the report](../design/sampler-convergence-measures.md#models-for-the-report)). It says how the
regimes and datasets were chosen, and what 0.9-34 did. Two independent referees read it on a neutral brief
before it reaches the maintainer. Whatever the maintainer decides, the finding is then a design note in the usual
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
- **A burn-in-only arm.** Dropped (dec-B621): step 2 is the two mechanisms head to head, and no arm switches
  samplers within a run. The k study's open question about the start goes to the convergence measure, which
  reads each sampler's stream from the start with no burn-in (dec-B622).
- **Held-out scores as a gate, with the datasets named.** Adopted, amended. It is a no-harm guard in step 2 at
  0.01 nats a row, not a gain criterion: under logistic the held-out probabilities converge where k does not
  ([binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md)), so a gain there is not expected. It runs
  on R2, six datasets none of the probes has seen, named by rule above; about 1 core-hour.
- **What "wins" means.** Adopted: noise-corrected distance from nominal on either side, by class, pooled over
  rows and per fit, a gain counting only under both weightings and a harm under either.
- **The no-harm margin against nominal, not the told fit.** Adopted: every bar compares the two arms'
  distances from 0.90, and no told fit enters.
- **How equal time is measured.** Adopted, under Common to steps 1 to 3.
- **A number for "cost too high".** Adopted, as rule 3 of step 2, with the two bars ruled (dec-B620).
- **Whether the probit half stays out.** Answered: out (dec-B617).
- **The note from dec-B555.** Taken as framing. Step 2 asks for a gain only on the cells where the k study showed
  harm, all on probabilities, and only no harm on effects, since a gain on effects is not what that evidence
  predicts.
- **Whether the move's code size is a gate.** Declined as a gate, reported. The code is written (about 1,600
  lines with its tests, less the probit half), and the shape, consistency and what users want decide, not the
  work.

## Calls for the maintainer

Ruled, 2026-11-07. The numbers are kept so that the entries and this plan's text still meet; call 4 fell away.
- **Call 1, the probit half.** Out of pricing and of 1.0-0, probit keeping the step it ships; pricing is logistic
  only (dec-B617).
- **Call 2, the move's reach under logistic.** The forests that share one prior total as well as a single forest
  (dec-B618); hazard.logistic and count-weighted logistic fits in reach, timed, their timings entering the cost
  bar, and not scored as a regime of their own (dec-B619).
- **Call 3, what "cost too high" is.** Two bars: no harm failing only in B while the gain holds, or any timed
  logistic cell above 1.25 times the shipped library (dec-B620).
- **Call 5, a switch if the move is on by default.** None, the engine keeping its toggle for its own tests;
  provisional, put to the maintainer again beside step 2's reading and its times before the move lands on
  (dec-B623).
- **Call 10, bands in step 2.** No harm bar and no label: the report gives each sampler's empirical band coverage
  against nominal, with its interval and replications, and the difference with its interval (dec-B624, revised
  by dec-B625). Whatever a step shows, its results reach the maintainer as a research report (dec-B624).

Open.

**Call 6. Where the move waits.** Options: (a) on its rebuilt branch until pricing decides, landing once in the
form decided (recommended, dec-B494's rule that a feature merges complete); (b) on main now with the engine
toggle off and no R surface, so pricing reads main itself, at the cost of an engine change on main whose fate is
undecided. Either way the bitwise check at landing applies. 0.9-34 shipped no engine toggle of this kind.

**Call 7. The held-out datasets.** The candidate list and the rule that splits it into R2 and R3 are above; you
may add, strike or reorder before step 2 runs, and nothing is fitted on them before then. No 0.9-34 study used
held-out real datasets for a sampler change.

**Call 8. The September prior study's revisit under logistic.** TODO binary-hyperprior-revisit waits for k's
scale move to land (dec-B390), and its first condition is that k mixes. Options: (a) read that condition free
from step 2's move fits at the default length, k's split-Rhat below 1.05 and its effective draws in the hundreds
(recommended); if it holds and the move lands for logistic, a logistic rerun of the study goes to you as its own
item after landing; (b) leave the revisit as it stands, a probit item, with nothing read here. The study was
probit, and with the probit half out (dec-B617) its probit condition is not met by this plan either way. 0.9-34
had no logistic family.

**Call 9. Negative binomial's item.** The move's plan leaves negative binomial to its own item, and there is
none. Options: (a) name TODO k-mixing-pg-families as that item, which already holds the Polya-Gamma families' k
mixing (recommended); (b) open a new item. 0.9-34 had no negative binomial family.

**Call 11. How the convergence measure is read.** The measure reads, for each sampler, how many sweeps and
seconds its chains need from the start a user gets before they draw from the posterior (dec-B622). Step 2's arms
both start there and keep the 500 burn-in, so whatever the start costs in calibration is already inside step 2's
readings; the measure says how much of it is the start, and whether 500 sweeps suit either sampler. Options:
(a) read beside step 2's readings: the report gives each sampler's curves, t_sweeps, t_seconds, d(500) and
d(1,000) by regime, start and class, and the recommendation is framed by step 2's three readings alone
(recommended). (b) a fourth reading in the frame: the move is recommended on for every logistic fit only if, in
no regime, it is shown to settle later in seconds than the shipped sampler (its t_seconds interval wholly above
the shipped one's); otherwise step 3 is recommended. The compute is the same either way, about 8 core-hours
expected and 22 at most. For (a): the distance read against the truth, Reading 1, and its 0.02 threshold are
the study's own, with no published model, and under dec-B624 the maintainer decides from the report with the
curves in it either way, so a fourth reading adds a frame the survey cannot back. For (b): it holds the move to
reaching the posterior no slower in time, a claim the equal-time readings make only at the default length. A
recommendation on the default burn-in itself is not offered: dec-B622 keeps it, and the report's d(500) shows
whether it should be reopened. 0.9-34 had no logistic family and no move; its bart() burned in 100 sweeps with
200 trees, 1.0-0's 500 with 75.

## Left out

- **Probit**, every probit-only regime (S9, D2, E1, E2) and the probit half (dec-B617).
- **Negative binomial.** Its item is TODO k-mixing-pg-families (call 9). The negative binomial default study
  measured k under a scratch build of that half, and its default is still open.
- **The multiple-treatments designs.** Their miss is shared by both samplers and points to the trees and the model,
  and every design with a real effect failed its first long run.
- **The colleagues' case.** It has not been shared; it has a slot in the k study's plan (dec-B529).
- **Real outcomes scored by long-run agreement.** They carry no truth, and their 22 datasets are spent; real
  outcomes enter only as R2's held-out score and R3's trigger test.
- **Count-weighted logistic and hazard.logistic** are reached by the move, gated by its exact gate's counts arm and
  timed, but not scored as regimes (dec-B619): the likelihood is the same, and the counts were not a k study
  regime.
- **S3 in the convergence measure**: its streams at 10,000 rows would cost about 6 to 15 core-hours, and it is in
  neither of the sets the measure reads.
- **A coupling of the two samplers' chains**, the published route to an upper bound on the distance from
  stationarity: none exists for a BART sampler, and one would be new work with its own proof
  ([sampler-convergence-measures.md](../design/sampler-convergence-measures.md)).
- **A change to the default burn-in**: dec-B622 keeps it; the convergence measure's d(500) is reported.
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
- The convergence measure: its review checks Reading 1's argument (the truth an exact posterior draw given the
  data, independent of the 19 chains, so its rank uniform at stationarity) against the harness; a smoke that feeds
  the scorer independent exact draws from a known conjugate posterior in place of the chains, with truths drawn
  from its prior, reads the truth's rank uniform and the range's share at 0.90 within its noise; Reading 2's
  model-equality check passes on each stored dataset before use, and the datasets it drops are listed.
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
- The convergence measure's spreads: a 19-chain range at one sweep is assumed spread like step 2's intervals, to be
  measured at the first look.
- Its probe: shipped sampler only, on the laptop with a library built from main on 2026-11-03, three stored
  datasets; the move's time from the start is not measured, only its autocorrelation at stationarity.
- That main fits the k study's stored datasets under the model their references were made under; the harness
  checks it dataset by dataset.
- k's starting value of 2 at the default start, read from the engine and from the probe's first sweeps, not from
  any documented default.

## Amendments

**2026-11-07, the maintainer's rulings dec-B617 to dec-B626.**
- Calls 1, 2, 3, 5 and 10 are marked ruled with their entries (dec-B617 to dec-B620, dec-B623 to dec-B625);
  call 4 fell away (dec-B621). Calls 6 to 9 stay open, and call 11 is new.
- The burn-in-only library, its patch and its arm are gone from step 0, step 2's arms and readings, Sizing and
  stops, the report and the proposals (dec-B621); step 2's fitting falls by about a third.
- The rules state what a reading shows and recommend; every step's results reach the maintainer as a research
  report, which gives empirical against nominal coverage, widths, root-mean-square error and bias for every class
  and band (dec-B624). The sequence, step 2's "After the reading", step 3's outcomes, step 4 and The report are
  rewritten so; step 3 runs only if the maintainer sends pricing there. A harm shown no longer ends step 2.
- Bands carry no bar: the 0.035 worsening is gone from rule 2, Sizing and stops and call 10, and each sampler's
  band coverage is reported with its interval (dec-B625).
- The report's contents are taken from published comparisons of MCMC samplers, read when it is written, anything
  without a published model named as the study's own (dec-B626).
- New: [The convergence measure](#the-convergence-measure) (dec-B622), after a survey,
  [sampler-convergence-measures.md](../design/sampler-convergence-measures.md); its rows in the sizing table and
  its audit; its checks under Verification; call 11 on how it is read.
- Totals: about 38 core-hours expected through step 2 (was 37), 44 with step 3 (was 43), about 70 at most (was
  56), stop 66 (was 60), the convergence measure's second looks cut first.

**2026-11-06, the blind critique.** The critique of the first draft made 22 points; each is listed by its number
with what changed.

1. Taken. The shipped arm at equal time takes ceil(1000 r) - 500 draws after the 500 burn-in, so the move's whole
   cost is charged; the burn-in-only comparator likewise at r_b.
2. Taken. No harm is now a bound shown: harmless when the upper bound of the worsening is below 0.02, harm
   when the worsening reaches the bar with its lower bound above 0, one-value estimands read for
   shown harm only, unresolved classes to the maintainer. The sizes were re-audited at a true difference of 0 on
   the worse of the two weightings (S2 to 480, S5 to 270, S6 and E3 to 220, S7 and S8 to 300); S3 is not sized to
   resolve and says what would. On the critique's confirmation, bands became one-value readings for shown harm
   only (call 10), the expected count of unresolved classes is stated, and an empty spared set meets step 3's
   trigger condition.
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
