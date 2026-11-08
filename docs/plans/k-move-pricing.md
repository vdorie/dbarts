# k-move-pricing: pricing the collapsed scale move under logistic

Status: planned, in critique (2026-11-06); amended after its blind critique, and on 2026-11-07 for the
maintainer's rulings dec-B617 to dec-B626 and the convergence measure, then for dec-B633 to dec-B635, then for
dec-B651, the real-data guard and the dataset rule, for its referee, then revised on that referee's report
([Amendments](#amendments)).
Nothing has run. The sequence is dec-B530's; every rule below is written before any fit of the step it reads.
The rules state what a reading shows and recommend; the maintainer decides from a research report (dec-B624).

agent: one opus implementer for step 0 (the move rebuilt on main), one opus review of it that runs the move's
mutants; one opus implementer for the harness, steps 1 to 3 and the convergence measure, run on the bench box;
one review of the harness and its smoke (the gate-only tier of [Process by risk](README.md#process-by-risk)),
which also checks each of the convergence measure's readings against its published form; one opus review of the
step 2 reading before its report is written; two independent referees on each report
([The report](#the-report)); one independent referee of the dataset rule and list before they are sealed
([k-move-pricing-datasets.md](k-move-pricing-datasets.md)).
rng: the study changes no package code. Step 0 is NEUTRAL with the engine's toggle off and NEUTRAL with it on for
every fit outside the move's reach; where the move reaches (a logistic fit drawing k under a finite prior scale)
it is POSTERIOR-CHANGING by the README's rule once it lands on, the stationary distribution unchanged by
construction.
window: after the k study, before 1.0-0 (dec-B530). Engine slices stay serial.
budget: step 0, the move's own diff (about 300 engine lines, 950 of tests/cpp, 290 of tinytest, 255 of its
exact gate) replayed onto main less the probit half, plus the conflict work below and about 250 lines of
tests/cpp rewritten for the move's new placement; the pricing harness about 1,300 lines of R, most of it reused
from the k study's, the convergence measure's streams and scoring about 300 of them, the real-data guard's
about 250. Compute about 49 core-hours expected through step 2 with the convergence measure, about 70 with step
3, about 152 at most with every contingent arm, every escalation and every second look, stop at 105
([Sizing and stops](#sizing-and-stops)).

## Goal

After this plan runs, the maintainer has chosen one of four outcomes for the logistic collapsed scale move,
from a research report on each step read (dec-B624), each outcome recommended by a rule written here first: on
for every logistic fit it can reach, with no trigger; on where a check on the data fires; an opt-in; or parked.
Each rests on the move as it would ship on main's engine, on calibration of what users estimate at equal
running time, and, read beside it, on how long each sampler takes to reach its stationary distribution from the
start a user gets (dec-B622).

## Summary for the maintainer

Pricing decides whether the logistic family's new move for the leaf scale is turned on for every logistic fit,
turned on only where a check on the data calls for it, offered as an option, or parked. The move does not change
the answer the sampler aims at, only how quickly a run gets there, so the question is whether, at the same running
time, intervals and estimates come out better with it where the shipped sampler struggles, and no worse anywhere
else. The move is first rebuilt on today's engine and shown to change nothing it should not; then shown correct
where the scale is small; then the two samplers go head to head at equal time on simulated truths placed on real
covariates, beside a measure of how long each takes to settle from its start; a check on the data follows only if
you send pricing there. Each step reaches you as a short research report.

Real outcomes enter on datasets nobody here has looked at. Twelve of them are scored by how far each sampler's
estimates fall from a very long run of the same model, in place of predicting held-out outcomes: since both
samplers aim at the same answer, that measures directly the only thing the move can change. For the same reason the
screen of datasets keeps only what bears on such a comparison, and 21 candidates from published benchmarks pass
where 7 passed before; the last three of the 24 needed come from a seeded draw from the UCI repository, screened
twice before the list is fixed. How well each sampler predicts the held-out outcomes is still reported beside it.
The very long runs are accepted only when both samplers' runs agree, and on the strongest-signal and the largest
data they often do not. The guard reads the move harmless only when every accepted dataset shows it harmless on its
own and the average over them stays within the margin, as published comparisons across datasets test; a harm on
any one dataset blocks. Twelve was set the way a study sets its sample size: the smallest number at which, allowing
for how uncertain the acceptance rates are, the guard is expected to reach a conclusion at least eight times in ten
when the move does no harm (about 8.3 in 10 at twelve, 7.5 at ten). The check on the data, if you send pricing
there, uses every dataset the sources yield, about 35 to 45; so few of them can show a need for the move that it is
built to catch a wrong trigger, not to prove a right one, and its evidence that a trigger finds the fits that need
the move comes from the simulated fits. The rule, the list and the sizes went to independent referees; they come to
you before anything is fitted on them. About 49 core-hours through the head-to-head comparison, 70 with the check on
the data.

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
  reports of its kind (dec-B626); the convergence measure rebuilt from published comparisons' readings and read
  beside step 2 (dec-B633); the September study's first condition read from step 2's move fits against published
  thresholds (dec-B634); the held-out datasets put to independent agents (dec-B635); every study brought, in
  order of priority, to what published research of its kind would do (dec-B637); and the agents left to
  redesign how the real datasets enter, the rule and list still coming back as dec-B635's decision (dec-B651).
- The dataset rule and list: [k-move-pricing-datasets.md](k-move-pricing-datasets.md). How published comparisons
  score a sampler against a reference posterior: posteriordb (Magnusson, Torgander, Burkner, Zhang, Carpenter and
  Vehtari 2024), Margossian and coauthors (2024), Zhang, Carpenter, Gelman and Vehtari (2022).
- How published comparisons measure the time to reach the stationary distribution:
  [sampler-convergence-measures.md](../design/sampler-convergence-measures.md).
- What the move is: leaves times a factor r and k divided by it, r drawn by one slice step from its conditional
  with the latents integrated out, taken on a coin of one in four, then the latents drawn at the scaled fit. It is
  a generalized Gibbs step on the scale group (Liu and Sabatti 2000), whose second theorem lets a step that leaves
  the factor's conditional invariant stand in for an exact draw only when the step is itself invariant along the
  group's orbits; the move's exact gate and step 1 test the result. Survey:
  [k-tail-samplers.md](../design/k-tail-samplers.md); logistic mixing:
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
  lindner, ACTG 175 and lalonde, and IHDP, ACIC 2016 and 2019 and Twins as designs. Steps 1 and 2 use some of them
  for covariates only, under drawn truths. Real outcomes enter only through two sets of datasets not yet looked at,
  R2 (step 2's real-data guard, 12 datasets) and R3 (step 3's test, every other passer), drawn by the rule in
  [k-move-pricing-datasets.md](k-move-pricing-datasets.md) (dec-B635, dec-B651). Its frame is published and curated
  by others: TabArena's binary tasks, then OpenML-CC18's, then the binary datasets TabArena's curators set aside
  for reasons that do not bear on a sampler comparison, then a seeded UCI census only if those leave the pool
  short. Its screen keeps what bears on a comparison of two samplers of one posterior: the outcome as its source
  defined it, real data, no source shared with a seen dataset and one dataset per source, at least 100 rows and a
  0.01 minority share, a license permitting research use, and at most 2,000 columns. On documentation 21 pass, 14
  of them with at least 1,000 rows, so the census runs, on documentation before the seal, until 24 pass; on one
  screen its first 80 draws give 4 passers and 3 more whose share a second screen may settle, and the second
  screen, read so far on three, passes Toxicity and Bone marrow transplant: children, so 23 pass and it reads on.
  Each passer's separation score is computed from the data alone and kept sealed in the harness until step 3
  starts, since it is the trigger's own feature; R2 takes one dataset from each of 10 blocks by score, by a seed
  fixed in the rule, and R3 is every other passer, no passer held in reserve (sizes set by method, in call 7). The
  rule and list are sealed after their referee and before any file is loaded, and the assignment is committed
  before step 2 runs.
- **Equal running time.** Per regime, the move's time ratio r (move over shipped) on the box, idle: a full default
  fit as the user runs it, three replications a regime, three alternating repetitions per library, per-round
  minima, the load stamped, each less the same library's time for a fit of one sweep on the same data, the time
  that does not scale with sweeps (the model frame, cut points, start-up); r is the median of the nine pairs of the
  remainders, the ratio of the two libraries' seconds per sweep (call 20). The share of the fixed time in a whole
  fit is reported. The shipped arm at equal time takes ceil(1000 r) - 500 draws a chain after the same 500 burn-in,
  so its total sweeps take the time the move's take, which a user running both pays beside the same fixed time; a
  ratio of whole fits would pull r toward 1 and give the shipped arm fewer sweeps than equal time buys. This keeps
  the pilot's frozen choice of spending the extra time on draws ([What was
  compared](../design/binary-interval-pilot.md#what-was-compared)). Its default reading is the first 500 draws of
  the same run, by the prefix property. The sampler loop (E3) is timed by total fit time less a one-sweep fit's,
  its per-sweep R work counted as per-sweep. These timings are also the cost benchmark of step 2, together with
  time-only cells at 10,000 and 50,000 rows, 75 and 200 trees, both signals, from the move plan's large-cell
  generator, timed as the move plan's early measurement 1 did: both libraries restored from one stored state before
  every timed segment of 100 sweeps, ten interleaved repetitions in alternating order. Person-period cells of
  hazard.logistic and a count-weighted cell are timed the same way, their timings entering the cost bar (dec-B619).
  The same timings give each library's seconds per sweep, which the convergence measure's readings in seconds use.
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
   d(shipped at equal time) is *harmless* when its upper bound at one-sided 0.025 a look is below 0.02 for the row
   classes, individual effects and subgroups (the house's noticeable difference, the k study's one bar). *Harm* is
   shown when the worsening is at least that bar with its lower bound above 0 at z 2.81. One-value estimands
   (average effects, ratios) are read for shown harm only, as dec-B533 reads them on calibration. Bands carry no
   bar and no label (dec-B625): the report gives each sampler's empirical band coverage against nominal, with its
   interval and the replications behind it, and the difference between the samplers with its interval. The ratios
   of root-mean-square error and interval score, move over shipped, are harmless when their upper bound is below
   1.02. On R2's real outcomes, each sampler's error at equal time against a reference run of the same posterior is
   harmless when the upper bound of the move's worsening is below the error of 100 independent posterior draws,
   read over the datasets with an accepted reference on every row estimand and class, with a confirmed harm on any
   one dataset a harm; otherwise the guard is unresolved and reads as an unresolved class ([The real-data
   guard](#the-real-data-guard), dec-B651). The move is recommended on only when every class is harmless or
   improves. Any class unresolved at full counts is reported with the counts. This is the burden the binary
   interval study put on its guards: a class passes on a bound shown, never on a harm not seen ([Decision
   rule](binary-interval-study.md#decision-rule)). A harm shown frames the recommendation only once confirmed: its
   class is rerun on fresh seeds named before the rerun, at the count it was read at (the look's count; the guard's
   as [The real-data guard](#the-real-data-guard) gives it), and the harm stands when the rerun shows it too, as
   step 1 confirms its failures; otherwise the class is unresolved. At a true 0 about 50 class and weighting cells
   each at most 0.005 over two looks bound an unconfirmed false harm somewhere in A near 0.09, and a confirmed one
   near 50 x 0.005 x 0.0025, about 0.0006; a harm in A recommends parking (call 23).
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
- Gain, with a harm shown only in B or confirmed only on R2's real outcomes, or with the cost bar: the cost is too
  high, and step 3 is recommended (call 25).
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
k study; and k's rank-normalized split R-hat and bulk and tail effective draws from the move's fits at the default
length in every regime, against R-hat below 1.01 and at least 400 effective draws (Vehtari et al. 2021). That is
the September prior study's first condition read under logistic: if it holds and the move lands for logistic, a
logistic rerun of that study goes to the maintainer as its own item after the landing (dec-B634, call 8).

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

## The real-data guard

**Question.** On real outcomes, at equal running time, does the move make any estimate a user reads worse than
the shipped sampler does? Real outcomes carry no truth, so each fit is scored against a long reference run of the
posterior both samplers target, as comparisons on posteriordb's reference posteriors are scored (Magnusson et al.
2024) and as Margossian and coauthors (2024) score a run's error. This replaces the held-out log score as the
guard (dec-B651): the move leaves the posterior unchanged, so any difference between the samplers on a real
dataset is in how well a finite run approximates that posterior, which a reference reads directly and a held-out
score reads only through the outcomes' own noise. The held-out score is still reported, framing nothing (below).
Calibration, coverage against a known truth, stays on step 2's drawn truths (call 17).

**Datasets.** R2, twelve datasets (sized in [Sizing the real-data stages](#sizing-the-real-data-stages)), drawn by
[k-move-pricing-datasets.md](k-move-pricing-datasets.md): per dataset a
seeded subsample of at most 5,000 rows, 80 percent fitted (at most 4,000, the plan's cap) and 20 percent
prediction rows, whose outcomes stay sealed until the guard's reading is committed. Missing predictors are kept,
as bart()'s default `na.keepPredictors` keeps them.

**The reference.** Per dataset, 16 chains, 8 of each library, each started from its own draw from the model's
prior (trees, leaves and k, as the convergence measure's dispersed start; Vehtari et al. 2021 ask for diffuse
starts), 5,000 sweeps of burn-in, then recorded at thin 20 in blocks of 4,000 sweeps: the k study's form, with
its length guard and its escalation ([References](binary-k-mixing-need.md#references) and its length guard in
[Pre-registration](binary-k-mixing-need.md#pre-registration)), which step 3 uses. The first job records 16,000
sweeps, the checks read at 12,000 and 16,000; a reference is accepted only when checks 2 to 5 below pass at two
consecutive block ends, so that a lucky single block cannot pass, and check 1, the draw count, at the later one,
the end it is accepted at (12,000 sweeps hold 9,600 draws); a failure adds a block, to 20,000. Past 20,000
only the shipped half doubles, to 40,000, 80,000 and 100,000, thinned to keep its draw count, and only while
effective draws are the only checks failing and the first pass's effective draws per sweep project a pass at
100,000; the move half stays at its length. Accepted, the k study's acceptance as written, with posteriordb's draw
count:
1. at least 10,000 draws, posteriordb's first criterion (12,800 at 16,000 sweeps);
2. the two libraries' halves agree on every scored quantity's mean and 5th and 95th percentiles (the latent fit
   and the probability at each fitted and prediction row), and on k's 5th, 10th, 30th, 50th, 70th, 90th and 95th
   percentiles, within 4.5 combined Monte Carlo standard errors: a reference counts only when both samplers' long
   runs agree and converge (dec-B462);
3. rank-normalized split R-hat over all 16 chains at most 1.02 at the 99th percentile of the row quantities;
4. tail effective draws at least 400 at the 1st percentile of the row quantities (the thresholds of Vehtari et
   al. 2021);
5. k's effective draws at least 100. k's own R-hat is reported and not required: requiring it selects the
   datasets where the shipped sampler already mixes (the k study's call); the binary interval pilot, which
   required it, accepted 18 of its 52 references ([binary-interval-pilot.md](../design/binary-interval-pilot.md)).

The k study's limits on one-value quantities (R-hat at most 1.01, tail effective draws at least 400) gate nothing
here, since the guard scores none: its one-value quantities are k, which check 5 reads, and the average probability
over the prediction rows, which frames nothing and has its R-hat and effective draws reported beside it. Over 16
chains a one-value R-hat of 1.01 is a demand for precision (below), on a quantity that would frame nothing.
posteriordb's approximately independent draws, a mean lag-1 autocorrelation under 0.05, is reported per estimand
and not required (call 13); its two criteria specific to Hamiltonian Monte Carlo do not apply. The reference is
made by the two pricing libraries, so Reading 3's check against another engine's model does not arise.

**How many are accepted, written before any run.** R-hat at 1.01 over 16 chains is mostly a demand for precision:
the k study measured that a stationary chain passes it a quarter of the time at 1,200 effective draws, two-thirds
at 1,600 and reliably only at 2,000 to 2,400. Under this acceptance 12 of the k study's 22 real prediction
references failed: 10 of the 12 scoring at least 0.8 and 2 of the 10 below; 6 of the 7 with at least 1,000 rows and
6 of the 15 smaller. The k study records each of those datasets read under both links or under neither, so a
failure there is a failure under either link; logistic alone passes at least as often, and in that respect the
projection is conservative. R2 is drawn by score blocks; with acceptance read from that record by score and size,
and 14 of the pool's 24 at least 1,000 rows, about 4.2 of its twelve references are expected to pass, agreement's
chance failures below counted. The rates are soft, and the sizing carries their uncertainty ([Sizing the real-data
stages](#sizing-the-real-data-stages)). Agreement adds chance failures the k study's sizes hid: a sound reference
is tested at 4.5 standard errors on about 4 independent summaries a row (the probability's percentiles are the
latent fit's, transformed), which with exact standard errors fails by chance with probability about 0.1 to 0.2 at
5,000 rows, 0.03 at 1,000 and under 0.01 at 200, and more with noisy quantile standard errors (the k study's record
gives a sound run 2 to 22 percent at its sizes). The sizing counts this as a loss of 0.15 of the accepted
references on datasets of at least 1,000 rows, the middle of that range.

**Reading over datasets.** On each row estimand and class the guard reads harmless only when both of two readings
hold. First, every accepted dataset reads harmless on its own Welch bound (Reading, below): an intersection-union
test, valid at each dataset's own level with no correction and no assumption across datasets, which reads each
dataset as step 2's rule 2 reads each regime, so a harm on one dataset is not averaged away by gains on others.
Second, the claim over datasets, framed as Boulesteix, Hable, Lauer and Eugster (2015) frame a comparison over real
datasets: each accepted dataset a draw from a population of datasets, its D an estimate carrying its own error, the
upper one-sided 0.025 bound of a one-sample t over the accepted datasets' D, unweighted, below the margin of 1.
Every row estimand and class must read so, an intersection of tests that needs no correction. The t assumes D
normal across datasets, which a few datasets cannot check (Demsar 2006 asks for about 30 datasets or reason to
believe normality), so its level holds for symmetric spreads. Where a mean D of 1 sits in a few datasets, the t
alone reads harmless far more often than 0.025 (the referee's simulation: about 0.1 at a gamma spread of shape 1,
about 0.2 at shape 0.3), and the per-dataset readings catch the harmed datasets that are drawn; when none is drawn,
every dataset's own claim is still right and the claim over datasets is not, which is why the bound on the share
harmed is given beside it. Beside, framing nothing: the exact one-sided 0.025 upper bound on the share of such
datasets harmed (with all J read harmless, 1 - 0.025^(1/J): 0.71 at 3, 0.52 at 5, 0.37 at 8), and the Wilcoxon
signed-ranks test of D below 1, which Demsar recommends over the t and which can reach one-sided 0.025 only from 6
datasets.

The population the claim covers: R2 takes one dataset from each of equal blocks of the pool, so its accepted
datasets are close to an equal-probability sample of the pool's datasets whose references would be accepted, and
the claim is for that finite set; carrying it to datasets like these beyond the pool rests on Boulesteix and
coauthors' idealization that datasets are drawn at random. References fail where the shipped sampler is slowest, so
the guard reads the datasets both samplers can reach, the fits that do not need the move, which are the ones a
no-harm guard is for. Acceptance requires both libraries' long runs to agree, so a dataset where the move's own
long run misbehaves drops out rather than reading as a harm; the report shows, for each failed reference, which
half failed which check. The datasets are shown together as a forest of intervals. There is no floor on the count:
the intervals set how many datasets a reading needs. At a true difference of 0, two accepted datasets read harmless
about 5 times in 100 (the t on one degree of freedom is wide), three about 96 times in 100 and four or more about
98; at a true difference of 1 on every dataset, never. What the guard concludes at each count:
- none accepted: nothing; unresolved, the two samplers' fits shown side by side with the held-out score;
- one: that dataset's own reading, harmless or not, and no claim over datasets (a t needs two); unresolved;
- two or more: harmless when every dataset and the t read so, otherwise unresolved, each dataset's own reading
  shown;
- at any count, a confirmed harm on a dataset is a harm, and a harm shown and not confirmed leaves the guard
  unresolved.

Unresolved, rule 2 of step 2 reads it as an unresolved class: the move is not recommended on by it. A dataset whose
reference fails is not replaced. Where a reference fails, the report shows the two samplers' fits side by side,
with the held-out score, and draws no conclusion from them.

**What the reference's tolerance allows.** The halves may disagree by up to 4.5 combined standard errors, about
0.2 to 0.28 posterior standard deviations at 1,000 to 500 effective draws a half, so the pooled reference can sit
up to about 0.1 to 0.14 standard deviations toward one half: 1 to 2 units of the error below, the size of the
margin. That offset cancels from D when both samplers' fits are unbiased for the posterior; it favors a sampler
whose short fits lean the way of its own long half, the shipped one (making the guard conservative) or the move
(making it lenient). So D is also read against each library's half alone, at no compute. A harmless reading is
flagged when it does not hold against the shipped half, the half that can only favor the shipped sampler; a harm
stands only when it also holds against the move's own half, the half that can only favor the move (call 22).

**Scored fits.** Per dataset, 40 fits of each library from the default start on independent seeds, bart()'s
defaults with 4 chains: the move at the default length, the shipped sampler at equal time, ceil(1000 r) - 500
draws after the same 500 burn-in. r is the dataset's own ratio of seconds per sweep, measured as step 2 measures a
regime's, on nine alternating repetitions of a default fit per library, each less a one-sweep fit's time ([Common
to steps 1 to 3](#common-to-steps-1-to-3), call 20), the timed fits predicting at the prediction rows as the scored
fits do, since prediction's cost scales with the draws kept. Each fit predicts at the prediction rows. 80 fits a
library on a dataset whose own spreads call for them (Power, below).

**Estimands.** At every fitted and every prediction row: the latent fit's posterior mean, its posterior standard
deviation, its 5th and 95th percentiles (the ends of the 90 percent interval step 2 reads), and the posterior mean
of the probability, the prediction a user reads. posteriordb names means, variances and tail quantiles as the
expectations of interest. Rows are read together and by step 2's classes, interior and near 0 or 1, labeled by the
k study's thresholds from the reference's posterior mean of the probability, a function of the data and the
reference and never of a scored fit. Fitted and prediction rows are read together in each grouping, so a dataset
has 15 tests (five row estimands, three groupings: all rows, interior, near 0 or 1); the prediction rows alone are
reported beside. Beside them, framing nothing: k's posterior mean and 5th and 95th
percentiles, and two joint estimands, the average probability over the prediction rows and the reference mass
inside each fit's 90 percent band, which a sampler with the right margins and the wrong dependence across rows
would miss (call 24).

**Error.** For each fit, estimand and row, the squared error against the reference, in units of the expected
squared error the same estimate would carry from 100 independent draws of the reference, found by drawing 100 of
the reference's pooled draws 2,000 times. For a posterior mean the unit is the posterior variance over 100,
Margossian and coauthors' precision target for one estimator ("the expected squared error attained with 100
independent samples"); its use as a margin on a difference between two samplers, and for the standard deviation
and the interval ends, is the study's own. A fit's error on a class is the average over its rows. A default fit
with E effective draws errs about 100 / E in these units on a mean. The reference's own error adds the same amount
to both samplers' errors, so it cancels from their difference.

**Reading.** On each dataset, for each row estimand and class, D is the mean error of the move's fits minus that of
the shipped fits. Harmless when the upper bound of D, one-sided 0.025 on Welch's t, is below 1; harm shown when D
is at least 1 with its lower bound above 0; one look. A fit's error is a sum of squared errors, skewed, nearly one
chi-square where the rows move together; at a true D of 1 Welch's bound then reads harmless more often than 0.025,
about 0.05 at a gamma shape of 0.5 and 0.035 to 0.045 at shape 2 with 20 fits, about 0.036 at shape 2 with 40, and
nominal at shape 8. A percentile bootstrap is narrower and worse, and requiring it beside Welch is Welch alone, so
it is not used; a studentized bootstrap helps only at moderate skew (call 26). Each library's skewness of the
error over its fits is read in the power check's pass, and the report states each test's level at the skew
measured. That is each dataset's own reading, which the guard requires of every accepted dataset (above).

**A harm shown is confirmed** before it frames anything: as many more fits of each library as it was read at (40 or
80) on fresh seeds named before the rerun, against the same reference, read on the estimand and class that showed
it, as step 1 confirms its failures. The harm stands when the rerun shows it too and it holds against the move's
own half as well (above); otherwise it is reported unconfirmed and the dataset's reading is unresolved, which still
blocks, so confirmation never makes a harmless reading. About 0.2 core-hour a dataset at 40. A false harm at a true
0 has probability about 0.0004 a test at s = 0.94 and 0.0125 at 1.41 over about 90 correlated tests, and needs both
independent runs once confirmed, about 0.0004 squared a test. A confirmed harm on R2 blocks the move on for every
logistic fit and recommends step 3 ([Step 2: it wins at equal time](#step-2-it-wins-at-equal-time), After the
reading, call 25): the datasets whose references pass are mostly fits with no strong signal, which is what B stands
for.

**How a slower sampler can pass**, stated so that no reading is taken for more than it is. The margin is absolute,
so a move with a quarter of the shipped sampler's effective draws on the rows (errors 0.2 against 0.8, D = 0.6)
passes wherever s is small; the error ratio, reported and flagged past 1.25, is where that shows. k is outside the
guard. Datasets whose references fail are not read, though the held-out score is. Joint estimands frame nothing.

**Beside, framing nothing**: each sampler's error by estimand with its interval; the ratio of the two mean errors,
move over shipped, with its 95 percent bootstrap interval on every estimand and class, any whose lower bound passes
1.25 flagged (dec-B620's cost bar restated as accuracy: for converged chains the ratio at equal time is the inverse
ratio of effective draws per second, which posteriordb calls "the central implementation-dependent performance
metric"); the
1-Wasserstein distance between each fit's draws and the reference at each row, scaled by the reference's standard
deviation (posteriordb's second measure; the per-row form and its scaling are the study's own, since Zhang et al.
2022 compute it on the whole parameter vector and scale by a baseline's distance); the reference mass inside each
fit's 90 percent interval, which step 3 reads; a count over datasets of the sign of D (Demsar 2006); each sampler's
effective draws per second on the row estimands; and every reading again without separation scores of 0.995 or
more.

**The held-out score, framing nothing.** Once the guard's reading is committed, the prediction rows' outcomes are
unsealed and each fit's held-out log score and Brier score are read; the report gives, per dataset, the mean
difference move minus shipped with its Welch interval, on every R2 dataset, those whose references failed
included, where it is the only real-outcome reading. The question is whether two samplers of one posterior
predict differently on one fixed split, so the fits are its replications and no correction for repeated splits
applies.

**The margin.** One unit, the error of 100 independent posterior draws, is a Monte Carlo precision target: at
equal time the move may err by at most what 100 independent draws would add, an error under a tenth of a posterior
standard deviation on a mean, small beside the posterior's own spread, as the house's bars are absolute (dec-B455).
It is not matched to the house's 0.02 coverage bar, and its coverage equivalents differ by estimand, for a normal
posterior and its 90 percent interval: on the mean, one unit displaces the interval by 0.1 standard deviations and
leaves coverage at 0.898 (a 0.02 loss takes 0.35, about 12 units); on the standard deviation, one unit is about 7
percent too narrow, coverage 0.8735, a loss of 0.027 (0.02 is about 0.6 units); on the interval ends, one unit is a
shift of 0.21 standard deviations, and both ends shifted inward cover 0.848, a loss of 0.052 (0.02 is about 0.18
units). So on the ends the margin is about 2.5 times looser than the coverage bar. These treat the unit as all
bias; loss of coverage is quadratic in a mean's displacement, so the mean's 0.898 holds however the error splits
between bias and noise, but on the standard deviation and the ends a unit of zero-mean noise costs much less
(0.8977 and 0.8924), and the ends' 0.848 is the worst case, all bias and inward. Calibration itself is read on
step 2's drawn truths (call 14).

**Power.** A fit's error on a row estimand spreads over fits with standard deviation s in these units, so with n
fits a library D's standard error is s (2 / n)^(1/2), 0.316 s at 20 and 0.224 s at 40. At a true difference of 0
one test resolves with probability 0.99 at s = 0.71, 0.91 at 0.94 and 0.59 at 1.41 at 20 fits; 0.997 at 0.95 and
0.88 at 1.41 at 40; 0.994 at 1.41 at 80. The guard is a conjunction of 15 tests a dataset, correlated within a
dataset: at 5 effectively independent tests a dataset it needs about 0.993 to 0.9965 a test, which 20 fits give
only up to s of about 0.67 to 0.71, so every dataset takes 40 fits from the start. Then a dataset's own conjunction
resolves at s = 0.95 with probability about 0.98. Interval ends rest on tail effective draws, have the largest s
and will usually decide it. On a mean, s is at most 141 / E when the rows' errors move together and less when they
do not. R2 spans about 120 to 4,000 fitted rows, so one dataset's s says little about another's, and each dataset
is its own internal pilot: after its 40 fits a library, read after step 3 is frozen (below), each library's spread
and skewness of the error over its own fits are read, never their difference, and if any row estimand's s passes
0.95 in any grouping, that dataset takes 80 fits a library, about 0.2 core-hour more. At 80 a test resolves 0.994
at s = 1.41. The pilot is blind to D only in that it computes no difference: for errors whose spread tracks their
mean, each library's spread nearly reveals its mean, which is harmless only because the rule is mechanical and its
sole output is the fit count. Re-sizing a dataset on its own fits moves its test's level by at most a few
thousandths (the referee's simulation: 0.024 to 0.030 against 0.025, worst with the threshold at the true s). A
harm of twice the margin shows with probability near 1 at s up to 1. Over datasets, each dataset's D carries at
most about 0.22 units of its own error at 40 or 80 fits, so at a true difference of 0 on every dataset the reading
over datasets is limited by the count accepted, not by the fits; a spread of true differences across datasets
around 0 is what would widen it ([Sizing the real-data stages](#sizing-the-real-data-stages)).

**Expected reading, written before any fit.** The chance that the guard reads harmless at a true difference of 0,
averaged over the uncertainty in the acceptance rates, is about 0.83 at twelve datasets and 0.75 at
ten; the table and its inputs are in [Sizing the real-data stages](#sizing-the-real-data-stages). Where tail
effective draws are few, s passes 0.95, the dataset goes to 80 and the figures hold about as shown.

**Shared with step 3.** The reference, its acceptance and the default fits (20 a library in step 3, 40 here, of
which step 3 reads the first 20) are step 3's design for R3, so R2's runs give step 3's reading on its twelve
datasets with no further compute. They join step 3's test only when the whole of step 3 is committed before any R2
reading is opened: the selected trigger itself, its feature and threshold, computed from step 2's simulated fits by
step 3's rule, which is mechanical and leaves no discretion, so knowing R2's and R3's blocks from the committed
assignment cannot steer it; the need-the-move rule; the refutation rule; and R3's size. Step 2's report waits for
that commit whether or not step 3 ever runs, since it shows R2's reference mass, step 3's own reading, before the
maintainer decides whether step 3 runs. Any later amendment of step 3 drops R2 from its test. Step 3 then tests on
every passer the frame yields.

**Storage.** Full draws are kept until every reading that needs joint draws (the band's mass, the average
probability, the Wasserstein distance) has been scored, about 0.6 GB a reference for the latent fit at 5,000
rows, then compacted as the k study compacts its references.

**What dropping the held-out score as the guard loses** (call 12). A held-out score needs no reference, so it reads
datasets whose references fail, the hard ones, where the samplers differ most, in the units of predictive
accuracy users judge models by, comparable with the September study and with Chipman, George and McCulloch
(2010); it is kept as a report for that reason. As the guard it would be weaker: a no-harm guard asks about fits
that do not need the move, where references exist; held-out probabilities under logistic converge where k does
not ([binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md)), so the score was expected to show
nothing either way; and it barely sees the interval ends and spread, which the reference reading scores.

## The convergence measure

**Question.** From the start a user gets, how many sweeps, and how many seconds, does each sampler need before
its chains draw from the stationary distribution, the posterior? Step 2's arms keep bart()'s fixed 500 burn-in;
this measure runs each sampler with no burn-in and reads the whole stream from its start (dec-B622). It reads
how much of the harm is the start, which the k study left open, and whether 500 sweeps suit either sampler
under logistic. dec-B622 keeps the fixed burn-in as the default, so the plan proposes no change to it. The answer
is read as curves over a series of run lengths, in sweeps and in seconds, with no convergence bar.

**Its readings are published ones (dec-B633).** Each reading below is a form that published comparisons of samplers
use, as the survey found them ([sampler-convergence-measures.md](../design/sampler-convergence-measures.md)):
- the coverage and length of intervals and the error of what users estimate, at fixed run lengths from different
  starts, with the running time (He and Hahn 2023);
- R-hat as a trajectory, in its rank-normalized split form, with bulk and tail effective draws (Brooks and Gelman
  1998; Vehtari, Gelman, Simpson, Carpenter and Burkner 2021), and effective draws per second (Polson, Scott and
  Windle 2013; Pratola 2016);
- the error of estimates against a trusted reference as the run lengthens (Margossian and coauthors 2024).

None of these comparisons gates a conclusion on a convergence bar, and none is set here. No coupling of a BART
sampler exists, so the survey's coupling bounds are not read. A choice below with no published model is named as
the study's own.

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
  that reaches small k from above, as from the default start, reaches it as fast from either side, and it is the
  guard Reading 2 needs against chains that agree because they started together. On the first 25 replications of
  each regime.

**Streams and fits.** On each replication, 20 independent chains of each library from the start, as one stream
with no burn-in, every other setting bart()'s default (75 trees, chi(1.5, 2), `n.threads = 1`). The chains are
read as the fits a user runs: five fits of bart()'s 4 chains, or in the effect regimes two fits of 10, as
bartCause runs them. A fit at length L holds every draw of its chains from sweep 1 to sweep L, nothing
discarded, so every reading is of the whole stream from its start (dec-B622). The replications are step 2's own,
the first 50 of each regime in its order, then up to 100: the same data, truths and row labels, so no selection
runs again and every label comes from step 2's fits, never from the truth. Each chain's draws of the quantities
read are held for the stream, each fit is summarized at each length, and only the summaries are kept. A library's
streams in a regime run to 2,000 sweeps, and on to 8,000 where, at the first look, k's R-hat over a replication's
20 chains at 2,000 sweeps is at or above 1.01 in more than half the replications (the rule, a choice of how long
to run, is the study's own). Past 8,000 a curve is reported as it stands, beside the k study's long runs, which
at 5,000 sweeps of burn-in and 20,000 recorded did not reach agreement between the two samplers in most logistic
fits with a small true k.

**Run lengths.** Every reading is taken at L = 25, 50, 100, 200, 300, 500, 700, 1,000, 1,500 and 2,000 sweeps, and
at 3,000, 4,000, 6,000 and 8,000 where a stream runs on. The steps grow by factors of 1.4 to 2, so the curves are
sampled about evenly on a logarithmic axis, as the run lengths of Brooks and Gelman's batches and of Margossian and
coauthors' figures are. The series holds 100, the length of He and Hahn's started chains; 500, bart()'s burn-in;
and 1,000, bart()'s whole default run of 500 burn-in and 500 draws, so the fit at 1,000 with nothing discarded
stands beside step 2's default fit, which discards 500. It starts at 25, where a fit of 4 chains holds 100 draws,
the fewest from which the ends of a 90 percent equal-tailed interval are read. The points themselves are the
study's own. Each length is also read in seconds: L times the library's seconds per sweep in that regime, from
step 2's timing on the idle box, the move's whole cost charged, with each stream's own clock recorded beside it.

**Quantities.** Every class step 2 scores, labeled as step 2 labels it: k; the latent fit at the scored rows,
interior and near 0 or 1; the band on probabilities, from 250 sweeps on, where a fit of 4 chains holds 1,000
draws (the cutoff is the study's own: shorter fits hold too few draws for a simultaneous interval); the log risk
ratio in S2; individual effects, subgroup effects, the average effects overall and on the treated, the log risk
ratio and the log odds ratio in S7 and S8.

**Reading 1, what users estimate at each run length (He and Hahn 2023).** For each library, start, class and
length, each fit's 90 percent equal-tailed interval, formed as step 2 forms it: its empirical coverage of the
simulated truth against nominal, pooled over rows and per fit as step 2 weighs them, with its interval and the
replications behind it; its length; the root-mean-square error and bias of the fit's posterior mean against the
truth; and the running time at that length. This is He and Hahn's table 4, read at every length of the series
rather than at their two, and from the default start and from draws from the prior, where they compared the
default start with chains started at draws of an approximation. Under step 2's truths, drawn from the model's
prior, a fit that draws from the posterior covers at nominal on average over replications, whatever the data
(step 1's argument), so coverage under nominal at a length reads a fit that has not arrived, or one too short for
its intervals' ends: Reading 2's effective draws say which. The curves are drawn against sweeps and against
seconds, both libraries on one axis, as Biswas, Jacob and Vanetti (2019, figure 4) set two samplers side by side.

**Reading 2, agreement and efficiency.** Rank-normalized split R-hat (Vehtari et al. 2021) at each length, for k
and each scored quantity: over each fit's chains, the diagnostic its user would read at that length; over a
start's 20 chains; and over both starts' 40 where both ran. It is plotted against sweeps and seconds with its
pooled and within-chain variances, as Brooks and Gelman (1998) plot it. They compute each point on the second half
of the run so far; here each point is computed over the whole fit at that length, so that it is the R-hat that fit
would report with no burn-in (dec-B622), a departure that is the study's own. The first length from which R-hat
stays below 1.01, Vehtari and coauthors' threshold, is reported in sweeps and in seconds. Beside it, each fit's
bulk and tail effective draws at each length, against Vehtari and coauthors' 400, and its effective draws per
second, the efficiency reading of Polson, Scott and Windle (2013) and Pratola (2016); step 2's reading of
effective draws per second after burn-in gives the same once the chains have arrived. R-hat reads agreement, not
arrival. Chains that share a start or a slow direction can agree before they arrive: in Brooks and Gelman's
Weibull example R-hat fell under 1.1 while every chain had explored only one part of the space. bart()'s chains
all share the default start, and k's slow drift from its starting value is shared; the draws from the prior are
the dispersed starts that guard against this, and R-hat over both starts' chains reads whether chains from the
default start agree with chains from elsewhere.

**Reading 3, error against the reference, where one exists (Margossian et al. 2024).** The k study stored its
simulated datasets with their truths and, for each long run that passed acceptance (the two samplers' long runs
agreeing, dec-B462), the pooled draws as 1,000 quantiles a row and k's draws in full. Among the regimes read here
these are 2 of S1's 12 long runs, 4 of S2's 12 and 9 of S4's 24, 15 in A, and 3 of 6 in S5, 6 of 6 in S6, 5 of 6
in S7, 4 of 6 in S8 and 8 of 12 in E3, 26 in B. On each, 20 chains of each library from the default start, read as
fits as above, and at each length each fit's squared error against the reference in the posterior mean of k and
of each scored row's latent fit, scaled by the reference's posterior variance, per dataset and averaged over each
regime's datasets, with the error of 100 independent draws marked, as Margossian and coauthors mark it. The same
for the ends of the 90 percent interval is the study's own addition, since the intervals are what step 2 reads.
The reference is missing where the shipped sampler is slowest: the accepted datasets in A are the easier ones, a
logistic long run with a true k below 0.5 having passed 1 time in 16 in the k study, and long runs that failed
are not used, since a reference from the move's half alone is the sampler under test (dec-B456). So Reading 3
checks Reading 1 only on the easier datasets. The references were made on the k study's engine; the
posterior does not depend on the engine as long as main fits the stored data under the same model, which the
harness checks on each dataset before use (the cut points and the prior main builds against those the k study's
harness builds at its recorded commit), dropping a dataset that differs.

**Replications and noise.** Intervals bootstrap the replications, 4,000 resamples, every reading of a replication
resampled together. A fit at a length is a fit of step 2's kind, so planning spreads are step 2's per replication
(Sizing and stops), measured at the first look; the shorter fits' spreads are assumed equal. At 50 replications
the per-fit standard error of the interior's coverage at a length is then about 0.011 to 0.031 at spreads of
0.076 to 0.22, about half that pooled; at 100, 0.008 to 0.022. The five fits of a replication share its data and
truth, so they narrow only the Monte Carlo part. k is one value a replication, its coverage indicator spread 0.30,
so its coverage is known to about 0.04 at 50 and 0.03 at 100; k's reading leans on Readings 2 and 3, and it is a
diagnostic beside the row classes, which are what users estimate. R-hat and effective draws are read per
replication and reported by their distribution over replications.

**Stop.** A first look at 50 replications a regime, both libraries, default start. A regime goes on to 100, He
and Hahn's count of replications, where a row class's standard error of coverage at 500 or 1,000 sweeps exceeds
1.5 times its planned value, as step 2 restates its power; otherwise it stops at 50. No replication is added past
100. The dispersed start stays at 25, its streams lengthened by the same rule.

**How it is read.** Reported beside step 2's readings, framing no recommendation (dec-B633, call 11): the report
gives each reading by regime, start and class.

## Step 3: a check on the data, if the maintainer sends pricing there

Step 3 runs only if the maintainer, having read step 2's report, sends pricing here (dec-B624); step 2's
readings recommend it when the cost is too high.

**Question.** Can a rule that reads only the data find the fits that need the move, so that only they pay for it?
A rule on the data alone keeps the kernel exact. A fit that crosses it through `setData` changes kernel
mid-chain, and both kernels are exact.

**Form, fixed on the simulated fits of step 2 before any R3 dataset is read. The selected trigger itself, with the
rest of step 3 (the need-the-move rule, the refutation rule and R3's size), is committed before
any R2 reading is opened, so that R2's datasets can join the test ([The real-data guard](#the-real-data-guard)).**
The candidates are the separation score at or above t (t from 0.70 to 0.95 by 0.01), the smaller class's count at
or below e (e in 10, 20, 40, 80), since the score misses rare outcomes (the k study's [Could a rule find the fits
that need it?](../design/binary-k-mixing-need.md#could-a-rule-find-the-fits-that-need-it)), and the union of one of
each. A replication's gain is d(shipped at equal time) minus d(move) on the row classes. For each candidate the
gain is pooled over the replications it fires on and, apart, over those it spares, so noise in single replications
cancels rather than being truncated into gain. The trigger is the candidate that fires on the fewest of B's fits
while the upper bound on pooled gain over its spared set, one-sided 0.025 by a bootstrap over replications (4,000
resamples, every cell of a replication resampled together, as step 2's standard errors are), is below 0.02 in every
A cell, an empty spared set meeting the condition. Ties are broken in this order: the single-feature form over the
union; then the candidate firing on more of A's replications, the more sensitive; then a score candidate over a
count candidate; then the smaller t; then the larger e. The selection is mechanical: it reads only step 2's
simulated fits and leaves no discretion. The k study's 0.8 was chosen after its data and is not carried over.

**Test on R3.** For each R3 dataset, every passer not in R2 ([Sizing the real-data
stages](#sizing-the-real-data-stages)): a long run of 8 chains of each library, 5,000 sweeps of burn-in and 16,000
to 20,000 recorded at thin 20, as the k study's references ran, started and checked as the real-data guard's
references are, their length guard and escalation included; then 20 default-length fits of each library. R2's twelve
datasets, whose runs are the same, join the test when the whole of step 3 was committed before any R2 reading was
opened, so the test reads every passer; any later amendment of step 3 drops them from it.

A dataset needs the move when its long runs pass the k study's acceptance (both libraries' halves agree, k's
quantiles included, dec-B462) and then, scored against the pooled reference, 20 shipped fits at equal time hold
a reference mass on the row classes at least 0.02 from 0.90 (its t bound excluding 0), while 20 move fits at the
default length are within 0.02. It does not need the move when the shipped fits' upper bound is below 0.02. This
reads the default-length harm users meet, not long-run disagreement, which the trees' mixing drives as much as
k's (dec-B555). Datasets whose references fail acceptance, or that read neither way, are unresolved and counted
apart. Expected, from the k study: most high-score datasets fail acceptance, as 10 of its 12 did. The trigger is
computed on each dataset's training rows.

**Reading, fixed now.** The real-data test can refute a trigger and cannot confirm one: over the whole frame about
2 datasets are expected to resolve as needing the move ([Sizing the real-data
stages](#sizing-the-real-data-stages)). The trigger is refuted when it spares a dataset that needs the move, the
miss confirmed as the guard confirms a harm (20 more fits of each library on fresh seeds named before the rerun,
against the same reference, reading the dataset as needing again), or when the exact one-sided 0.025 upper bound
on its specificity, the share of the datasets not needing the move that it spares, is below one half: below half,
a second kernel would buy back too little of the cost to be worth carrying. Otherwise it is not refuted. Its
sensitivity and specificity on real data are reported at every count with their exact binomial intervals, with no
floor on either count, and the report says that sensitivity cannot be confirmed on real data at this size. The
evidence for sensitivity is step 2's simulated regimes, where needing cases are many: the report gives the
trigger's sensitivity and specificity there with their intervals, read on the fits that selected it and so labeled
(call 34).

**After the reading.** Step 3's results go to the maintainer as a research report ([The report](#the-report)),
the trigger's firing and sparing shown dataset by dataset, and the maintainer decides (dec-B624). A trigger
not refuted recommends the move on where it fires, resting on step 2's sensitivity and the real data's failure to
refute it; a refuted trigger recommends step 4.

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
| step 2, R2 references | the posterior each real-data fit is scored against | 12 datasets, 16 chains of 21,000 to 25,000 sweeps, the shipped half to 100,000 where only effective draws fail | about 8.0 with the overhead and two escalations; 17.6 at most | accepted or not, by the k study's acceptance and length guard and posteriordb's draw count | a failed reference leaves its dataset unresolved | 12, assurance 0.83 that the guard reads harmless at a true 0 | none |
| step 2, R2 scored fits and timing | no harm on real outcomes | 40 fits a library a dataset, 80 where its own s passes 0.95, as many more where a harm shows, nine timing pairs and nine one-sweep pairs | about 3.6; 22.8 at most | D's se 0.224 s at 40, 0.158 s at 80; s unmeasured | 1, the error of 100 independent draws | 40 at s up to 0.95, else 80 | each dataset its own internal pilot |
| step 2, timing | r, and the cost bar | 9 regimes, 8 large cells, hazard and counts cells | about 1 as 100-sweep segments (3.5 as full fits), box idle | ratio to about 0.01 | 1.25 on a cell | 9 pairs; 10 segments | none |
| step 2, effective draws a second | reported | 10 a regime, S1, S2, S4, S6, both libraries | 0.3 | | decides nothing | | none |
| convergence, default start, Reading 1 | coverage, length, error and time of a user's fit at each run length with no burn-in, both libraries, 8 regimes | 50 a regime (50, 100), 20 chains as five fits of 4 (two of 10 for effects), 2,000 sweeps, 8,000 where k's R-hat stays at or above 1.01 | 3.8 at 2,000 sweeps; 5.2 with the shipped sampler to 8,000 in A; 18 with every second look and every shipped stream to 8,000 | per fit, interior coverage se 0.011 to 0.031 at 50, 0.008 to 0.022 at 100; k 0.04 and 0.03 | decides nothing; reported beside step 2 (dec-B633) | 50; 100, He and Hahn's count, where a first look's se exceeds 1.5 times planned | first look, unless a standard error exceeds 1.5 times planned |
| convergence, dispersed start | does a sampler arrive as fast from a prior draw | 25 a regime, 20 chains, as above | 1.9; 2.6 with A's shipped streams to 8,000; 4.5 at most | about 1.4 times the default start's at 50 | read for its direction | | none |
| convergence, Reading 3 | error against the reference as the run lengthens, on accepted long runs | 41 stored datasets, 20 chains a library | 0.4; 0.8 at most | per dataset | | | none |
| convergence, Reading 2 | what R-hat and effective draws say at each length | the same streams | 0 | | frames nothing | | none |
| step 3, R3 long runs | which datasets need the move | every passer not in R2, about 23 to 33, 8 chains a library, the shipped half to 100,000 where only effective draws fail | about 19 with the overhead and five escalations; 41 at most | per dataset | | the whole frame, measured; a refutation test, its error rates stated | none |
| step 3, R3 default fits | the default-length harm on each | 20 shipped at equal time and 20 move fits a dataset, as many more to confirm a miss | about 2.3 | t bound a dataset | 0.02 | 20 | none |
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

**Audit of the convergence measure.** Its streams cost 20 chains times the stream's length at the k study's
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
  for the shipped sampler, so its streams are expected to stop at 2,000 in most regimes. The cap is 8,000;
  20,000, the k study's length, would add about 3 core-hours expected and up to 20, and on the hardest datasets the
  k study's long runs already read past it.
- Replications. 50 give each sampler's coverage at a length to the standard errors in the table, which shows a
  difference the size of the k study's equal-time readings in A (shipped 0.78 to 0.85 per fit against the move's
  0.89 to 0.91) at about 3 standard errors or more. He and Hahn read 100; the second look reaches it where the
  assumed spreads prove too small. The dispersed start at 25 is read for its direction, its intervals about 1.4
  times the default start's at 50.
- Chains. 20 a replication, read as five fits of bart()'s 4 chains (two of bartCause's 10 in the effect
  regimes), so that every reading is of a fit a user runs. The 19 of the earlier design served a reading dec-B633
  dropped; the twentieth chain adds 5 percent to the measure's compute. Fewer fits a replication would cost less
  and lose little, since the fits share their replication's data and truth, but R-hat over 20 chains is the more
  sensitive of the two.
- S3 is left out (about 6 core-hours at 2,000 sweeps, 15 with the shipped sampler's longer streams).

### Sizing the real-data stages

Each real-data stage is sized by a sample-size justification of the kinds Lakens (2022) sets out. R2's is an
a-priori power analysis whose uncertain inputs are carried by averaging the probability of the primary reading's
conclusion over their uncertainty at a fixed true difference of 0: an expected power over nuisance inputs, closest
to the expected power Kunzmann and coauthors (2020) review, which averages power over the effect's prior. Its target
is 0.8, which Kunzmann and coauthors give, with 0.9, as the common choice for power in trials; Lakens writes that
"the default recommendation to aim for 80% power lacks a solid justification" and that weighing costs "would
probably make studies designed to achieve a power of 90% or 95% more common for articles reporting a single study"
(call 27). R3's is the frame measured whole, the first of Lakens's justifications, with a resource constraint the
reason the frame is not wider. Simulated, 8,000 draws a size, a Monte Carlo error of about 0.005.

**Inputs.** Acceptance from the k study's record of its 22 real prediction references, by separation score and
size, each cell under its Jeffreys posterior: at a score of 0.8 or more, 0 of 5 datasets of at least 1,000 rows
accepted and 2 of 7 smaller; below 0.8, 1 of 2 and 7 of 8. The share scoring 0.8 or more by size, from the same
record: 5 of 7 large, 7 of 15 small. Every one of those failures failed row R-hat, which this acceptance keeps, so
none would pass here. The pool's sizes are known: 14 of its 24 have at least 1,000 rows, against 7 of the k study's
22, which is why acceptance by size matters. Agreement's chance failure is counted again on large datasets, a loss
of 0.15 of their accepted references, though the record's rates already hold part of it at their own sizes, so the
inputs are conservative there (call 28). Each accepted dataset's spread s is uniform from 0.5 to 1.41 on each of
five effectively independent row estimands and classes, 40 fits a library and 80 where any s passes 0.95; R2 is
drawn by blocks from the pool of 24.

**R2.** The primary reading is the guard reading harmless ([The real-data guard](#the-real-data-guard), Reading
over datasets: every accepted dataset on its own bound and the t over them) when the true difference is 0 on every
dataset. Given 3 or more accepted references it reads so nearly always, given 2 about 5 times in 100, so the
assurance is close to the chance of at least 3 accepted: a feasibility calculation for the count of usable
datasets, the test's power given that count hardly entering. With the sensitivity analyses Lakens asks for beside
it:

| R2 | assurance | without the extra agreement loss | t alone | references accepted, expected | at least 3 accepted | power at the inputs' 0.2 quantile | true differences spread 0.25 by estimand | spread 0.25 by dataset | acceptance by score alone (the earlier inputs) |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.75 | 0.78 | 0.75 | 3.6 | 0.76 | 0.63 | 0.46 | 0.53 | 0.85 |
| 12 | 0.83 | 0.84 | 0.84 | 4.2 | 0.84 | 0.74 | 0.53 | 0.61 | 0.89 |
| 14 | 0.90 | 0.92 | 0.92 | 5.0 | 0.92 | 0.85 | 0.58 | 0.67 | 0.94 |
| 16 | 0.93 | 0.93 | 0.95 | 5.7 | 0.95 | 0.89 | 0.60 | 0.67 | 0.95 |
| 24 | 0.96 | 0.97 | 0.99 | 8.4 | 0.99 | 0.96 | 0.55 | 0.63 | 0.97 |

R2 is 12, the smallest size whose assurance reaches 0.8; at 0.9 it would be 14 (call 27). Its power at a
pessimistic draw of the inputs, the 0.2 quantile of power over their uncertainty, is 0.74 (Kunzmann and coauthors'
quantile approach evaluates power at a quantile of one uncertain input; with several, a quantile of the power stands
in for it). Where the samplers' true differences spread across datasets, the between-dataset term of Boulesteix and
coauthors' frame, no size reaches 0.8 under the primary reading: every added dataset is one more that must read
harmless alone, so the assurance peaks near 0.6 at 16 and falls after. That is the cost of asking every dataset to
be harmless; it is reported and not sized on (call 29). Where the pool scores higher than the record's datasets
did, the guard is also more often unresolved.

**R3.** Every passer not in R2, the census read to its end under both screens if step 3 runs, before any R3 file is
loaded: step 3's test then reads every passer of the frame, so the frame is measured whole, and the frame is what
the published sources this rule reads hold. At the census's yield so far, 3 to 5 passers in its first 80 draws, the
frame holds about 35 to 45 passers, so R3 is about 23 to 33 datasets beside R2's 12. With R2's inputs, the
agreement loss carried as there, and further: an accepted dataset scoring 0.8 or more needs the move with
probability uniform on 0 to 1 (the k study accepted two such references), one below with probability drawn from
Beta(0.5, 9.5), mean 0.05 (below 0.8 the k study's shipped interior was within 0.004); a needing dataset resolves
as needing 9 times in 10 at 20 fits, another as not needing 8 times in 10. Over the whole frame about 1.9 to 2.4
datasets are expected to resolve as needing the move, none at all about one time in four, and about 8.5 to 10.7 as
not needing: the datasets that need the move are the ones whose references most often fail.

So the real-data test is built to refute a trigger, not to confirm one (Step 3, Reading). A confirmatory pass would
need far more: the earlier pass rule (every needing dataset caught, at least 3 resolved on each side, half the
others spared) reads with assurance about 0.31 on 40 datasets and reaches 0.8 only at about 150, and a bound on
sensitivity against a stated value is sized by the count of datasets with the condition over its prevalence
(Buderer 1996, as Bujang and Adnan 2016 give and tabulate it), hundreds at a prevalence near 0.05. The 0.8
convention would not carry over either: Kunzmann and coauthors warn that "conventions on power levels cannot be
transferred to thresholds for probability of success without adjustment" for the prior chance of a relevant effect,
and here the chance that any dataset needs the move is itself uncertain, with real mass near 0. Lakens's points for
a size set by resources, answered:
- the decision: you decide on step 3's report whatever it shows, and no later study is planned to pool with it;
- the intervals expected: with both needing datasets caught, sensitivity's exact 95 percent interval runs about
  0.16 to 1; with 7 of 9 others spared, specificity's about 0.40 to 0.97;
- the refutation rule's error rates at 36 and 45 passers: a trigger that is right (it fires on every needing
  dataset and spares three in four of the others) is refuted essentially never; one that misses one needing
  dataset in ten is refuted 0.16 to 0.19 of the time, one that misses half 0.54 to 0.60, and one that spares only
  three in ten of the others 0.12 to 0.16. A bad trigger escapes about half the time, and the report says so;
- the evidence for sensitivity: step 2's simulated regimes, where needing cases are many (call 34).

**Reserves.** None, by order rather than by count. The pre-fit check, where a file fails (its counts, its parse,
its checksum), runs on every passer before the draw, and a shortfall it leaves is met by reading the census deeper
under both screens, so no reserve is needed for it. Each file is archived privately, outside the tree, at its first
download, as the licenses allow; a file that fails to reload after the draw is restored from the archive and
checked against its recorded sha256, or fetched again from its recorded address, and one that cannot be is
reported unresolved, as a failed reference is, and not replaced, so no replacement selects on anything after the
draw (call 32). The chance that a drawn dataset ends unresolved for a lost file is taken as negligible and is not
in the sizing. The documentation screens found two entries whose files are not the described data (spectf and
cylinder-bands); the pre-fit check's own rate is unknown until it runs, which is why the census meets it rather
than a fixed count.

**Sources.**
- Boulesteix, A.-L., Hable, R., Lauer, S. and Eugster, M. J. A. (2015). A statistical framework for hypothesis
  testing in real data comparison studies. The American Statistician; read as Technical Report 136, Department of
  Statistics, University of Munich, 2013.
- Buderer, N. M. F. (1996). Statistical methodology: I. Incorporating the prevalence of disease into the sample
  size calculation for sensitivity and specificity. Academic Emergency Medicine 3(9), 895-900; doi
  10.1111/j.1553-2712.1996.tb03538.x. Not read: cited as Bujang and Adnan give it.
- Bujang, M. A. and Adnan, T. H. (2016). Requirements for minimum sample size for sensitivity and specificity
  analysis. Journal of Clinical and Diagnostic Research 10(10), YE01; doi 10.7860/JCDR/2016/18129.8744.
- Demsar, J. (2006). Statistical comparisons of classifiers over multiple data sets. Journal of Machine Learning
  Research 7, 1-30.
- Kunzmann, K., Grayling, M. J., Lee, K. M., Robertson, D. S., Rufibach, K. and Wason, J. M. S. (2020). A review of
  Bayesian perspectives on sample size derivation for confirmatory trials. arXiv:2006.15715.
- Lakens, D. (2022). Sample size justification. Collabra: Psychology 8(1), 33267; doi 10.1525/collabra.33267.

**Audit of the real-data guard.** Its costs use the plan's own rate for R3: 12 references of 16 chains at 25,000
sweeps in about 5 core-hours, about 3.75 milliseconds a chain-sweep at up to 4,000 fitted rows, and 40 default fits
in about 0.08 core-hours; every figure below carries the plan's 35 percent for scoring and shared memory bandwidth.
A reference costs about 0.56 core-hours at 25,000 sweeps and 0.47 if it is accepted at 21,000. An escalation of the
shipped half to 100,000 sweeps adds 8 chains of 80,000, about 0.9; it runs only where effective draws are the only
failing checks, the rarer failure in the k study (most of its failures were on row R-hat and on agreement), so
about two are expected in R2 and five in R3, and every reference at most. R2's twelve references cost about 8.0
core-hours expected, 17.6 at most. The scored fits at 40 a library, about two datasets raised to 80, and the timing
fits cost about 3.6; at most, with every dataset at 80, a confirmation of 80 on every dataset and every fit at the
reference's rate, 22.8. Together about 11.6 core-hours expected, 40.4 at most, against the held-out score's 1.7.
R3's references, at about 28 datasets, cost about 19 expected, 41 at most, and its default fits about 2.3; a
confirmation of a miss about 0.1. On one of the box's 2-core legs a reference takes about 13 minutes of wall time,
an escalation about 25 more. If step 3 runs, R2's runs serve it at no cost. Memory: about 0.6 GB a reference before
compaction. Power is in [The real-data guard](#the-real-data-guard); the one figure not measured is the spread s,
which each dataset's own fits measure.

**Totals.** About 49 core-hours expected through step 2 with the convergence measure: step 0 about 5, step 1 about
7, step 2 about 28.3 (about 21 if its regimes stop at their first look), the convergence measure about 8. About 70
with step 3, on every passer the frame yields, its references escalating as the guard's do (about 58 if the census
is not read past the sealed pool). About 152 at most, with step 1 at its longest thin and its third arm, the timing
run as full fits, the real-data guard at its upper cost, a confirmation of one harm in step 2 (about 4, S8 at its
full count), the convergence measure at every second look with every shipped stream at 8,000, and every R2 and R3
reference escalating (about 36 of the 148, though an escalation runs only where effective draws are a reference's
only failing check); about 15 more if the move also needs its long streams everywhere. The real-data guard adds
about 8 expected and 32 at most over the held-out score it replaces (dec-B651). Against the sizes before dec-B667
(6 and 12), R2 at 12 adds about 5.8 expected and 20.2 at most, and R3 on the whole frame about 12.2 and 24.8, the
second only if step 3 runs. That is within dec-B530's rough guess of tens of core-hours. Against the amendment
before dec-B621 (37 expected, 43 with step 3, 56 at most), step 2's fitting fell by a third, 25 to about 17.5, with
the burn-in-only arm dropped (dec-B621), and the convergence measure adds about 8 expected and 23 at most
(dec-B622); its rebuild from published readings (dec-B633) keeps its streams and changes only the twentieth chain,
5 percent of its compute (8.2 against 7.8 expected, 23.3 against 22.1 at most). The stop is 105, 1.5 times the
expected total with step 3. Under it, cuts go in this order: the convergence measure's second looks, the dispersed
start's long streams, the contingent arm's thin to its pilot's T, S3 to its first look, the other step 2 regimes to
their first looks, R3 to the sealed pool's passers (the census not read further). R2 is not cut: its size is its
justification. A smaller design, B at 100 a regime and S2 at 200 (about 7 core-hours less), would resolve B's
classes at a true 0 only about 50 to 60 percent of the time and S2's interior per fit 25 percent, so more classes
would reach the maintainer unresolved.

Where it runs: every timing, every fit of steps 1 to 3 and every stream of the convergence measure on the bench
box, at most four jobs at once, the two libraries' streams interleaved, the box idle for timing
([The x86 leg](../../benchmarks/README.md#the-x86-leg)); step 0's gates on the laptop, with its x86 leg on the
box.

## The report

Whatever a step shows, its results go to the maintainer as a research report, and the maintainer decides
(dec-B624). Step 2's report also carries step 1's results and the convergence measure; a step 1 that fails
reaches the maintainer on its own, and step 3, if it runs, has its own report. Each is to the standard for
changing a default ([Studies](README.md#studies)) (dec-B559): the maintainer set that standard for default priors
and said mixing might need it too, and a change of kernel for every logistic fit is a default users meet. Each
report is a short PDF, at most 5 pages of text with graphs as needed, self contained. Every report gives each
sampler's empirical coverage against nominal, with its interval and replications, widths, root-mean-square error
and bias, for every class and band, and the difference between the samplers with its interval (dec-B625). Past
that, its contents are taken from published comparisons of MCMC samplers, read when the report is written:
their readings, tables, figures and diagnostics. The report follows them and cites them, and names as the
study's own anything it shows with no published model (dec-B626). The survey notes which published comparisons
can serve ([Models for the report](../design/sampler-convergence-measures.md#models-for-the-report)); for the
convergence measure and step 2's readings of k they give these, each tied to its model:
- a table by regime of each sampler's coverage, interval length, root-mean-square error and running time, from
  each start, at 500 and 1,000 sweeps with no burn-in and at the stream's end, the form of He and Hahn's (2023)
  table 4;
- curves of coverage, length and error against sweeps and against seconds, both samplers on one axis, as Biswas,
  Jacob and Vanetti (2019, figure 4) set two samplers' readings against iteration;
- R-hat with its pooled and within-chain variances against the run length (Brooks and Gelman 1998, figure 1), from
  both starts, beside autocorrelations, trace plots of k and histograms of k against the reference for the two
  samplers (van Dyk and Meng 2001, figures 3 to 6), the first length from which R-hat stays below 1.01 marked;
- bulk and tail effective draws and local efficiency plots (Vehtari et al. 2021), and a table of effective draws
  and effective draws per second by regime and sampler (Polson, Scott and Windle 2013; Pratola 2016, table 2);
- squared error against the reference as the run lengthens, the error of 100 independent draws marked
  (Margossian et al. 2024, figure 1), on the stored datasets with accepted long runs, the datasets without one
  named;
- k's rank-normalized split R-hat and bulk and tail effective draws from step 2's move fits at the default length,
  against R-hat below 1.01 and 400 effective draws (Vehtari et al. 2021), the September study's first condition
  (dec-B634).

It says how the regimes and datasets were chosen, and what 0.9-34 did. Two independent referees read it on a
neutral brief before it reaches the maintainer. Whatever the maintainer decides, the finding is then a design note
in the usual form.

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
- **Held-out scores as a gate, with the datasets named.** Adopted, then amended (dec-B651): the guard on real
  outcomes reads each sampler's error at equal time against a reference run of the same posterior in place of a
  held-out score, a no-harm guard and not a gain criterion
  ([The real-data guard](#the-real-data-guard)). It runs on R2, twelve datasets none of the probes has seen, drawn
  by the rule in [k-move-pricing-datasets.md](k-move-pricing-datasets.md); about 11.6 core-hours, 40.4 at most. The
  held-out log score and Brier score are reported beside it, framing nothing.
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
- **Call 6, where the move waits.** Settled by the standing rule that a feature merges complete (dec-B494): (a),
  the move stays on its rebuilt branch until pricing decides and lands once in the form decided, the bitwise
  check at landing applying.
- **Call 8, the September prior study's revisit under logistic.** (a) (dec-B634): step 2 reads k's
  rank-normalized split R-hat and bulk and tail effective draws from its move fits at the default length, against
  R-hat below 1.01 and at least 400 effective draws (Vehtari et al. 2021), in place of the earlier 1.05 and "in the
  hundreds"; if they hold and the move lands for logistic, a logistic rerun of the study goes to the maintainer as
  its own item after the landing. The probit revisit's condition under dec-B390 is a separate question.
- **Call 9, negative binomial's item.** The orchestrator's internal call, (a): TODO k-mixing-pg-families, which
  already holds the Polya-Gamma families' k mixing, is negative binomial's item.
- **Call 11, how the convergence measure is read.** (a) (dec-B633): reported beside step 2's readings, framing no
  recommendation; the measure is rebuilt from published comparisons' readings, and the study's own reading of the
  truth's rank among chains and its 0.02 bar are dropped.

Open.

**Call 7. The held-out datasets.** Open for the rule and the list (dec-B635, dec-B651); the sizes ruled by method
(dec-B667). The rule and the list are in [k-move-pricing-datasets.md](k-move-pricing-datasets.md), with the
per-entry screen in [k-move-pricing-screen.csv](k-move-pricing-screen.csv): a published frame, a screen of six
criteria that bear on a comparison of two samplers of one posterior, 23 passers on documentation after the
census's second screen, 14 with at least 1,000 rows, the census read on until 24 pass, and a seeded stratified
draw of R2 and R3 fixed in it. After the survey, the two audits and the first referee, a planner rebuilt how the
datasets enter; independent referees read the rule and list, and they were revised on their reports (in
[Amendments](#amendments)). No 0.9-34 study used held-out real datasets for a sampler change.

The sizes are not a choice put to you: you ruled that they be set as a research study sets them (dec-B667). Each is
set by a sample-size justification, in [Sizing the real-data stages](#sizing-the-real-data-stages), and a referee
read the first sizing and its revision answers that report. R2 is 12, the smallest size at which the guard reads
harmless at a true difference of 0 with probability 0.8 averaged over the uncertain inputs (0.83), acceptance read
from the k study's record by score and size. The guard reads harmless only when every accepted dataset is harmless
on its own bound and a one-sample t over them is below the margin, the frame Boulesteix and coauthors give
real-data comparisons, in place of the floor of 3 accepted references. R3 is every passer not in R2, the census
read to its end if step 3 runs, so the frame is measured whole; its test refutes a wrong trigger and does not
confirm a right one, since only about 2 datasets in the frame are expected to show a need for the move, and step
3's evidence for sensitivity comes from step 2's simulated regimes. Reserves are replaced by the order of the
pre-fit check and the census, with each file archived at its first download. The sealed rule and list then come
to you as dec-B635's decision.

### Calls made in the dec-B651 amendment and its revision

The planner's, each with its alternative; the dataset rule's own calls are in its section 8.
- **Call 12, the guard's measure.** Each sampler's error at equal time against a reference of the posterior,
  posteriordb's form, in place of the held-out log score as the guard; the held-out log and Brier scores are
  reported beside it on the prediction rows, unsealed after the guard's reading is committed, framing nothing
  (the referee's recommendation). Alternatives: drop the held-out score, leaving datasets whose references fail
  with no real-outcome reading; or make it the guard, as first planned.
- **Call 13, the reference's acceptance.** The k study's as written, with its length guard (the draw count read
  at the accepting block end, the other checks at both) and its escalation of the shipped half, and posteriordb's
  10,000 draws; posteriordb's lag-1 independence and k's R-hat reported, not required; the k study's one-value
  limits gating nothing, since the guard scores no one-value quantity; agreement at 4.5 standard errors on every
  row, whose chance failure at 5,000 rows (about 0.1 to 0.2) is stated. Alternatives: agreement read at the 99th
  percentile of rows, as R-hat is, which removes the chance failures and lets a few rows disagree unseen, or at a
  threshold scaled to each dataset's count of tests; the amendment's first form, R-hat below 1.01 at every row and
  on k with bulk and tail effective draws of at least 1,000, which a stationary reference over 16 chains passes
  reliably only at about 2,000 effective draws and which requires k to converge, selecting the datasets where the
  shipped sampler already mixes; posteriordb's independence literally, by thinning about nine times further for
  k, about 20 core-hours more for R2; or a reference from the move's long runs alone where the shipped half
  disagrees, which would reach the hard datasets and which dec-B462 does not allow.
- **Call 14, the error's unit and the margin.** Error in units of the error of 100 independent posterior draws,
  harmless below 1 on every row estimand, a Monte Carlo precision target and not a match to the 0.02 coverage bar
  (equivalents in [The real-data guard](#the-real-data-guard)); the error ratio reported with its interval, a lower
  bound past 1.25 flagged; k beside. Alternatives: margins per estimand matched to the coverage bar, about 0.6
  units on the standard deviation and 0.2 on the interval ends, at which 40 fits resolve only where s is under
  about 0.3 on the ends, so few datasets would read; a ratio bar of 1.25, which at 20 fits bounds a true ratio of
  1.1 near 1.5 and would seldom resolve; or k inside the guard, which would leave most datasets unresolved on k
  alone.
- **Call 15, the fits.** 40 a library on every dataset from the start, since 20 meet the conjunction only up to s
  of about 0.7; each dataset its own internal pilot, raised to 80 when any s of its own passes 0.95, the
  threshold at which 40 give the conjunction 0.9 over 6 datasets; one look. Alternatives: 20, raised to 40 past
  s = 0.67 and to 80 past 0.95, which saves about 0.6 core-hour where the spreads are small and re-sizes every
  dataset on its own fits; or one pilot dataset sizing all, which reads one dataset's s for datasets of 120 to
  4,000 rows.
- **Call 16, sharing with step 3.** The whole of step 3 (the selected trigger, the need-the-move and refutation
  rules, R3's size) committed before any R2 reading is opened, step 2's report waiting for it, so R2 joins step 3's
  test and the test reads every passer; step 3's references take the guard's acceptance. Alternative: R2 and R3
  kept apart, step 3 testing R3 alone, about 1.3 to 1.8 datasets then expected to resolve as needing the move, step
  2's report not waiting.
- **Call 17, calibration stays on drawn truths.** Coverage of f is read only where the truth is known, step 2's
  nine regimes, which already place drawn truths on fifteen real datasets' covariates. Alternative: a further B
  regime with drawn truths on R2's twelve covariate sets, about 5 core-hours at 200 replications over the twelve,
  or 30 at 200 on each; it would widen the covariates under drawn truths, mostly in the band with no strong signal,
  which B already resolves, and would put R2's covariates to a second use.
- **Call 18, R2 of every size.** The 1,000-row floor kept a test fold of 200 rows, and the guard reads no test
  fold. Alternative: keep the floor, R2 then reaching larger real fits than step 2's regimes, most of which are
  under 1,000 rows.
- **Call 19, the convergence measure.** R2's references are not added to Reading 3. Alternative: 20-chain streams
  to 2,000 sweeps on R2's twelve datasets, about 1.0 core-hour, which would give Reading 3 real outcomes.
- **Call 20, equal time per sweep.** r is the ratio of the two libraries' seconds per sweep, each whole fit less a
  one-sweep fit's time, in step 2's regimes and in the guard alike, and the fixed share reported. A user running
  either fit pays the same fixed time, so equal time gives the shipped arm the sweeps the move's sweeps take; a
  ratio of whole fits pulls r toward 1 and favors the move. Alternatives: r from 100-sweep segments restored from a
  stored state, as the cost benchmark times them, which reads sweeps past burn-in only, where trees are larger; or
  the whole-fit ratio with the fixed share only reported.
- **Call 21, the reading over datasets (dec-B667, revised on its referee).** Harmless only when every accepted
  dataset is harmless on its own Welch bound, an intersection-union test that reads datasets as rule 2 reads
  regimes, and the one-sided 0.025 bound of a one-sample t over the accepted datasets' D is below 1 (Boulesteix
  and coauthors' frame), on every row estimand and class; the exact bound on the share harmed and the Wilcoxon
  signed-ranks test beside. No floor: the count needed follows from the intervals, about 3 at a true 0.
  Alternatives: the t alone, as first revised, which can average a harm on a few datasets away and whose level
  fails where the harm sits in a few (about 0.1 to 0.2 at strongly skewed spreads of true D, the referee's
  simulation); every dataset alone with no claim over datasets, which reads at a true 0 about as often and states
  nothing about datasets like these; the floor of 3 with every accepted dataset harmless, which has no source; or
  the Wilcoxon test as the claim over datasets, which Demsar (2006) recommends over the t and which cannot reach
  one-sided 0.025 below 6 datasets, 6 or more accepted at R2 of 12 with probability about 0.23.
- **Call 22, the reference's halves.** D is also read against each library's half alone: a harmless reading not
  holding against the shipped half is flagged, framing nothing, and a harm stands only when it also holds against
  the move's own half, each the half that can only favor the other reading. Alternatives: harmless required
  against both halves, which halves each reference's draws behind the reading and costs power; or each library's
  fits scored against the other library's half only.
- **Call 23, harms confirmed.** A harm shown, in the guard or in step 2's regimes, frames the recommendation only
  when a rerun on fresh seeds at the count it was read at shows it again, as step 1's failures are confirmed;
  about 0.2 core-hour a guard dataset at 40 fits, a regime's look (0.1 to 4) in step 2. Confirmed, a false harm
  somewhere in A is about 0.0006 against 0.09 unconfirmed. Alternative: no rerun, the report stating the expected
  number of false harms at a true 0 and the family's level, where a harm in A recommends parking.
- **Call 24, joint estimands beside.** The average probability over the prediction rows and the reference mass
  inside each fit's band are reported, framing nothing. Alternative: both in the guard's conjunction, which would
  catch a sampler with the right margins and the wrong dependence across rows, and lower the guard's power.
- **Call 25, where a confirmed harm on R2 sends pricing.** It blocks the move on for every logistic fit and
  recommends step 3, as a harm only in B does: the real datasets whose references pass are mostly fits with no
  strong signal. Alternative: read it as a harm in A, which recommends parking.
- **Call 26, Welch alone.** The guard's bound is Welch's t, its level at the measured skew stated (about 0.035 to
  0.05 at the strongest skew simulated, at 20 to 40 fits). Alternatives: a studentized bootstrap, which helps at
  moderate skew (0.029 at shape 2) and not at strong (0.065 at 0.5); or the percentile bootstrap beside Welch, as
  the revision first had it, which is narrower and so reads as Welch alone.
- **Call 27, the level.** 0.8 on the primary reading, which Kunzmann and coauthors give with 0.9 as common.
  Alternative: 0.9, which Lakens's weighing of costs for a single study favors: R2 at 14 (0.90), about 1.9
  core-hours more expected and 6.7 at most.
- **Call 28, acceptance and agreement's chance failures in the sizing.** Acceptance from the k study's record by
  score and size, and agreement's chance failure counted again as a loss of 0.15 of accepted references on
  datasets of at least 1,000 rows, which double-counts part of what the record holds. Alternatives: the record
  alone, which gives R2 at 12 as well (0.84); or acceptance by score alone, the first sizing's inputs, which gives
  R2 at 10 (0.85) and ignores that 6 of the record's 7 large datasets failed.
- **Call 29, the primary reading's scenario.** A true difference of 0 on every dataset, the guard's question, as a
  non-inferiority design is sized; a spread of true differences across datasets (0.25 units) and a pool scoring
  higher than the record's are sensitivity analyses, reported. Alternative: size on the spread, which no R2 meets
  under the primary reading (about 0.6 at most), or under the t alone (R2 at 14, 0.83 by estimand), which would
  size on a reading the plan does not make.
- **Call 30, R3 measures the frame.** Every passer not in R2, the census read to its end under both screens if
  step 3 runs, before any R3 file is loaded: the whole frame, about 23 to 33 R3 datasets, about 12.2 core-hours
  more than R3 at 12. Alternatives: the sealed pool's other passers only (12), about 12 core-hours less and no
  further census, about 1.3 datasets expected to resolve as needing; or step 3's real-data test dropped, the
  trigger read on step 2's simulated regimes alone.
- **Call 31, step 3's reading.** A refutation test: a confirmed miss, or specificity's exact upper bound below one
  half, refutes the trigger; sensitivity and specificity estimated with exact intervals, no floor; the evidence for
  sensitivity from step 2's regimes. Alternatives: the earlier pass rule with its floor of 3 on each side, which has
  no source and reads with assurance about 0.31 on 40 datasets; or a confirmatory pass on the exact lower bound of
  sensitivity above a stated value (Clopper-Pearson, sized as Bujang and Adnan size it), about 0.18 on 40, which
  the frame cannot power.
- **Call 32, no reserves; files archived.** The pre-fit check runs before the draw and the census meets its
  shortfall; each file is archived privately at its first download, and a file lost after the draw is restored or
  fetched again, else reported unresolved. Alternatives: a fixed count of reserves from the pool, which the
  pre-fit check's unknown rate cannot size; or no archive, the re-fetch depending on the repositories' addresses.
- **Call 33, a donor type that is neither numeric nor factor.** The loader keeps the file's type: numbers stay
  numeric, text becomes a factor (the dataset rule's section 4; it decides no passer). Alternative: read each such
  type as the screener judged it (a date-time as a number, an ordinal as numeric), which is a choice of this study's
  made entry by entry.
- **Call 34, where step 3 reads sensitivity.** On step 2's simulated fits, those that selected the trigger, the
  report labeling it so. Alternative: select the trigger on half of step 2's replications and read its sensitivity
  on the other half, which removes the selection's optimism and halves the data the selection reads.
## Left out

- **Probit**, every probit-only regime (S9, D2, E1, E2) and the probit half (dec-B617).
- **Negative binomial.** Its item is TODO k-mixing-pg-families (call 9). The negative binomial default study
  measured k under a scratch build of that half, and its default is still open.
- **The multiple-treatments designs.** Their miss is shared by both samplers and points to the trees and the model,
  and every design with a real effect failed its first long run.
- **The colleagues' case.** It has not been shared; it has a slot in the k study's plan (dec-B529).
- **The k study's real outcomes.** They carry no truth, and their 22 datasets are spent; real outcomes enter only
  through R2's reference-scored guard and R3's trigger test.
- **A held-out score on real outcomes as the guard** (dec-B651, call 12): replaced by the reference-scored guard,
  the held-out log and Brier scores reported beside it, framing nothing.
- **Coverage on real outcomes**: no truth is known there; calibration is read on drawn truths only (call 17).
- **Count-weighted logistic and hazard.logistic** are reached by the move, gated by its exact gate's counts arm and
  timed, but not scored as regimes (dec-B619): the likelihood is the same, and the counts were not a k study
  regime.
- **S3 in the convergence measure**: its streams at 10,000 rows would cost about 6 to 15 core-hours, and it is in
  neither of the sets the measure reads.
- **A coupling of the two samplers' chains**, the published route to an upper bound on the distance from
  stationarity: none exists for a BART sampler, and one would be new work with its own proof
  ([sampler-convergence-measures.md](../design/sampler-convergence-measures.md)).
- **A change to the default burn-in**: dec-B622 keeps it; the convergence measure reports each sampler's fits at
  500 and 1,000 sweeps with nothing discarded.
- **A convergence bar, and the truth's rank among chains at each sweep**: dropped (dec-B633). No published
  comparison the survey found gates a conclusion on a convergence bar or reads the truth's rank sweep by sweep.
- **Multi-forest fits with amplitudes, BCF and the latent BCF families**: out of the move's reach, as of the probit
  step's.
- **The September prior study's rerun** (TODO binary-hyperprior-revisit, dec-B390): it was a probit study, and with
  the probit half out pricing does not supply the probit mixing fix its revisit condition asks for. It stays its own
  item; step 2 reads its first condition under logistic (dec-B634).
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
- The convergence measure: its review checks each reading against its published form (He and Hahn's table 4,
  Vehtari and coauthors' R-hat and effective draws, Margossian and coauthors' error against a reference) and that
  every fit at a length holds its chains' draws from sweep 1, nothing discarded; a smoke that feeds the scorer
  independent exact draws from a known conjugate posterior in place of the chains, with truths drawn from its
  prior, reads coverage at 0.90, R-hat near 1, effective draws near the draws held and the error near that of as
  many independent draws, at every length within noise; Reading 3's model-equality check passes on each stored
  dataset before use, and the datasets it drops are listed.
- The real-data guard: the dataset rule's seal precedes any file load, its committed loader drops the outcome
  before anything is printed, the assignment and the whole of step 3 are committed before any R2 reading is
  opened, the guard's reading before the prediction rows' outcomes are unsealed, and step 2's report after step
  3's commit (each by commit order); every reference's acceptance figures at both block ends, its escalation and
  its lag-1 report are given per dataset; the internal pilot's record holds each library's spreads and no D; a
  smoke that scores independent exact draws from a known conjugate posterior as both "samplers" and as the
  reference reads D near 0, D against either half near 0 and each error near its 100 / E, and one that scores 100
  independent reference draws reads an error near 1; the equal-time draws are checked against r from the timing
  record.
- `Rscript tools/check-doc-freshness.R .` on this file and the indexes.

## Not verified

- That the rebuild's three semantic changes are as small as read: the two conflicting hunks in chain.hpp and
  the placement before the augmentation draw were read on main, not built.
- The move's cost on main with the redraw dropped; the k study's 1.04 to 1.14 is the old placement's.
- That a thin of 500 or less gives 50 effective draws on each of step 1's four functionals under chi(1.5, 0.1),
  and the selection costs of the tilted arm's S2 rate test.
- The candidate datasets' sizes, rates and files: from suite records and repository pages, checked against the
  files at the pre-fit check; the census entries that complete the pool rest on two screens not yet made.
- The real-data sizing's inputs: the k study's acceptance by score and size and its shares of high scores by
  size, read as this pool's; agreement's chance loss of 0.15 on large datasets; the spread of s over datasets; and,
  for step 3, how often an accepted dataset needs the move (assumed, uniform above a score of 0.8, mean 0.05 below)
  and the census's yield, projected from its first 80 draws to about 35 to 45 passers in the whole frame.
- The real-data guard's spread s over fits, and so its power; the per-sweep times at R2's sizes, taken from the
  plan's R3 figures; how many R2 references pass acceptance, projected from the k study's rates by score and size
  on its own datasets; how often a reference escalates; the share of fixed time in a default fit; the pilot's
  effect on a test's level and Welch's level under skew, from the referee's simulations, not this plan's.
- Per-replication spreads in step 2 between the move and the shipped arm at equal time; the sizes use the k study's
  spreads against its told fit.
- The convergence measure's spreads: a fit's coverage at the short lengths is assumed spread like step 2's
  default-length fits, to be measured at the first look.
- Its probe: shipped sampler only, on the laptop with a library built from main on 2026-11-03, three stored
  datasets; the move's time from the start is not measured, only its autocorrelation at stationarity.
- That main fits the k study's stored datasets under the model their references were made under; the harness
  checks it dataset by dataset.
- k's starting value of 2 at the default start, read from the engine and from the probe's first sweeps, not from
  any documented default.

## Amendments

**2026-11-07, the sizing revised on its referee (dec-B667).** An independent referee reproduced the first
sizing and asked for revision. Taken: the guard reads harmless only when every accepted dataset is harmless on its
own bound and the t over them is below the margin (call 21); acceptance read from the k study's record by score and
size, R2 then 12 (call 28); step 3 a refutation test, its floor of 3 removed and its sensitivity read on step 2's
regimes (calls 31 and 34); R3 justified as the frame measured whole, its figures restated with the agreement loss
and a yield of 35 to 45 passers (call 30); files archived at first download (call 32); Kunzmann and coauthors'
expected power named as the quantity computed, Lakens quoted, Buderer (1996) cited through Bujang and Adnan. In the
dataset rule: a donor type that is neither numeric nor factor keeps the file's type (call 33), Occupancy Detection
then failing 4 first, Bone marrow transplant: children's listing stated, Mammographic Mass passing both screens with
its class-count doubt recorded, so the pool is 24, and the frame count corrected. Costs about 49 core-hours
through step 2, 70 with step 3, 152 at most, stop 105.

**2026-11-07, the real-data sizes set by method (dec-B667).** The maintainer ruled that the sizes be set as a
research study sets them. R2 is 10, the smallest size at which the guard's primary reading reaches its conclusion
at a true difference of 0 with assurance 0.8, the acceptance rates' uncertainty carried (Lakens 2022; assurance as
Kunzmann and coauthors review it); R3 is every passer not in R2, the census read to its end if step 3 runs, a
stated resource constraint, since step 3's pass rule would need about 130 datasets. The floor of 3 accepted
references and the harmed-share bound are replaced by a one-sample t over the accepted datasets (Boulesteix,
Hable, Lauer and Eugster 2015), each dataset read alone beside it and a confirmed harm on any one a harm. Reserves
are replaced by the pre-fit check's place before the draw and the census. Call 7 is rewritten as the sizing; calls
16 and 21 revised; calls 27 to 32 new; the sizing table, the expected reading and the costs restated (about 47
core-hours through step 2, 70 with step 3, 148 at most, stop 105). The census's second screen passed Toxicity and
Bone marrow transplant: children and failed Occupancy Detection; the pool is 23, Mammographic Mass next for its
second screen.

**2026-11-07, the second screen of the seven single-screened passers.** APSFailure, Internet-Advertisements and
skin-segmentation pass; nomao, clean1 and NHANES_age_prediction fail 1, cylinder-bands fails 5, and they leave the
pool. The screener's five readings of the rule are stated in the dataset rule and applied to every entry (its
calls 11 to 15): a doubt fails only on a kept criterion or the pre-fit check; factors count as bart() expands
them; a share is established only by the suite, the donor, the originating paper or a repository record of the
same file, which fails rice and sources six other shares anew; a suite's value edits are recorded; the loader sets
predictors' types as the donor lists them, in place of the screener's failing of a departure. The pool is 21, 14
with at least 1,000 rows, short of 24, so the census runs on documentation before the seal (dataset rule call
16): its first 80 draws, re-read, give 4 passers on one screen and 3 whose share a second screen may settle. Call
7's options are restated at a pool of 24, with their feasibility at 21; the totals stand.

**2026-11-07, dec-B651's second referee round** (minor revision). Each point, with what changed; none is declined.
- The length guard reads the draw count at the accepting block end only, the other checks at both, so a reference
  can be accepted at 16,000 recorded sweeps (12,000 hold 9,600 draws).
- The fits are sized to the conjunction: 40 a library from the start, each dataset its own blind pilot raised to
  80 past s = 0.95 (call 15); the expected reading restated from the joint chance of enough references and a
  resolving conjunction.
- The percentile bootstrap is dropped, since requiring it beside Welch is Welch alone; each library's skew is read
  and each test's level stated at it (call 26).
- A harm stands only when it also holds against the move's own half, and a harmless reading is flagged against the
  shipped half (call 22); a confirmation reruns at the count the harm was read at, the confirmed family rate
  (about 0.0006) is stated, and a confirmed harm on R2 recommends step 3 (call 25).
- The ways a slower sampler can pass are stated in the guard; the margin's ends equivalent is marked the worst
  case; the floor of 3 is marked the study's own.
- The one-value estimands are named: the guard scores none, so the k study's one-value limits gate nothing;
  criterion 2 reads "every scored quantity"; agreement's chance failure at R2's sizes is stated, and the yield's
  basis over both links.
- Step 3's trigger bound is one-sided 0.025 by a bootstrap over replications, and its ties are broken across
  features.
- E3 and the guard's timing take the fixed time out and time prediction as scored.
- In the dataset rule, criterion 3 keeps the earliest entry of a source that passes every other criterion (spect,
  rice), and the screen file records Pumpkin_Seeds's page as read.
- Compute with the escalations and the overhead counted: R2 about 5.8 core-hours, 20.2 at most; R3's references
  about 8.1; totals about 43 through step 2, 52 with step 3, 107 at most, stop 78.
- R2's and R3's sizes are laid out in call 7 as options with costs, the planner's recommendation marked.

**2026-11-07, dec-B651 revised on its referee's report** (major revision, core design sound). Each point by the
report's section, what changed; none is declined, and where a fix differs from the report's the call says so.
- 1.1 Taken. The k study's length guard (checks at two consecutive block ends of 4,000 sweeps, the first job
  16,000 so the draw count stays at posteriordb's 10,000) and its escalation of the shipped half to 100,000 where
  effective draws are the only failing checks are restored.
- 1.2 Taken. The acceptance is the k study's as written (call 13); posteriordb's 10,000 draws kept and its lag-1
  reported; k's R-hat reported, not required; the expected yield stated, about 2.7 of 6 accepted.
- 1.3 Taken. The tolerance is stated in the guard's units (1 to 2), and D is also read against each half, a
  harmless reading not holding against both flagged (call 22).
- 1.4 Taken. At least 3 accepted references, below which the guard is unresolved and blocks in rule 2 (call 21).
- 2.1 Taken. The margin's use and the per-row Wasserstein reading are marked the study's own; the average
  probability and the band's reference mass are reported (call 24).
- 2.2 Taken in another form. r is the ratio of seconds per sweep, each fit less a one-sweep fit, in step 2 and the
  guard, the fixed share reported (call 20); a percentile bootstrap bound stands beside Welch's.
- 2.3 Taken. The margin is justified as a precision target, with its three coverage equivalents corrected (0.898,
  0.8735 and 0.848; 0.02 is 12, 0.6 and 0.18 units); per-estimand margins are call 14's alternative.
- 2.4 Taken. The aggregate power is the conjunction's (0.24 over 3 datasets and 0.06 over 6 at s = 0.94 and 20
  fits); an internal pilot written now raises every dataset to 40 fits where s exceeds 0.94 (0.96 and 0.92); the
  0.39 bound is for accepted datasets like these.
- 2.5 Taken. The error ratio is reported with its interval on every estimand and class, a lower bound past 1.25
  flagged; the ways a slower sampler can pass are stated.
- 2.6 Taken, and extended to step 2's regimes. A harm shown is confirmed on fresh seeds before it frames anything
  (call 23).
- 3 Taken. The held-out log and Brier scores are reported on the prediction rows after the guard's reading is
  committed, framing nothing; the argument from repeated splits is withdrawn.
- 4 Taken, in [k-move-pricing-datasets.md](k-move-pricing-datasets.md): the leak sentence corrected; the kept
  criteria named population criteria; missing predictors kept, as bart()'s default does; one license standard;
  suite-curated row edits allowed and recorded; Parkinsons's listing 195; spect's two files concatenated; no file
  of the data committed or shipped; the per-entry screen committed as
  [k-move-pricing-screen.csv](k-move-pricing-screen.csv); the seven single-screened passers named for a second
  screen before the seal.
- 5 Taken. The whole of step 3 is committed before any R2 reading is opened and step 2's report waits for it; the
  trigger's selection is stated mechanical, its ties broken.
- 6 Taken, in the dataset rule: the reserve rule, the row split and the score's folds in code with their seeds, the
  full census checksum, R2's and R3's sizes in call 7 to settle before the seal, and the pre-fit check before any
  score.
- 7 Taken. Sizing restated: R2 about 4.9 core-hours, 13.4 at most; R3's references about 6.3, 17.6 at most; totals
  about 42 through step 2, 49 with step 3, 98 at most, stop 73. Draws compacted after scoring.

**2026-11-07, dec-B651: the real-data guard and the dataset rule, for their referee.**
- Step 2's held-out log score guard (6 datasets, 20 splits, 0.01 nats) is replaced by
  [The real-data guard](#the-real-data-guard): each sampler's error at equal time against a reference of the
  posterior, posteriordb's form; references of 16 chains, 8 a library, from prior draws, accepted on posteriordb's
  draw count and R-hat, bulk and tail effective draws of at least 1,000 and dec-B462's agreement; the latent fit's
  mean, standard deviation and interval ends and the probability's mean at fitted and prediction rows, k beside;
  error in units of the error of 100 independent draws (Margossian et al. 2024), harmless below 1 on every row
  estimand and dataset; 20 fits a library; power stated in the unmeasured spread s.
- The candidate list of 26 UCI sets and its blocks of three are replaced by
  [k-move-pricing-datasets.md](k-move-pricing-datasets.md): the survey's published frame, a screen of six criteria
  (outcome, real data, provenance, size, license, feasibility) with the predictive criteria dropped and each
  dropping argued, a re-screen of the 137 documented entries (26 pass in stages 1 to 3, 17 with at least 1,000
  rows; the census does not run), licenses read where the second auditor had not, R2 drawn from every size, and
  the draw written out in R with its seeds. The orchestrator's count of 51 and 27 is corrected there.
- Step 3's form is committed before any R2 reading, its references take the guard's acceptance, and R2 joins its
  test, which then reads 18 datasets.
- Sizing: R2 about 3.5 core-hours, 5.4 at most, in place of 1; totals about 40 expected through step 2, 46 with
  step 3, 75 at most, stop 69.
- New: the Summary for the maintainer, call 7's text, calls 12 to 19 with their alternatives, a Verification item,
  Left out and Not verified lines, the agent line's dataset referee.

**2026-11-07, dec-B633 to dec-B635 and the settled calls.**
- The convergence measure is rebuilt from published comparisons' readings (dec-B633). Dropped: the truth's rank
  among 19 chains at each sweep, its distance from 0.90, the 0.02 bar for when that distance had settled, and the
  sweeps and seconds defined by it, with the rank histogram. In their place, on the same streams and starts:
  Reading 1, He and Hahn's coverage, length, error and running time of a user's fit at each length of a series
  from 25 to 8,000 sweeps with nothing discarded, against sweeps and seconds; Reading 2, rank-normalized split
  R-hat trajectories with bulk and tail effective draws and effective draws per second, and the first length from
  which R-hat stays below Vehtari and coauthors' 1.01; Reading 3, Margossian and coauthors' error against the
  reference on the stored datasets with accepted long runs, the reference's absence where the shipped sampler is
  slowest stated. The streams take a twentieth chain so that they read as five fits of bart()'s 4 chains (two of
  10 for effects), 5 percent more compute (8.2 core-hours expected, was 7.8; 23.3 at most, was 22.1); the plan's
  totals stay about 38 and 44 expected, about 71 at most (was 70), stop 66. Bands are read from 250 sweeps.
  Streams run on to 8,000 on k's R-hat at 2,000, and the second look goes to 100, He and Hahn's count, where a
  standard error exceeds 1.5 times planned. The report lists the figures and tables each reading gives, with
  their published models; Verification, Left out, Not verified, the agent line and the sizing table follow.
- Call 11 ruled (a), dec-B633. Call 8 ruled (a), dec-B634, its reading put in step 2 and The report against R-hat
  below 1.01 and 400 effective draws. Call 6 settled (a) by dec-B494. Call 9 marked the orchestrator's internal
  call, (a). Call 7 stays open, noted with dec-B635.
- Liu and Wu (1999), Liu and Sabatti (2000) and Choi and Hobert (2013), read in full, replace the earlier account
  of them from other papers, in the survey and in Context's description of the move.

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
