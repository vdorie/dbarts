# Benchmarks and equivalence harness

Benchmark and equivalence tooling for the bartcore engine (background:
docs/design/core-generalization.md). Nothing here ships with the package
(excluded via .Rbuildignore). All of it runs against the *installed* dbarts,
except the kernel benchmark, which links the in-tree static library; build
first with `R CMD INSTALL .` from the package root.

## kernels/ - C kernel microbenchmarks

Times the misc kernel vocabulary (docs/design/kernel-vocabulary.md) across
sizes and instruction sets. This is the per-operation floor the generalized
core must hit.

    cd benchmarks/kernels && make && ./bench > baseline-kernels.csv

Notes: `partitionIndices` rows include a per-rep index-array memcpy; subtract
the `memcpyBaseline` rows. Instruction-set sweep covers scalar C through the
host's maximum.

## R/bench-sampler.R - end-to-end timing (the zero-regression gate)

Times sampler throughput (continuous/binary, several sizes), the
embedded-Gibbs pattern (setOffset + run(0, 1)), and single-column
setPredictor updates (accept and reject paths timed separately with
deterministic workloads, since the accept rate of random replacements
depends on the chain state).

    Rscript benchmarks/R/bench-sampler.R record baseline-sampler.csv
    # ... install candidate build ...
    Rscript benchmarks/R/bench-sampler.R compare baseline-sampler.csv

`compare` exits nonzero if any metric is more than 5% slower. Record and
compare on the same quiet machine; append `quick` only for smoke tests.
Sub-millisecond metrics drift a few percent between invocations on a
laptop; confirm a marginal flag by re-running before chasing it.

Three more grids share the same record/compare/print grammar, each behind its
own opt-in flag (or `BENCH_BIGGRID=1` / `BENCH_CALLBACK=1` /
`BENCH_RESTORE=1`) and default
baseline name, and leave the grid above untouched: `biggrid` times n in
{1e4, 1e5, 1e6} x numTrees in {75, 200} (not for routine/CI use - the full
grid at full reps can run upwards of an hour; `biggrid quick` restricts it
to the smallest cell as a smoke test), `callback` times per-draw
callback overhead (none/noop/running-mean at two sizes), recorded at
benchmarks/baselines/bench-sampler-callback-1456e999.csv, and `restore`
times `setState` of a stored state, without force and with it, beside a
sweep at n = 1000 and 1e5 with 200 trees; it has no recorded baseline and
gates nothing.

    Rscript benchmarks/R/bench-sampler.R biggrid record biggrid.csv
    Rscript benchmarks/R/bench-sampler.R callback compare \
      benchmarks/baselines/bench-sampler-callback-1456e999.csv

## R/equivalence.R - statistical equivalence

Verifies two builds target the same posterior. Data are fixed per scenario;
only the MCMC seed varies (20 seeds), so per-summary Welch z-statistics
between builds should look standard normal. |z| > 4 on any summary fails.
Identical output (same RNG stream) is reported as an exact match, so this
also detects unintended RNG shifts from refactors of the current engine.

    Rscript benchmarks/R/equivalence.R record baseline-equivalence.rds
    # ... install candidate build ...
    Rscript benchmarks/R/equivalence.R compare baseline-equivalence.rds

Without a flag, `compare` exits 0 on a moved RNG stream whose posterior
agrees, so a clean exit is not a bitwise pass. `--bitwise` makes it one: any
compared scenario whose draws are not identical fails, whatever its |z|, and
so does a baseline scenario the run did not produce.

Baselines are RNG- and build-dependent, so a *bitwise* compare is a
same-machine check. A baseline is still reusable off-host through the sibling
harnesses' `--cross-host` flag (below), which exempts the channels that cannot
reproduce across machines and gates the rest under a two-tier verdict.
Comparisons must use the same settings (seeds, iterations) as the recorded
baseline; the script enforces this.

`compare` prints a coverage line (scenarios compared / skipped) and warns -
or, with `--strict-coverage`, fails - when the installed engine offers
scenarios the baseline predates. baselines/MANIFEST records each baseline's
role (current, historical, or historical-classic), recording commit,
machine, and scenario list. The scheduled workflow
(.github/workflows/equivalence.yaml) runs `compare` in the statistical
mode against the current baseline. Bitwise exactness is not local-only:
cpp-tests.yaml runs all three compares per push on the reference build,
each with `--bitwise`, pinned to macos-latest arm64 because the baselines
are recorded there, so off that architecture bitwise remains a same-machine
check.

R/bcf-equivalence.R and R/multinomial-equivalence.R are sibling harnesses
for the two multi-forest samplers, with their own current baselines named
by the MANIFEST (which is authoritative; hashes rotate at every
re-record); the three together are the "equivalence trio" the plan docs
gate on. The two sibling harnesses source their draws-axis reductions and
cross-host verdict logic from R/equivalence-common.R; each keeps its own
channel taxonomy and scenario list.

Both take `--cross-host` to compare a baseline recorded on another machine.
The point-in-time snapshot channels (a forest's raw fit/amplitude/varcount
query, a transaction's accept/reject verdict) are exempt: the scenario data
run through the platform libm before any engine code does, so they can never
reproduce bitwise off-host. Every draws-axis channel stays gated. Tier 1 is
the gate - continuous channels within a tight relative bound, integer split
counts exactly - and is a stream-lock detector, not a posterior test. Tier 2
runs only when tier 1 fails, and only adjudicates: its draws axis is one
autocorrelated chain rather than the independent seeds equivalence.R reduces
over, so even with an ESS-adjusted denominator its |z| = 4 bar tolerates a
per-cell shift of over a posterior sd. A tier-2 pass says the failure is not
gross, never that the two builds agree.

Both also take `--bitwise`, same-host compare only and refuses to combine
with `--cross-host`. Without it, a `statChannels` (drawn) mismatch falls
through to the same Welch-z adjudication as equivalence.R's default mode,
which a small RNG-stream move can still pass; `--bitwise` removes that
fallback, so any such mismatch fails outright, whatever its |z|. Point-in-time
snapshot channels already gate their mismatch directly with no fallback in
either mode, and a baseline scenario the run did not produce already fails
either way, so `--bitwise` only closes the drawn channels' statistical escape
hatch - the same effect equivalence.R's `--bitwise` has (c8ff7726). This is
the mode cpp-tests.yaml uses for both siblings' per-push compares.

## R/forests-blocks-equivalence.R - forests that share one total against blocks()

Two spellings of one model, live on one build: `y ~ forest(x1) + forest(x2) +
forest(x3)`, forests with no basis sharing one prior total, against one forest
with `blocks()` at the same trees per group. Six scenarios (three forests
under gaussian, and with k drawn under probit and under gaussian, two forests
of 38 and 37 trees, a stated count beside two that split the rest under gaussian
and under probit), compared by Welch z over 40 seeds a spelling (80 in full) on
the posterior mean of the fit at each row and of sigma or k; a z past 4 fails, a
lone one under 5 after a second set of seeds. In exact-gates.yaml's list, as a
seeded statistical check. Quick mode takes about 25 seconds, full about 70.

    Rscript benchmarks/R/forests-blocks-equivalence.R quick

## R/classic-compare.R - the 0.9-34 comparison (measurement, not a gate)

Widens the cross-release evidence past equivalence.R's nine-scenario
classic record: 26 scenarios written in the 0.9-x vocabulary, so ONE script
runs under an installed dbarts 0.9-34 and under this one. It records rather
than compares live - 0.9-34 is a released package, not a build of this tree,
so neither side can host the other - and `compare` then takes two recordings.

    R_LIBS=<lib-0.9-34> Rscript benchmarks/R/classic-compare.R record old.rds
    R_LIBS=<lib-1.0-0>  Rscript benchmarks/R/classic-compare.R record new.rds mixture=classic
    Rscript benchmarks/R/classic-compare.R compare old.rds new.rds

Every prior and control setting whose DEFAULT moved is pinned on both sides;
`mixture=classic` additionally pins 0.9-x's tree-move mixture, so the compare
isolates the engine, and `mixture=default` leaves 1.0-0 at its own new one,
which measures the kernel change instead. Same Welch-z and disjoint-range
verdict as equivalence.R, at 20 seeds and 1000 draws after 500 burn-in.
`CLASSIC_COMPARE_SCENARIOS` restricts the run, `CLASSIC_COMPARE_SEED_OFFSET`
shifts the seed block (for re-checking a marginal flag), `CLASSIC_COMPARE_CORES`
sets the worker count, and `merge` joins chunked recordings. The 0.9-34 side
is kept as baselines/classic-compare-0.9-34.rds; the 1.0-0 side is the compare
side and is not. A full recording takes under two minutes on seven cores;
re-record the 0.9-34 side only when a scenario is added. Findings:
docs/plans/classic-compare.md - 21 scenarios agree at the null's own rate,
and the five that do not are the zero-weight fit, the two crossvalidation
rows, the unequal-cut-point probe of the change move and, since the
quantile grid was spread over the whole column, the quantile row.

## R/classic-timing.R - 0.9-34 wall time (measurement, not a gate)

The timing half of the same comparison: classic-compare.R records posteriors
and times nothing, and bench-sampler.R times this tree against itself, so
neither answers "how much faster than the released package". This one fits
one design cell through the BayesTree-style door under whichever dbarts
`R_LIBS` points at and prints the fit's wall-clock seconds. It cannot run in
CI - it needs two INSTALLED releases - and it needs an idle machine.

Four cells, Friedman data at ten predictors, 200 trees, 1000 draws after
500: `a` n = 1000 one chain, `b` n = 10000 one chain, `c` n = 1000 four
chains on four threads, `d` n = 1000 probit. A fresh process per fit, so
only one dbarts ever loads; alternate which library goes first on each
repetition so drift over the run falls on both alike:

    for cell in a b c d; do
      for rep in 1 2 3 4 5; do
        if [ $((rep % 2)) -eq 1 ]; then libs="$LIB_NEW $LIB_OLD";
                                   else libs="$LIB_OLD $LIB_NEW"; fi
        for lib in $libs; do
          R_LIBS="$lib" Rscript benchmarks/R/classic-timing.R "$cell"
        done
      done
    done

Take the median over the repetitions per (cell, library) and report the
ratio. Recorded run and numbers: docs/plans/classic-compare.md, "Wall time"
- 1.7x at n = 1000, 2.6x at n = 10000, 1.7x on four chains and 1.6x on a
probit fit, on a four-core x86-64 box. Re-time on the x86 leg below rather
than on the arm64 development machine, and check `/proc/loadavg` first.

## R/binary-hyperprior.R - the binary k prior study (measurement, not a gate)

Re-evaluates the binary (probit) leaf hyperprior default, chi(1.5, 2),
over a grid of priors and a case set far wider than the study that set the
default: twenty-four chi(df, scale) arms (df in {1, 1.25, 1.5, 2, 3} crossed
with scale in {1, 2, 5, Inf}, plus df in {1.5, 3} crossed with scale in
{0.5, 0.25}) and four fixed-k arms (k in {1, 1.5, 2, 3}), over 162 simulated
cells (six data-generating processes by three sample sizes by three predictor
counts by three base rates) and twenty-two real datasets scored by repeated
80/20 splits - six from R and its recommended packages, sixteen from the UCI
Machine Learning Repository, fetched on demand and cached, by
benchmarks/R/uci-binary.R. Arms are paired inside a case and repetition (same
data, same seed). Scores, all on held-out rows: log score and Brier against
the outcome; against the known truth on the simulated cells, the coverage and
width of the 90 percent interval for the true probability and the RMSE of the
posterior mean probability; plus the sampled k and the fit time.

    Rscript benchmarks/R/binary-hyperprior.R blocks           # list the blocks
    Rscript benchmarks/R/binary-hyperprior.R sim:weak:500 DIR # one block
    Rscript benchmarks/R/binary-hyperprior.R real:biopsy DIR
    Rscript benchmarks/R/binary-hyperprior.R all DIR          # every block
    Rscript benchmarks/R/binary-hyperprior.R summarize DIR    # the tables

One block is one invocation and writes one rds into DIR, so a full run splits
across a session; a large simulated block narrows further by appending a
predictor count (`sim:weak:2000:50`). A full run is 40 blocks (18 simulated,
22 real), 73,248 fits at the default chain length; the real blocks vary
widely in cost, since the UCI datasets run up to 48,842 rows (capped
at 4,000 training rows per split past 5,000 rows), so time it on the machine
you'll run it on rather than assuming a figure here. `BINARY_HYPERPRIOR_CORES`,
`_REPS` and `_SPLITS` set the worker count, the simulated repetitions and the
real-data splits, and `_BURN`, `_DRAWS` and `_CHAINS` the MCMC length a fit
gets (the plan doc's three chain lengths are 500/500/1 "short", 2000/2000/4
"long" and 8000/8000/8 "probe"); `quick` is a smoke run and its files are
marked so summarize will not mix them with a real one. No baseline and no
pass/fail exit: the verdict is written by a person. Findings:
docs/plans/binary-hyperprior.md - keep chi(1.5, 2); nothing in the grid clears
the doc's bar for moving the default, and the improper scale and every fixed
k are worse at every chain length tried. Poor mixing explains most of the
coverage shortfall for a sampled k and almost none of it for a fixed k, but a
residual shortfall persists even at the doc's longest chains.

## R/*-exact.R, *-balance.R - deterministic exact-posterior gates

The exact-posterior gates (aft-exact, aft-hetero-pit, backfit-exact,
bcf-exact[-weak, -restricted], bcf-latent-exact, blocks-exact, categorical-exact,
hazard-exact, heteroscedastic-exact, hurdle-exact, linear-exact,
mask-redraw-exact, multinomial-exact, negbin-exact, ordinal-exact,
outer-step-exact, probit-k-scale-exact, split-probs-exact, t-exact, logistic-reference, monotone-reference) and the detailed-balance gates (bd-balance = birth/death,
change-balance = change, swap-balance = swap, perturb-balance = perturb,
rule-gibbs-balance = rule_gibbs) each drive a long fixed-seed MCMC run and
compare the engine's draws to an analytic or brute-force-enumerated target with
a z-score / tolerance bound computed IN-SCRIPT (no recorded baseline), then
quit(status=1L) on deviation. Because the target is derived rather than a
recorded draw, they are deterministic regression detectors that are
host-portable (unlike the equivalence bitwise check). Two further gates,
hazard-reduction and hurdle-reduction, compare draws bitwise against a
hand-built reference fit instead of an analytic target, take no `quick`
argument, and otherwise run and exit like the rest. negbin-mixing is a
fixed-seed mixing gate rather than an exactness one: on two cells whose
posterior on r spreads over several grid values (r0 = 8 at n = 400, r0 = 10 at
n = 500) with the default forest, two chains whose shapes are set apart before
the run (1 and 50, through the stored state: the sampler has no setter) must
leave their starts, agree by split-Rhat, cover the true r and cover fresh
counts at 90% by randomized PIT; a cell with no spread fails as undefined, and
`quick` runs one dataset a cell and full three, with a further check
on the 95% set over the six runs (the header of the script says what each
catches).

mask-redraw-exact is the gate for a mask redrawn every sweep by a larger
sampler: a two-part mixture over ten rows and one tree, whose joint posterior
over trees and masks is enumerated, against a combined sampler that alternates
a sweep with a redraw of every row's membership and `$setActiveRows`. It
compares the long-run membership probabilities and the fitted means, on a
gaussian arm, an unordered-factor arm and a probit arm. It fails if the trees
a sampler may hold depend on the mask, and its probit arm fails if a row
switched back in keeps the latent it left with (docs/design/empty-leaf-veto.md,
"What counts as empty: membership").
Quick mode is about 45 seconds and bounds each of 44 z statistics at 4.5, a
false failure once in 1400 fresh seeds for a correct sampler; an arm can be
run alone by naming it.

    Rscript benchmarks/R/mask-redraw-exact.R quick probit

outer-step-exact is the gate for where a sweep draws its augmentation variables
against an outer step that reads the fit: 32 rows in four cells of one ordinal
predictor and one tree, the tree space enumerated, an offset gamma * v with
gamma on the grid -1, -0.5, ..., 4 under a uniform prior, redrawn after every
sweep and installed with `$setOffset`, and k drawn under its chi(1.5, 2) prior
(the probit rescaling step runs), the leaves and k integrated by quadrature.
Four arms: probit with the outer step integrating z out (pnorm of the fit) and
conditioning on z, and logistic integrating omega out (plogis) and conditioning
on omega. It compares P(gamma = g) at each grid point holding at least 1e-4 of
the exact mass (eight of eleven) and the fit at each cell, at |z| 4.5 on 12
statistics an arm. An outer step that conditions on the latent is exact
whatever the placement of the draw, as long as the latent it reads is a draw
given the fit; one that integrates the latent out is exact only if the latent
is drawn before the next sweep's trees read it. The gate checks where the draw
sits and nothing else: a stale working pointer, a missed cache invalidation and
a latent split from the parameter drawn with it leave its arms exact or out of
its reach. Quick mode is 4 chains of 250,000 sweeps and about 2 minutes (35,
23, 35 and 28 seconds an arm); full mode 4 chains of 2,000,000. An arm can be
run alone by naming it.

Largest |z| in quick mode, P(gamma) over the fit, on the shipped build of main
(the draw after the trees), a build with the draw first, and a build with the
draw after the rescaling step:

| arm | draw last (main) | draw first | draw after rescaling |
|---|---|---|---|
| probit-integrated | 119.9 / 22.2 | 2.2 / 2.0 | 109.8 / 16.1 |
| probit-conditional | 1.4 / 1.5 | 1.6 / 1.9 | 1.7 / 0.8 |
| logistic-integrated | 17.8 / 8.7 | 1.4 / 1.6 | 1.4 / 1.6 |
| logistic-conditional | 18.2 / 3.8 | 1.0 / 0.9 | 1.0 / 0.9 |

Main fails three arms, the two that integrate the latent out and the logistic
one that conditions on omega (where setOffset leaves the logistic working
response at the old offset). The draw-first build passes all four; the draw
after the rescaling step fails the probit integrated arm alone, the one arm
that sees the rescaling step (logistic has none).

    Rscript benchmarks/R/outer-step-exact.R quick probit-integrated

blocks-exact is the gate for the `blocks()` prior, each group an equal part of
its forest's prior variance whatever its tree count: one forest of three trees
in two groups of two and one, each tree a root or its one split, enumerated
over its 8 structures with the leaves integrated in closed form, at k = 2 and
under k ~ chi(1.5, 2) by an outer quadrature. It compares the split
probabilities, the fit at each cell and k at |z| 4.5; the per-tree prior,
every leaf at the forest's leaf scale, misses the fits by over 300 standard
errors. Quick mode is about 10 seconds.

    Rscript benchmarks/R/blocks-exact.R quick

split-probs-exact is the gate for a zero split probability: one tree over two
0/1 columns, each used up below a split on it, and a third of probability
zero that carries most of the signal, whose posterior over the nine trees left
is enumerated with no growth where no column of positive probability is left.
It compares the split-count classes and the fitted means, and fails on any
split on the third column. The `vars` arm leaves the third column out by a
column list instead and is held to the same enumeration. A sampler that
proposes the first available column where none of positive probability is
left misses by hundreds of standard errors. Quick mode is about a second an
arm and bounds each of 10 z statistics at 4.5, a false failure once in 5400
fresh seeds for a correct sampler.

    Rscript benchmarks/R/split-probs-exact.R quick vars

probit-k-scale-exact is the gate for the probit rescaling step: one frozen
tree under probit with k ~ chi(1.5, 2), whose posterior of k and leaf means are
exact by one quadrature a leaf, on four arms (every leaf separated, none, an
offset, a mask), and three arms of two frozen single-tree forests with no basis
sharing one drawn k, the second a stump integrated in an outer quadrature: alone,
beside a held factor stump the step reads as offset, and with a user offset;
and two where a sharer leaves by `$setForest(f, sd = )` before the run, its stump
then read as offset beside those left at their creation parts. It bounds
P(k < q) at the exact deciles and the unseparated leaves' means at |z| 4.5, and the separated arm's batch-means error of P(k < q) at
0.01 in full mode (0.014 quick), which a sampler without the step does not
meet, so a step that silently declines fails. It cannot see k left unscaled
after a correct draw; the engine's mapping test catches that. Quick mode is
about 2.5 minutes, full about 7.

    Rscript benchmarks/R/probit-k-scale-exact.R quick

bd-balance takes a second arm, `zeroweight`, which the workflow runs after its
loop: the same enumeration under a weight vector that zeroes two adjacent
cells, installed on whatever tree the chain has grown, against the exact
posterior in which a leaf of only zero-weight rows is legal and scores nothing.
A third, `zerocreate`, is the same gate with the weights given at creation and
the zeroed cells' response 3 higher, so the exact posterior calibrates against
the range of the rows in the likelihood, which a sampler that read every row
misses.

    Rscript benchmarks/R/bd-balance.R quick zerocreate

negbin-exact takes a third arm, run inside the script: a fixed shape with 4 rows
in of 25 a cell and the rows out 30 higher, left out by 0/1 weights at creation.
The reference takes its likelihood and its center over the rows in, and the gate
also requires the fit to miss the reference whose center is over every row.

    Rscript benchmarks/R/change-balance.R        # full
    Rscript benchmarks/R/change-balance.R quick  # fast smoke

frequency-exact is the gate for frequency weights (docs/plans/frequency-weights.md,
S7). Seven arms hold a weighted fit to the fit on its rows written out, draw for
draw over 300 sweeps within 1e-12: gaussian, a linear leaf, two forests,
useQuantiles with fewer cuts than unique values, a variance forest, nbinom at a
fixed and a drawn shape, and multinomial; an eighth holds cbind(k, f) at a row
weight to its two adjacent logistic rows of counts within 1e-10. Three are in
distribution, at |z| 4.5 by batch means: sigma on one leaf at real weights
against its exact posterior, Student-t's nu on fixed residuals at whole, real
and zero counts against the grid's exact probabilities, and Student-t against
its expanded rows by Welch z on the fits, sigma and nu. Quick mode takes about
4 seconds, full about 20; each arm takes its name.

    Rscript benchmarks/R/frequency-exact.R quick sigma nu

.github/workflows/exact-gates.yaml runs all of them in `quick` mode on push /
pull_request (one install, looped, one ::error:: per failing gate); dispatch it
with mode=full for the long grid. exact-gates-weekly.yaml runs the full mode of
every gate here weekly on a longer limit, reading the main job's gate list from
exact-gates.yaml and running the monotone enumeration gate in full beside it
(bcf-latent-exact stays quick). exact-gates.yaml also carries the two `--cross-host`
compares, outside that loop: the baselines are recorded at full settings and
the settings guard refuses to compare them against a `quick` run. Contrast the
STATISTICAL gates (sbc.R, equivalence.R z-mode), which can false-alarm at the
nominal level and stay schedule / workflow_dispatch only.

## R/sbc.R - simulation-based calibration

sbc.R draws theta0 from the sampler's own prior, simulates data, refits and
ranks theta0 among the posterior draws; ranks are uniform for a calibrated
sampler. .github/workflows/sbc.yaml runs a matrix of arms, M = 84 functionals,
each band at alpha = 0.05 / 84, so a correct matrix is green about 95% of the
time. It is a statistical test: it runs on a schedule, by hand, and on a push
that changes its own file, not on other pushes.

The probit-k arm draws k under chi(1.5, 2) and ranks it beside the other
functionals. k is flagged when its ecdf difference leaves the band or its
chi-square p (20 bins over the 100 ranks 0 to 99) is below the band's alpha,
6e-4: k is tested by two rules, while M stays 84. Detection floor, at
R = 600, L = 99, thin 700: an error of one half in the shape of the k^2
conditional (a missing Jacobian, a leaf too many, the prior's df off by one) is
flagged about 80% of the time; an error that scales with the number of leaves
(every k draw 5% high) is flagged at R = 50. Smaller errors pass.

With the probit rescaling step k mixes fast enough that at thin 700, after a
10,000-sweep burn, every dataset of the arm keeps about 50 or more effective draws of
the 99 kept (fewest 49 to 55 on the mixing census by two estimates), so the arm
is read as any other: k's end bins near 30 and its chi-square p meaningful.
Before the step a correct sampler flagged this arm about one run in three
(docs/design/probit-k-calibration.md).

By hand: `schedule` fires from main, the default branch, and a dispatch
runs from any branch the file is on. A push to main that changes sbc.yaml
already starts the matrix; otherwise run

    gh workflow run sbc.yaml

and read each job's report in the log of its "Run the SBC calibration check"
step (`gh run view <id> --log`, or the job's page). Green is the band verdict,
not the reading: for every arm look at the chisqP column. Of the 84 functionals expect about four
under 0.05 and one under 0.01 (they are not independent, so read nine or more
under 0.05, or any under 6e-4, as a finding to look at); discrete-selfcheck
must be green. The probit-k job takes
67 minutes on the laptop and up to twice that on the runner.

The frequency arm is run by dispatch, outside the matrix and its M: a gaussian
fit at real weights, Uniform(0.5, 2), whose theta0 draws sigma from the prior
tilted by the frequency likelihood's factor (sigma^2)^(-(W - n) / 2), W the
weights' total, and y0 at sd sigma / sqrt(w) (sbcFrequencySigmaDraw), so its
ranks are uniform exactly when the sampler reads the weights as counts. Each
functional is read at its own 5% band. At R = 200, L = 200, thin 30 it takes
about 11 minutes on the laptop.

    Rscript benchmarks/R/sbc.R frequency 200 200 30

## R/probit-k-mixing.R - k's mixing under probit (measurement, not a gate)

How fast k mixes under probit with the rescaling step: `census` gives k's and the average fit's integrated autocorrelation time over the SBC
probit-k arm's prior-drawn datasets (the arm's start, 30,000 sweeps of burn-in
and 200,000 recorded, about 25 seconds a dataset), `truth` the start-from-truth
drift of log k at R replications (2000 take about 5 minutes), `summarize` the
tables. Findings, the step-off figures included: docs/design/probit-k-scale-move.md.

    Rscript benchmarks/R/probit-k-mixing.R census 1 50 out-1.rds
    Rscript benchmarks/R/probit-k-mixing.R truth 19 2000 truth.rds
    Rscript benchmarks/R/probit-k-mixing.R summarize out-1.rds truth.rds

## R/constant-*.R - the engine-constants sweeps (measurement, not gates)

One script per fixed engine constant (docs/design/engine-constants.md), named
after the constant it measures: constant-categorical-cap,
constant-linear-leaf-covariates, constant-perturb-width,
constant-testfit-parallel-cutoff, constant-predict-parallel-cutoff,
constant-sparse-density-threshold, constant-gp-max-leaf-size,
constant-person-period-rows and constant-xint-caps. Each prints a table, takes
`quick` for a smaller grid, and has no baseline and no pass/fail exit.

    Rscript benchmarks/R/constant-predict-parallel-cutoff.R quick

Run the timing ones on a quiet machine, one at a time. constant-perturb-width
sets the window by `dbartsControl(n.perturb.cuts = )`, sweeping 1, 2 and 4 or
the one `width=<w>` names. The xint sweep is a probe of what each channel does
at its cap, not a timing run.

## R/move-census.R - the move census (measurement, not a gate)

Stage 0 of the tree-mixing falsifier (docs/design/tree-mixing-proposals.md
section 6.1): per-move acceptance on both denominators, the log-likelihood
difference among rejected structural proposals, change acceptance by node
depth, and the acceptance a same-variable cut move would have had at a
schedule of cut displacements. It reports; it never fails.

Unlike everything else here it needs a SPECIAL BUILD, because the records
come from scaffolding in the move kernels that compiles only under
-DBARTCORE_MOVE_CENSUS. The script's header carries the install commands and
the record format; the flag rides CPPFLAGS through R_MAKEVARS_USER and the
build goes to a private library, so the ordinary one is untouched.

    Rscript benchmarks/R/move-census.R                 # run, then summarize
    Rscript benchmarks/R/move-census.R summarize DIR   # existing files

## The x86 leg

The development machine is arm64 macOS, which masks Linux and x86 bugs:
a Linux-only build break, a CPUID-misdetected AVX2 kernel and an ABI
mismatch segfault each surfaced first on an x86 box and nowhere on
macOS. An engine landing gets an x86 leg alongside the local battery,
and it is the only place to time bench-sampler on x86 or to exercise a
new SIMD kernel.

The x86 host itself is local to the development environment and is
described nowhere in the tree; the pattern below assumes an ssh-reachable
x86-64 Linux box with R, the compilers and valgrind installed, and a
private library per run. Check its load before any timing run, and
never assume a core count: threading benches saturate at the physical
cores and regress past them.

Working pattern: ship the source (rsync excluding .git, or `git archive`
of the exact sha when the local tree is dirty), install with
`R CMD INSTALL --preclean -l ~/rlib-<tag>` into a private library, set
`R_LIBS=~/rlib-<tag>` on every R invocation, then run tests/cpp (plain
and ASAN), the full tinytest suite, and the equivalence trio.

Expected verdicts there: tinytest FAILURES == 0 (test-simd.R is bitwise
across dispatch levels 0/2/5/7/8 within that host; gate on failures, not
the total count, which differs from arm); on the shipped build every
equivalence scenario reports "max |z|" rather than "identical draws"
against an arm64-recorded baseline and passes in statistical mode
(docs/architecture.md, "Reproducibility contract"); the reference build
is where cross-ISA bitwise reproduction is asked for. Timing runs need
the host idle (`/proc/loadavg` first). Sanitizer binaries may need
`setarch $(uname -m) -R` on hosts with high ASLR entropy.

## tests/cpp - bartcore component tests

C++-level exact tests of the new engine's math against independently coded
references, plus end-to-end smoke runs:

    cd tests/cpp && make run
