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

Two more grids share the same record/compare/print grammar, each behind its
own opt-in flag (or `BENCH_BIGGRID=1` / `BENCH_CALLBACK=1`) and default
baseline name, and leave the grid above untouched: `biggrid` times n in
{1e4, 1e5, 1e6} x numTrees in {75, 200} (not for routine/CI use - the full
grid at full reps can run upwards of an hour; `biggrid quick` restricts it
to the smallest cell as a smoke test), and `callback` times per-draw
callback overhead (none/noop/running-mean at two sizes), recorded at
benchmarks/baselines/bench-sampler-callback-1456e999.csv.

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
(.github/workflows/equivalence.yaml) runs `compare` in this statistical
mode against the current baseline. Bitwise exactness is not local-only:
cpp-tests.yaml runs all three compares per push on the reference build,
pinned to macos-latest arm64 because the baselines are recorded there, so
off that architecture bitwise remains a same-machine check.

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
docs/plans/classic-compare.md - 22 scenarios agree at the null's own rate,
and the four that do not are the zero-weight fit, the two crossvalidation
rows and the unequal-cut-point probe of the change move.

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

Re-evaluates the binary (probit) end-node hyperprior default, chi(1.5, 2),
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
bcf-exact[-weak, -restricted], bcf-latent-exact, categorical-exact,
hazard-exact, heteroscedastic-exact, hurdle-exact, linear-exact,
multinomial-exact, negbin-exact, ordinal-exact, t-exact, logistic-reference,
monotone-reference) and the detailed-balance gates (bd-balance = birth/death,
change-balance = change, swap-balance = swap, perturb-balance = perturb,
rule-gibbs-balance = rule_gibbs) each drive a long fixed-seed MCMC run and
compare the engine's draws to an analytic or brute-force-enumerated target with
a z-score / tolerance bound computed IN-SCRIPT (no recorded baseline), then
quit(status=1L) on deviation. Because the target is derived rather than a
recorded draw, they are deterministic regression detectors that are
host-portable (unlike the equivalence bitwise check). Two further gates,
hazard-reduction and hurdle-reduction, compare draws bitwise against a
hand-built reference fit instead of an analytic target, take no `quick`
argument, and otherwise run and exit like the rest.

    Rscript benchmarks/R/change-balance.R        # full
    Rscript benchmarks/R/change-balance.R quick  # fast smoke

.github/workflows/exact-gates.yaml runs all of them in `quick` mode on push /
pull_request (one install, looped, one ::error:: per failing gate); dispatch it
with mode=full for the long grid. It also carries the two `--cross-host`
compares, outside that loop: the baselines are recorded at full settings and
the settings guard refuses to compare them against a `quick` run. Contrast the
STATISTICAL gates (sbc.R, equivalence.R z-mode), which can false-alarm at the
nominal level and stay schedule / workflow_dispatch only.

## R/constant-*.R - the engine-constants sweeps (measurement, not gates)

One script per fixed engine constant (docs/design/engine-constants.md), named
after the constant it measures: constant-categorical-cap,
constant-linear-leaf-covariates, constant-perturb-width,
constant-testfit-parallel-cutoff, constant-predict-parallel-cutoff,
constant-sparse-density-threshold, constant-gp-max-leaf-size,
constant-person-period-rows and constant-xint-caps. Each prints a table, takes
`quick` for a smaller grid, and has no baseline and no pass/fail exit.

    Rscript benchmarks/R/constant-predict-parallel-cutoff.R quick

Run the timing ones on a quiet machine, one at a time. Two constants are
compile-time with no knob: constant-perturb-width measures the build it is run
against and takes `width=<w>` only as a table label, so a width arm means a
scratch copy of the tree with the constant edited, installed into a private
library - never the working tree. The xint sweep is a probe of what each
channel does at its cap, not a timing run.

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
