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
An `engine=` flag is accepted and ignored: the installed package always
runs the bartcore engine. Sub-millisecond metrics drift a few percent
between invocations on a laptop; confirm a marginal flag by re-running
before chasing it.

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
mode against the current baseline; bitwise exactness stays a local,
same-machine check.

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

## R/*-exact.R, *-balance.R - deterministic exact-posterior gates

The exact-posterior gates (aft-exact, bcf-exact[-weak, -restricted],
categorical-exact, heteroscedastic-exact, linear-exact, multinomial-exact,
negbin-exact, ordinal-exact, t-exact, logistic-reference, monotone-reference)
and the detailed-balance gates (bd-balance = birth/death, change-balance =
change, swap-balance = swap, perturb-balance = perturb, rule-gibbs-balance =
rule_gibbs) each drive a long fixed-seed MCMC run and
compare the engine's draws to an analytic or brute-force-enumerated target with
a z-score / tolerance bound computed IN-SCRIPT (no recorded baseline), then
quit(status=1L) on deviation. Because the target is derived rather than a
recorded draw, they are deterministic regression detectors that are
host-portable (unlike the equivalence bitwise check).

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
