# Benchmarks and equivalence harness

Phase-0 tooling for the core-generalization work
(docs/design/core-generalization.md). Nothing here ships with the package
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
Passing `engine=new` runs the timed side on the bartcore engine, so
recording under classic and comparing under `engine=new` measures the
cross-engine gap on one build. Sub-millisecond metrics drift a few percent
between invocations on a laptop; confirm a marginal flag by re-running or
by interleaving the engines by hand before chasing it.

## R/equivalence.R - statistical equivalence

Verifies two builds target the same posterior. Data are fixed per scenario;
only the MCMC seed varies (20 seeds), so per-summary Welch z-statistics
between builds should look standard normal. |z| > 4 on any summary fails.
Identical output (same RNG stream) is reported as an exact match, so this
also detects unintended RNG shifts from refactors of the current engine.

    Rscript benchmarks/R/equivalence.R record baseline-equivalence.rds
    # ... install candidate build ...
    Rscript benchmarks/R/equivalence.R compare baseline-equivalence.rds

Baselines are RNG- and build-dependent; regenerate them rather than reusing
across machines. Comparisons must use the same settings (seeds, iterations)
as the recorded baseline; the script enforces this.

`compare` prints a coverage line (scenarios compared / skipped) and warns -
or, with `--strict-coverage`, fails - when the installed engine offers
scenarios the baseline predates. baselines/MANIFEST records each baseline's
role (current, historical, or historical-classic), recording commit,
machine, and scenario list. The scheduled workflow
(.github/workflows/equivalence.yaml) runs `compare` in this statistical
mode against the current baseline; bitwise exactness stays a local,
same-machine check.

Passing `engine=new` runs the comparison side against the bartcore engine
(src/bartcore/) through a standalone shim (tests/cpp/rshim.cpp, built on
demand by benchmarks/R/bartcore-shim.R); this is the phase-1 correctness
gate. A recorded old-engine baseline lives in baselines/.

## tests/cpp - bartcore component tests

C++-level exact tests of the new engine's math against independently coded
references, plus end-to-end smoke runs:

    cd tests/cpp && make run
