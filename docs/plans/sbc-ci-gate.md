# sbc-ci-gate

agent: sonnet
rng: neutral - no package code (R/, src/) touched; benchmarks/R/sbc.R changes
  are CLI-block only (exit-status wiring), no draw sequence is altered.
budget: ~10 lines in sbc.R's CLI block, ~40 lines of new workflow yaml, this
  plan file.

## Goal

Add benchmarks/R/sbc.R's gaussian calibration baseline as a scheduled,
non-blocking CI workflow: it reports weekly and on manual dispatch, and never
gates a pull request, because a statistical test at R=100 will occasionally
false-alarm at its nominal level.

## Context

sbc.R (docs/design reference: the file's own header) runs simulation-based
calibration against the installed package and prints a per-functional
PASS/FLAG table via sbcReport, but as an Rscript entry point it always exits
0 - there was no way for a CI job to fail on a FLAG. equivalence.yaml is the
existing model for a scheduled-only, non-PR-gating statistical check
(schedule + workflow_dispatch, no push/pull_request), so sbc.yaml is built
directly from it.

## Change

1. benchmarks/R/sbc.R, CLI block (`if (sys.nframe() == 0L)`): capture
   sbcReport's return value (`invisible(verdicts)`, already a
   PASS/FLAG character vector) and, when the env var SBC_FAIL_ON_FLAG is
   non-empty, `quit(status = 1L, save = "no")` if any verdict is "FLAG".
   Unset (the default), behavior is byte-for-byte unchanged: no quit() call,
   Rscript exits 0 as before. source()-based usage (`source("sbc.R"); res <-
   runSbc(...)`) never reaches this block (guarded by sys.nframe() == 0L), so
   it is untouched. The env-var style matches
   benchmarks/R/equivalence.R's existing `Sys.getenv("EQUIVALENCE_CORES",
   "")` / `nzchar()` convention rather than introducing positional-flag
   parsing into a file whose CLI is otherwise purely positional
   (config R L thin).
2. .github/workflows/sbc.yaml, modeled on equivalence.yaml: same
   checkout / setup-r / setup-r-dependencies / install steps,
   `permissions: read-all`, `runs-on: ubuntu-latest`. Triggers are
   `schedule` (weekly, Wednesday 03:00 UTC - distinct from equivalence's
   Monday 06:00) and `workflow_dispatch` only; no push/pull_request trigger,
   which is what keeps this off the PR-gating path while still turning the
   scheduled run red and visible on a FLAG. The job sets
   `SBC_FAIL_ON_FLAG: 'true'` and runs
   `Rscript benchmarks/R/sbc.R gaussian 100 200 30` (the decided
   configuration: gaussian family, R=100, L/thin at sbc.R's own defaults).
   Header comment states the statistical-check / non-gating rationale and
   flags that the BCF glue-on sigma config is under separate investigation
   and is not part of this workflow (gaussian only).

## Seeding

sbc.R's gaussian path is already fully seeded by the time the CLI reaches
runSbc: sbcConfig()'s default `configSeed = 1L` fixes the design (x, xTest,
yBuild) and sbcConfig calls `set.seed(configSeed)` before drawing them;
runSbc()'s default `seed = 20260709L` re-seeds before building the sampler
and before the replication loop. Both are literal constants, not derived
from wall-clock time, PID, or any other run-to-run varying input, and
n.threads is fixed at 1L, so there is no source of cross-run nondeterminism
in the gaussian CLI path. Confirmed empirically (see Verification): two back-
to-back `Rscript benchmarks/R/sbc.R gaussian 20 60 10` runs produced
byte-identical rank tables and verdicts (only a millisecond-scale timing
digit in the printed wall-clock line differed). No additional CI-specific
seeding was added; the existing defaults already give the weekly run
deterministic ranks while remaining fully sensitive to any draw-shifting
code change (a changed engine draw reshuffles the same fixed-seed replay).

## Verification

Private library /Users/vdorie/.claude/jobs/7fe13675/tmp/Rlib-land
(R_LIBS, bartcore build at the worktree's base commit 97b90a8); no
R CMD INSTALL run in this worktree (shared library owned by another job).

- Healthy baseline, the exact command the workflow runs:
  `Rscript benchmarks/R/sbc.R gaussian 100 200 30`. All 7 functionals
  (avg.f, f.star1..5, sigma) PASS. Wall-clock reported by the script:
  26.8s total (0.268s/rep). Measured end-to-end (R startup + package load +
  run): 33.6s. Exit code 0.
- SBC_FAIL_ON_FLAG=true against the same healthy run (R=20, smaller for
  speed): all PASS, exit code 0 - confirms the flag doesn't perturb a clean
  run.
- Failure path: a scratch copy of sbc.R (never committed) with the sigma
  rank line changed to `sum(as.numeric(res$sigma) < sig0 * 4)` to force
  miscalibration. `SBC_FAIL_ON_FLAG=true Rscript <scratch> gaussian 30 60
  10` produced two FLAGs (sigma, f.star1) and exited 1. The same scratch
  script run without SBC_FAIL_ON_FLAG set exited 0 despite the FLAGs -
  confirms the opt-in gate is what changes the exit code, not the FLAG
  itself.
- Determinism: two consecutive `Rscript benchmarks/R/sbc.R gaussian 20 60
  10` runs diffed byte-identical except a one-digit difference in the
  printed per-rep timing average; every rank value and verdict matched.
- air format --check benchmarks/R/sbc.R: clean (exit 0).
- lintr::lint("benchmarks/R/sbc.R"): 0 lints.
- actionlint / yamllint: neither is installed on this machine. Validated
  instead by parsing both equivalence.yaml and the new sbc.yaml with R's
  yaml::read_yaml and diffing their structure: the shared steps
  (checkout, setup-r, setup-r-dependencies, install) match field-for-field;
  sbc.yaml correctly omits equivalence's fork-parallelism step (a single
  serial R=100 run needs no core-count tuning) and adds
  SBC_FAIL_ON_FLAG to env. Both files parse as valid YAML with no syntax
  errors.

## Status

- 2026-07-10: CLI exit-status opt-in, sbc.yaml, and this plan landed on
  wt/sbc-ci-gate. Verification above run against the shared private library
  at the worktree's base commit; no package code touched (R/, src/, man/,
  inst/ untouched, per the item's scope).
- 2026-07-10: LANDED as ed34eea (squash of wt/sbc-ci-gate).
  Reviewer confirmed sbcReport already returns invisible(verdicts)
  (the CLI hook reads a real value) and independently re-ran the
  exact workflow command from a private library: all 7 functionals
  PASS, exit 0.
