# equivalence-ci

agent: sonnet
rng: neutral (test/CI only)
budget: ~150 lines

## Goal

The equivalence gate runs on a schedule in CI in its platform-
independent statistical mode; baseline coverage and provenance are
visible instead of being git archaeology.

## Context

- benchmarks/R/equivalence.R: compare mode is identical()-first, else
  Welch z per summary, FAIL |z| > 4. The bitwise baselines are
  ARM-recorded (TODO note), which is what has kept it out of CI.
- compare iterates names(baseline$results): scenarios added after a
  baseline was recorded are silently skipped - no coverage manifest.
- Three classic-era baseline .rds (19c499a, 5430fdb, fbd2168) can never
  be re-recorded (engine deleted); three bartcore-native ones exist.
  Docs call the classic files "the permanent cross-engine reference".
- benchmarks/baselines/ has no marker for which file is current.

## Constraints

- The CI job asserts the statistical mode only; bitwise exactness stays
  a local, same-machine check.
- Out of scope: changing scenario content; the speed gate
  (bench-sampler stays maintainer-run per its README).

## Steps

1. equivalence.R compare: warn (or fail with --strict-coverage) when
   the installed engine offers scenarios the baseline lacks; print a
   coverage line (scenarios compared / skipped, by name).
2. benchmarks/baselines/MANIFEST: one line per file - role (current |
   historical-classic), recording commit, machine/arch, scenario list.
   Mark the classic-era files historical evidence of the cutover.
3. Scheduled workflow (.github/workflows/equivalence.yaml, cron +
   manual dispatch): install package, run compare in z-mode against the
   current baseline, quick profile; fail on |z| > 4 or coverage
   regression.
4. Update benchmarks/README.md and the "permanent reference" phrasing
   in docs/design/public-surface.md to match (historical, not gate).

## Verification

- Local: Rscript benchmarks/R/equivalence.R compare <current>.rds
  passes with the coverage line printed.
- The workflow runs green via workflow_dispatch before merging.

Runner note (2026-07-06): scenario x seed units parallelized over
forks (bitwise-verified against the serial equivalence-235bebc.rds,
18/18 identical; ~11 min -> ~1.5 min). EQUIVALENCE_CORES overrides;
Windows serial. CI wiring should set it to the runner core count.
